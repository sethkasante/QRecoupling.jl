
# ========================================================================================
# Analytic Evaluator (Continuous & Complex Regimes)
# Maps CycloResults to generic complex parameters (q ∈ ℂ) or the unit circle q = exp(iθ).
# ========================================================================================

const UNIT_CIRCLE_CACHE  = LRU{BigFloat, Tuple{Vector{BigFloat}, Vector{BigFloat}}}(maxsize=100)
const ANALYTIC_CACHE     = LRU{Complex{BigFloat}, Tuple{Vector{BigFloat}, Vector{BigFloat}}}(maxsize=100)

function get_unit_circle_table(theta::BigFloat, D_max::Int)
    if haskey(UNIT_CIRCLE_CACHE, theta)
        lmag, lphs = UNIT_CIRCLE_CACHE[theta]
        if length(lmag) >= D_max
            return lmag, lphs
        end
    end
    lmag, lphs = build_unit_circle_table(D_max, theta)
    UNIT_CIRCLE_CACHE[theta] = (lmag, lphs)
    return lmag, lphs
end

function get_analytic_table(q::Complex{BigFloat}, D_max::Int)
    if haskey(ANALYTIC_CACHE, q)
        lmag, lphs = ANALYTIC_CACHE[q]
        if length(lmag) >= D_max
            return lmag, lphs
        end
    end
    lmag, lphs = build_analytic_table(D_max, q)
    ANALYTIC_CACHE[q] = (lmag, lphs)
    return lmag, lphs
end

function build_unit_circle_table(D_max::Int, theta::BigFloat)
    V_mag = zeros(BigFloat, D_max)
    V_phs = zeros(BigFloat, D_max) 
    V_valid = trues(D_max)
    
    half_theta = theta / 2
    pi_big = BigFloat(π) 
    
    for n in 1:D_max
        val_sin = 2 * sin(BigFloat(n) * half_theta)
        
        if iszero(val_sin)
            V_valid[n] = false
        else
            V_mag[n] = log(abs(val_sin))
            p_n = BigFloat(n) * half_theta + (pi_big / 2) 
            if val_sin < 0
                p_n += pi_big
            end
            V_phs[n] = p_n
        end
    end
    
    @inbounds for d in 1:D_max
        !V_valid[d] && continue
        for m in (2 * d):d:D_max
            if V_valid[m]
                V_mag[m] -= V_mag[d]
                V_phs[m] -= V_phs[d]
            end
        end
    end
    return V_mag, V_phs
end

function build_analytic_table(D_max::Int, q::Complex{BigFloat})
    V_mag = zeros(BigFloat, D_max)
    V_phs = zeros(BigFloat, D_max)
    V_valid = trues(D_max)
    
    # Rolling iterator avoids O(D log D) complex exponentiation overhead
    q_curr = q
    for n in 1:D_max
        term = q_curr - one(q)
        if iszero(term)
            V_valid[n] = false
        else
            V_mag[n] = log(abs(term))
            V_phs[n] = angle(term)
        end
        q_curr *= q 
    end
    
    @inbounds for d in 1:D_max
        !V_valid[d] && continue
        for m in (2 * d):d:D_max
            if V_valid[m]
                V_mag[m] -= V_mag[d]
                V_phs[m] -= V_phs[d]
            end
        end
    end
    return V_mag, V_phs
end

# ------------------------------------------------------------------------------
# High-Performance Continuous Projections
# ------------------------------------------------------------------------------

@inline function _project_to_complex_analytic(m::CycloMonomial, lmag_table::Vector{BigFloat}, lphs_table::Vector{BigFloat}, 
                                              z_lmag_base::BigFloat, z_lphs_base::BigFloat)
    m.sign == 0 && return zero(Complex{BigFloat})
    
    lm = m.z_pow * z_lmag_base
    lp = m.z_pow * z_lphs_base
    m.sign == -1 && (lp += BigFloat(π)) 
    
    @inbounds for (d, e) in m.exps
        lm += e * lmag_table[d]
        lp += e * lphs_table[d]
    end
    
    return exp(lm) * cis(lp)
end

export evaluate_unit_circle, evaluate_analytic


"""
    evaluate_unit_circle(res::CycloResult, theta::Real, ::Type{T}=Complex{BigFloat}; prec=512) where {T}

Evaluates a compiled `CycloResult` strictly along the unit circle where `q = exp(i * theta)`.

It tracks continuous floating-point phases and seamlessly extracts the 
analytic square root of the cyclotomic radical.
"""
function evaluate_unit_circle(res::CycloResult, theta::Real, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
    (res.radical.sign == 0 || res.m_min.sign == 0) && return zero(T)

    return setprecision(BigFloat, prec) do
        theta_big = BigFloat(theta)
        lmag_table, lphs_table = get_unit_circle_table(theta_big, res.max_d)

        z_lmag = zero(BigFloat) 
        z_lphs = theta_big / 2

        sum_val = one(Complex{BigFloat})
        curr_term = one(Complex{BigFloat})
        
        @inbounds for r in res.ratios
            r_val = _project_to_complex_analytic(r, lmag_table, lphs_table, z_lmag, z_lphs)
            curr_term *= r_val
            sum_val += curr_term
        end

        c_mmin = _project_to_complex_analytic(res.m_min, lmag_table, lphs_table, z_lmag, z_lphs)
        c_root = _project_to_complex_analytic(res.root, lmag_table, lphs_table, z_lmag, z_lphs)

        # Extract analytic square root of the radical
        rad_lm = res.radical.z_pow * z_lmag
        rad_lp = res.radical.z_pow * z_lphs
        res.radical.sign == -1 && (rad_lp += BigFloat(π))
        
        @inbounds for (d, e) in res.radical.exps
            rad_lm += e * lmag_table[d]
            rad_lp += e * lphs_table[d]
        end
        
        c_rad_sqrt = exp(rad_lm / 2) * cis(rad_lp / 2) 

        final = c_root * c_rad_sqrt * c_mmin * sum_val
        return T <: Real ? T(real(final)) : T(final)
    end
end


"""
    evaluate_analytic(res::CycloResult, q::Number, ::Type{T}=Complex{BigFloat}; prec=512) where {T}

Evaluates a compiled `CycloResult` for an arbitrary complex parameter `q`.

This function performs the full analytic continuation of the quantum 6j-symbol 
into the SL(2, ℂ) regime. It utilizes a rolling complex iterator to prevent 
O(N log N) exponentiation overhead during the table build phase.
"""
function evaluate_analytic(res::CycloResult, q::Number, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
    (res.radical.sign == 0 || res.m_min.sign == 0) && return zero(T)

    return setprecision(BigFloat, prec) do
        q_big = Complex{BigFloat}(q)
        lmag_table, lphs_table = get_analytic_table(q_big, res.max_d)

        r_q = abs(q_big)
        phi_q = angle(q_big)
        z_lmag = log(r_q) / 2
        z_lphs = phi_q / 2

        sum_val = one(Complex{BigFloat})
        curr_term = one(Complex{BigFloat})
        
        @inbounds for r in res.ratios
            r_val = _project_to_complex_analytic(r, lmag_table, lphs_table, z_lmag, z_lphs)
            curr_term *= r_val
            sum_val += curr_term
        end

        c_mmin = _project_to_complex_analytic(res.m_min, lmag_table, lphs_table, z_lmag, z_lphs)
        c_root = _project_to_complex_analytic(res.root, lmag_table, lphs_table, z_lmag, z_lphs)

        rad_lm = res.radical.z_pow * z_lmag
        rad_lp = res.radical.z_pow * z_lphs
        res.radical.sign == -1 && (rad_lp += BigFloat(π))
        
        @inbounds for (d, e) in res.radical.exps
            rad_lm += e * lmag_table[d]
            rad_lp += e * lphs_table[d]
        end
        
        c_rad_sqrt = exp(rad_lm / 2) * cis(rad_lp / 2) 

        final = c_root * c_rad_sqrt * c_mmin * sum_val
        return T <: Real ? T(real(final)) : T(final)
    end
end

