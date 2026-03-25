
# ---------------------------------------------------------------------------------------------
# Analytic Evaluator (Analytic continuation  & Complex regimes)
# Maps CycloResults to generic complex parameters (q ∈ ℂ) or the unit circle q = exp(iθ).
# ---------------------------------------------------------------------------------------------

const UNIT_CIRCLE_CACHE  = LRU{Tuple{DataType, Any}, Tuple{Vector, Vector}}(maxsize=256)
const ANALYTIC_CACHE     = LRU{Tuple{DataType, Any}, Tuple{Vector, Vector}}(maxsize=256)

function get_unit_circle_table(theta::T, D_max::Int, ::Type{T}) where {T}
    key = (T, theta)
    if haskey(UNIT_CIRCLE_CACHE, key)
        lmag, lphs = UNIT_CIRCLE_CACHE[key]
        if length(lmag) >= D_max
            return lmag::Vector{T}, lphs::Vector{T}
        end
    end
    lmag, lphs = build_unit_circle_table(D_max, theta, T)
    UNIT_CIRCLE_CACHE[key] = (lmag, lphs)
    return lmag, lphs
end

function get_analytic_table(q::Complex{T}, D_max::Int, ::Type{T}) where {T}
    key = (T, q)
    if haskey(ANALYTIC_CACHE, key)
        lmag, lphs = ANALYTIC_CACHE[key]
        if length(lmag) >= D_max
            return lmag::Vector{T}, lphs::Vector{T}
        end
    end
    lmag, lphs = build_analytic_table(D_max, q, T)
    ANALYTIC_CACHE[key] = (lmag, lphs)
    return lmag, lphs
end

function build_unit_circle_table(D_max::Int, theta::T, ::Type{T}) where {T}
    V_mag = zeros(T, D_max)
    V_phs = zeros(T, D_max) 
    V_valid = trues(D_max)
    
    half_theta = theta / 2
    pi_val = T(π) 
    
    for n in 1:D_max
        val_sin = 2 * sin(T(n) * half_theta)
        
        if iszero(val_sin)
            V_valid[n] = false
        else
            V_mag[n] = log(abs(val_sin))
            p_n = T(n) * half_theta + (pi_val / 2) 
            if val_sin < 0
                p_n += pi_val
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

function build_analytic_table(D_max::Int, q::Complex{T}, ::Type{T}) where {T}
    V_mag = zeros(T, D_max)
    V_phs = zeros(T, D_max)
    V_valid = trues(D_max)
    
    q_curr = q
    one_q = one(q)
    for n in 1:D_max
        term = q_curr - one_q
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

# --------------------------------------------
# High-Performance Continuous Projections
# --------------------------------------------

@inline function _project_to_complex_analytic(m::CycloMonomial, lmag_table::Vector{T}, lphs_table::Vector{T}, 
                                              z_lmag_base::T, z_lphs_base::T, ::Type{T}) where {T}
    m.sign == 0 && return zero(Complex{T})
    
    lm = m.z_pow * z_lmag_base
    lp = m.z_pow * z_lphs_base
    m.sign == -1 && (lp += T(π)) 
    
    # Pure SIMD math loop (No d >= h check because k doesn't exist here!)
    exps = m.exps
    @inbounds @simd ivdep for i in eachindex(exps)
        d, e = exps[i]
        lm += e * lmag_table[d]
        lp += e * lphs_table[d]
    end
    
    return exp(lm) * cis(lp)
end


function _eval_unit_circle_core(res::CycloResult, theta::T, ::Type{T}) where {T}
    lmag_table, lphs_table = get_unit_circle_table(theta, res.max_d, T)

    z_lmag = zero(T) 
    z_lphs = theta / 2

    sum_val = one(Complex{T})
    curr_term = one(Complex{T})
    
    @inbounds for r in res.ratios
        r_val = _project_to_complex_analytic(r, lmag_table, lphs_table, z_lmag, z_lphs, T)
        curr_term *= r_val
        sum_val += curr_term
    end

    c_mmin = _project_to_complex_analytic(res.base_term, lmag_table, lphs_table, z_lmag, z_lphs, T)
    c_root = _project_to_complex_analytic(res.root, lmag_table, lphs_table, z_lmag, z_lphs, T)

    rad_lm = res.radical.z_pow * z_lmag
    rad_lp = res.radical.z_pow * z_lphs
    res.radical.sign == -1 && (rad_lp += T(π))
    
    exps = res.radical.exps
    @inbounds @simd ivdep for i in eachindex(exps)
        d, e = exps[i]
        rad_lm += e * lmag_table[d]
        rad_lp += e * lphs_table[d]
    end
    
    c_rad_sqrt = exp(rad_lm / 2) * cis(rad_lp / 2) 

    return c_root * c_rad_sqrt * c_mmin * sum_val
end

"""
    evaluate_unit_circle(res::CycloResult, theta::Real, ::Type{T}=Complex{BigFloat}; prec=512) where {T}

Evaluates a compiled `CycloResult` strictly along the unit circle where `q = exp(i * theta)`.

It tracks continuous floating-point phases and seamlessly extracts the 
analytic square root of the cyclotomic radical.
"""
function evaluate_unit_circle(res::CycloResult, theta::Real, ::Type{T}=ComplexF64; prec=512) where {T}
    (res.radical.sign == 0 || res.base_term.sign == 0) && return zero(T)

    RealT = T <: Complex ? real(T) : T
    theta_T = RealT(theta)

    final = if RealT == BigFloat
        setprecision(BigFloat, prec) do 
            _eval_unit_circle_core(res, theta_T, RealT)
        end
    else
        _eval_unit_circle_core(res, theta_T, RealT)
    end
    
    return T <: Real ? T(real(final)) : T(final)
end



function _eval_analytic_core(res::CycloResult, q::Complex{T}, ::Type{T}) where {T}
    lmag_table, lphs_table = get_analytic_table(q, res.max_d, T)

    r_q = abs(q)
    phi_q = angle(q)
    z_lmag = log(r_q) / 2
    z_lphs = phi_q / 2

    sum_val = one(Complex{T})
    curr_term = one(Complex{T})
    
    @inbounds for r in res.ratios
        r_val = _project_to_complex_analytic(r, lmag_table, lphs_table, z_lmag, z_lphs, T)
        curr_term *= r_val
        sum_val += curr_term
    end

    c_mmin = _project_to_complex_analytic(res.base_term, lmag_table, lphs_table, z_lmag, z_lphs, T)
    c_root = _project_to_complex_analytic(res.root, lmag_table, lphs_table, z_lmag, z_lphs, T)

    rad_lm = res.radical.z_pow * z_lmag
    rad_lp = res.radical.z_pow * z_lphs
    res.radical.sign == -1 && (rad_lp += T(π))
    
    exps = res.radical.exps
    @inbounds @simd ivdep for i in eachindex(exps)
        d, e = exps[i]
        rad_lm += e * lmag_table[d]
        rad_lp += e * lphs_table[d]
    end
    
    c_rad_sqrt = exp(rad_lm / 2) * cis(rad_lp / 2) 

    return c_root * c_rad_sqrt * c_mmin * sum_val
end


"""
    evaluate_analytic(res::CycloResult, q::Number, ::Type{T}=Complex{BigFloat}; prec=512) where {T}

Evaluates a compiled `CycloResult` for an arbitrary complex parameter `q`.

This function performs the full analytic continuation of the quantum 6j-symbol 
into the SL(2, ℂ) regime. 
"""
function evaluate_analytic(res::CycloResult, q::Number, ::Type{T}=ComplexF64; prec=512) where {T}
    (res.radical.sign == 0 || res.base_term.sign == 0) && return zero(T)

    RealT = T <: Complex ? real(T) : T
    q_T = Complex{RealT}(q)

    final = if RealT == BigFloat
        setprecision(BigFloat, prec) do 
            _eval_analytic_core(res, q_T, RealT)
        end
    else
        _eval_analytic_core(res, q_T, RealT)
    end
    
    return T <: Real ? T(real(final)) : T(final)
end

