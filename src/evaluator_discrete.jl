# -------------------------------------------------------------------------------
# Discrete Evaluator: SU(2)_k Roots of Unity
# Maps symbolic CycloResults into high-precision floats at exact roots of unity.
# -------------------------------------------------------------------------------

# Cache the magnitude table
const GLOBAL_SIEVE_CACHE = LRU{Int, Vector{BigFloat}}(maxsize=4096)

export build_phi_table, get_phi_table

"""
    build_phi_table(D_max::Int, k::Int)

Computes the log-magnitude for cyclotomic polynomials Φ_d(q) using a 
Fast Möbius Transform. Phase tracking is analytically omitted as admissible 
SU(2)_k ratios are mathematically guaranteed to be positive reals.
"""
function build_phi_table(D_max::Int, k::Int)
    h = k + 2
    V_mag = zeros(BigFloat, D_max)
    V_valid = trues(D_max) 
    
    h_big = BigFloat(h)
    
    for n in 1:D_max
        if n % h == 0
            V_valid[n] = false
            continue
        end

        # We strictly take the absolute value, magnitudes only
        val_sin = 2 * sinpi(BigFloat(n) / h_big)
        V_mag[n] = log(abs(val_sin))
    end
    
    # Fast Möbius Transform for Cyclotomic factor log-magnitudes
    @inbounds for d in 1:D_max
        !V_valid[d] && continue
        for m in (2 * d):d:D_max
            if V_valid[m]
                V_mag[m] -= V_mag[d]
            end
        end
    end
    
    return V_mag
end

function get_phi_table(k::Int, D_max::Int)
    if haskey(GLOBAL_SIEVE_CACHE, k)
        lmag = GLOBAL_SIEVE_CACHE[k]
        if length(lmag) >= D_max
            return lmag
        end
    end
    
    lmag = build_phi_table(max(D_max, k + 2), k)
    GLOBAL_SIEVE_CACHE[k] = lmag
    return lmag
end

# -------------------------------------------
# High-Performance Discrete Projection
# -------------------------------------------

"""
    _project_to_real(m::CycloMonomial, lmag_table::Vector{BigFloat}, h::Int)

The ultra-fast hot loop. Evaluates the cyclotomic log-magnitudes. 
Gracefully detects topological zeros if d >= h.
"""
@inline function _project_to_real(m::CycloMonomial, lmag_table::Vector{BigFloat}, h::Int)
    m.sign == 0 && return zero(BigFloat)
    
    lm = zero(BigFloat)
    
    @inbounds for (d, e) in m.exps
        if d >= h
            e > 0 && return zero(BigFloat) # Topological zero (Admissibility violated or pole hit)
            e < 0 && throw(DomainError(h-2, "Topological pole: Level k=$(h-2)."))
            continue
        end
        lm += e * lmag_table[d]
    end
    
    # Exponentiate the magnitude and apply the analytical algebraic sign
    val = exp(lm)
    return m.sign == -1 ? -val : val
end

export evaluate_level

"""
    evaluate_level(res::CycloResult, k::Int, ::Type{T}=Float64; prec=512)

Evaluates the quantum symbol at the exact discrete SU(2)_k root of unity.
"""
function evaluate_level(res::CycloResult, k::Int, ::Type{T}=Float64; prec=512) where {T}
    (res.radical.sign == 0 || res.base_term.sign == 0) && return zero(T)

    return setprecision(BigFloat, prec) do
        h = k + 2
        lmag_table = get_phi_table(k, res.max_d)

        sum_val = one(BigFloat)
        curr_term = one(BigFloat)
        
        # The Hot Loop: Now purely vectorized float addition
        @inbounds for r in res.ratios
            r_val = _project_to_real(r, lmag_table, h)
            curr_term *= r_val
            sum_val += curr_term
        end

        c_base = _project_to_real(res.base_term, lmag_table, h)
        c_root = _project_to_real(res.root, lmag_table, h)
        
        # Evaluate radical directly (Delta^2 is rigorously strictly positive)
        p_lm = zero(BigFloat)
        @inbounds for (d, e) in res.radical.exps
            if d >= h
                return zero(T) # Topological zero safety net
            end
            p_lm += e * lmag_table[d]
        end
        c_rad_sqrt = exp(p_lm / 2) 

        final = c_root * c_rad_sqrt * c_base * sum_val
        return T <: Complex ? T(final, 0) : T(final)
    end
end

"""
    evaluate_level(m::CycloMonomial, k::Int, ::Type{T}=Float64; prec=512)
"""
function evaluate_level(m::CycloMonomial, k::Int, ::Type{T}=Float64; prec=512) where {T}
    m.sign == 0 && return zero(T)

    return setprecision(BigFloat, prec) do
        h = k + 2
        max_d = isempty(m.exps) ? 1 : maximum(keys(m.exps))
        lmag_table = get_phi_table(k, max_d)
        val_bf = _project_to_real(m, lmag_table, h)
        return T <: Complex ? T(val_bf, 0) : T(val_bf)
    end
end



"""
    cyclo_to_numeric(res::CycloResult, [T=Float64]; k=nothing, q=nothing, theta=nothing, prec=512)

Internal numeric projection engine. Maps the deferred graph into high-performance 
floating-point arithmetic (hardware or arbitrary precision).
"""
function cyclo_to_numeric(res::CycloResult, ::Type{T}=Float64; k=nothing, q=nothing, theta=nothing, prec=512) where {T}
    targets_defined = (!isnothing(k)) + (!isnothing(q)) + (!isnothing(theta))
    if targets_defined != 1
        throw(ArgumentError("Specify exactly one evaluation target: k (integer), theta (real), or q (complex)."))
    end

    if !isnothing(k)
        return evaluate_level(res, Int(k), T; prec=prec)
    elseif !isnothing(theta)
        return evaluate_unit_circle(res, theta, T; prec=prec)
    else
        return evaluate_analytic(res, q, T; prec=prec)
    end
end
