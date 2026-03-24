
# -------------------------------------------------------------------------------
# Discrete Evaluator: SU(2)_k Roots of Unity
# Maps symbolic CycloResults into high-precision floats at exact roots of unity.
# -------------------------------------------------------------------------------

const GLOBAL_SIEVE_CACHE = LRU{Int, Tuple{Vector{BigFloat}, Vector{Int}}}(maxsize=4096)

export build_phi_table, get_phi_table

"""
    build_phi_table(D_max::Int, k::Int)

Computes the log-magnitude and exact integer phase for cyclotomic polynomials Φ_d(q) 
using a Fast Möbius Transform with an O(1) Boolean mask.
Returns magnitude and phase (integer multiple of π/2h): |logΦ_d|, arg(Φ_d).
"""
function build_phi_table(D_max::Int, k::Int)
    h = k + 2
    V_mag = zeros(BigFloat, D_max)
    V_phs = zeros(Int, D_max) 
    V_valid = trues(D_max) 
    
    h_big = BigFloat(h)
    
    for n in 1:D_max
        if n % h == 0
            V_valid[n] = false
            continue
        end

        val_sin = 2 * sinpi(BigFloat(n) / h_big)
        V_mag[n] = log(abs(val_sin))
        
        # Exact integer phase: P_n = 2n + h (in units of π/2h)
        p_n = 2n + h
        if val_sin < 0
            p_n += 2h 
        end
        V_phs[n] = p_n
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

function get_phi_table(k::Int, D_max::Int)
    if haskey(GLOBAL_SIEVE_CACHE, k)
        lmag, lphs = GLOBAL_SIEVE_CACHE[k]
        if length(lmag) >= D_max
            return lmag, lphs
        end
    end
    
    lmag, lphs = build_phi_table(max(D_max, k + 2), k)
    GLOBAL_SIEVE_CACHE[k] = (lmag, lphs)
    return lmag, lphs
end

# -------------------------------------------
# High-Performance Discrete Projection
# -------------------------------------------

@inline function _project_to_real(m::CycloMonomial, lmag_table::Vector{BigFloat}, lphs_table::Vector{Int}, h::Int)
    m.sign == 0 && return zero(BigFloat)
    
    lm = zero(BigFloat)
    P_total = m.z_pow * 2 
    m.sign == -1 && (P_total += 2h) 
    
    @inbounds for (d, e) in m.exps
        if d == h
            e > 0 && return zero(BigFloat) 
            e < 0 && throw(DomainError(h-2, "Topological pole: Level k=$(h-2)."))
            continue
        end
        lm += e * lmag_table[d]
        P_total += e * lphs_table[d]
    end
    
    rem = mod(P_total, 4h)
    
    if rem == 0
        return exp(lm) 
    elseif rem == 2h
        return -exp(lm) 
    else
        error("Fatal: Non-real phase detected in topological symbol. P_total = $P_total, rem = $rem")
    end
end

export evaluate_level

"""
    evaluate_level(res::CycloResult, k::Int, ::Type{T}=Float64; prec=512)

Evaluates the quantum symbol at the exact discrete SU(2)_k root of unity.
Employs pure integer phase tracking to guarantee immunity against Float64 sign corruption.
"""
function evaluate_level(res::CycloResult, k::Int, ::Type{T}=Float64; prec=512) where {T}
    (res.radical.sign == 0 || res.base_term.sign == 0) && return zero(T)

    return setprecision(BigFloat, prec) do
        h = k + 2
        lmag_table, lphs_table = get_phi_table(k, res.max_d)

        sum_val = one(BigFloat)
        curr_term = one(BigFloat)
        
        @inbounds for r in res.ratios
            r_val = _project_to_real(r, lmag_table, lphs_table, h)
            curr_term *= r_val
            sum_val += curr_term
        end

        c_mmin = _project_to_real(res.base_term, lmag_table, lphs_table, h)
        c_root = _project_to_real(res.root, lmag_table, lphs_table, h)
        
        # Evaluate radical directly (Delta^2 is rigorously strictly positive)
        p_lm = zero(BigFloat)
        @inbounds for (d, e) in res.radical.exps
            d != h && (p_lm += e * lmag_table[d])
        end
        c_rad_sqrt = exp(p_lm / 2) 

        final = c_root * c_rad_sqrt * c_mmin * sum_val
        return T <: Complex ? T(final, 0) : T(final)
    end
end

"""
    evaluate_level(m::CycloMonomial, k::Int, ::Type{T}=Float64; prec=512)

Projects a single CycloMonomial into a numeric floating-point type at the specified root of unity.
"""
function evaluate_level(m::CycloMonomial, k::Int, ::Type{T}=Float64; prec=512) where {T}
    m.sign == 0 && return zero(T)

    return setprecision(BigFloat, prec) do
        h = k + 2
        # Determine the maximum polynomial degree needed for this monomial
        max_d = isempty(m.exps) ? 1 : maximum(keys(m.exps))
        
        # Fetch the cached LSE tables
        lmag_table, lphs_table = get_phi_table(k, max_d)
        
        # Project using your internal optimized projector
        val_bf = _project_to_real(m, lmag_table, lphs_table, h)
        
        return T <: Complex ? T(val_bf, 0) : T(val_bf)
    end
end