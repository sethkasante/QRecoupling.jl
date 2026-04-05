# ---------------------------------------------------------------------
# File: projection_discrete.jl
# Discrete Evaluator: SU(2)k Roots of Unity (q = exp(iπ/(k+2)))
# ---------------------------------------------------------------------

using LRUCache

# --- Magnitude Sieve Cache ---
const MAG_SIEVE_CACHE = LRU{Tuple{DataType, Int}, Vector}(maxsize=100)

"""
    get_mag_table(k::Int, max_d::Int, ::Type{T})
Unified cache accessor for the magnitude table log|Φ_d(q²)|.
"""
function get_mag_table(k::Int, max_d::Int, ::Type{T}) where T
    key = (T, k)
    h = k + 2
    # Ensure the table is large enough for the requested max_d
    table = get!(MAG_SIEVE_CACHE, key) do
        build_mag_table(max(max_d, h), k, T)
    end
    
    if length(table) < max_d
        table = build_mag_table(max_d, k, T)
        MAG_SIEVE_CACHE[key] = table
    end
    return table::Vector{T}
end

"""
    build_mag_table(D_max::Int, k::Int, ::Type{T})
Computes log|Φ_d(q²)| at q = exp(iπ/(k+2)). 
Uses log-domain Möbius Inversion (O(D log D)).
"""
function build_mag_table(D_max::Int, k::Int, ::Type{T}) where T
    h = k + 2
    table = Vector{T}(undef, D_max)
    h_inv = inv(T(h))
    
    # 1. Initialize with log|q²ⁿ - 1| = log|2*sin(πn/h)|
    @inbounds for n in 1:D_max
        if n % h == 0
            table[n] = -T(Inf) # Topological pole: Φ_h vanishes
        else
            # We use sinpi for maximum precision near the roots
            table[n] = log(2 * abs(sinpi(n * h_inv)))
        end
    end
    
    # 2. Multiplicative Sieve (Log-Domain Möbius Inversion)
    @inbounds for d in 1:D_max
        val = table[d]
        isinf(val) && continue
        # Forward sieve: subtract divisors from multiples
        for m in (2d):d:D_max
            if !isinf(table[m])
                table[m] -= val
            end
        end
    end
    return table
end

# ==============================================================================
# 1. MONOMIAL PROJECTOR
# ==============================================================================

"""
    project_discrete(m::CyclotomicMonomial, k::Int, ::Type{T}=Float64)
Evaluates [n]_q or dimensions at SU(2)_k.
Returns the real value (signs are tracked, complex phases are ignored).
"""
function project_discrete(m::CyclotomicMonomial, k::Int, ::Type{T}=Float64) where T
    m.sign == 0 && return zero(T)
    h = k + 2
    table = get_mag_table(k, m.max_d, T)
    
    lm = zero(T)
    @inbounds for (d, e) in m.phi_exps
        val = table[d]
        if isinf(val)
            e > 0 && return zero(T) # Topological zero
            e < 0 && throw(DomainError(k, "Topological Pole: Division by zero at level k"))
        end
        lm += e * val
    end
    
    # At roots of unity, q^P is a phase. For real symbols (3j/6j), 
    # q^P contributes to the sign only if P is a multiple of h.
    # However, in SU(2)k, the q_pow in prefactors usually results in a real sign.
    # Standard convention: use the monomial sign.
    return m.sign * exp(lm)
end

# ==============================================================================
# 2. DCR PROJECTOR (Hyper-Optimized)
# ==============================================================================

"""
    project_discrete(dcr::DCR, k::Int, ::Type{T}=Float64; ws=nothing)
The main TQFT evaluation engine. 
Calculates the real-valued symbol at the root of unity.
"""
function project_discrete(dcr::DCR, k::Int, ::Type{T}=Float64; ws=nothing) where T
    # 1. Fetch Basis Table
    table = get_mag_table(k, dcr.max_d, T)

    # 2. Evaluate Monomial Magnitudes
    @inline _log_mag(mono) = begin
        mono.sign == 0 && return -T(Inf)
        lm = zero(T)
        pe = mono.phi_exps
        @inbounds for j in eachindex(pe)
            d, e = pe[j]
            lm += e * table[d]
        end
        return lm
    end

    lm_root = _log_mag(dcr.root)
    lm_rad  = _log_mag(dcr.radical)
    lm_base = _log_mag(dcr.base)

    # Initial magnitude and parity
    curr_l = lm_root + (0.5 * lm_rad) + lm_base
    curr_s = dcr.root.sign * dcr.base.sign
    
    (isinf(curr_l) || curr_s == 0) && return zero(T)

    # 3. Workspace Setup
    nt = length(dcr.ratios) + 1
    _ws = (ws === nothing && T === Float64) ? _WS_F64 : ws
    if _ws === nothing || length(_ws.log_mags) < nt
        _ws = DiscreteBuffer{T}(max(nt, 1024))
    end
    
    l_mags = _ws.log_mags
    signs  = _ws.signs

    max_l = l_mags[1] = curr_l
    signs[1] = Int8(curr_s)

    # 4. HOT LOOP (Zero Allocations)
    for i in 1:length(dcr.ratios)
        r = dcr.ratios[i]
        
        # Inlined _log_mag for max performance
        lm_ratio = zero(T)
        r_exps = r.phi_exps
        @inbounds for j in eachindex(r_exps)
            d, e = r_exps[j]
            lm_ratio += e * table[d]
        end
        
        curr_l += lm_ratio
        curr_s *= r.sign
        
        l_mags[i+1] = curr_l
        signs[i+1]  = Int8(curr_s)
        
        # Track max for Log-Sum-Exp shift
        (curr_l > max_l) && (max_l = curr_l)
    end

    # 5. Final Summation (Log-Sum-Exp)
    fsum = zero(T)
    @inbounds for i in 1:nt
        fsum += signs[i] * exp(l_mags[i] - max_l)
    end

    # Handle vanishing sums (topological zeroes)
    abs(fsum) < 10 * eps(T) && return zero(T)

    return exp(max_l) * fsum
end