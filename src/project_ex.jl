# ---------------------------------------------------------------------
# File: projection_exact.jl
# Exact SU(2)k Evaluator using Nemo.jl Cyclotomic Fields.
# ---------------------------------------------------------------------

using Nemo
using LRUCache

# --- Basis Cache ---
# Stores (V_exact, V_inv) per level k. Rebuilds if max_d increases.
const EXACT_PHI_CACHE = LRU{Int, Any}(maxsize = 100)

"""
    _phi_exact_table(D_max::Int, k::Int, ζ::nf_elem)
Constructs the cyclotomic basis Φ_d(ζ²) and its inverses in Q(ζ).
"""
function _phi_exact_table(D_max::Int, k::Int, ζ::nf_elem)
    h = k + 2
    V_exact = Vector{nf_elem}(undef, D_max)
    V_inv   = Vector{nf_elem}(undef, D_max)
    
    # Initialize with (q²ⁿ - 1) where q = ζ
    @inbounds for n in 1:D_max
        V_exact[n] = (n % h == 0) ? zero(ζ) : ζ^(2n) - one(ζ)
    end
    
    # Sieve-based multiplicative Möbius inversion
    @inbounds for d in 1:D_max
        iszero(V_exact[d]) && continue
        for m in (2d):d:D_max
            !iszero(V_exact[m]) && (V_exact[m] = divexact(V_exact[m], V_exact[d]))
        end
    end
    
    # Precompute inverses for zero-division hot loops
    @inbounds for d in 1:D_max
        V_inv[d] = iszero(V_exact[d]) ? zero(ζ) : inv(V_exact[d])
    end
    
    return V_exact, V_inv
end

@inline function get_phi_exact_table(D_max::Int, k::Int, ζ::nf_elem)
    if !haskey(EXACT_PHI_CACHE, k) || length((EXACT_PHI_CACHE[k])[1]) < D_max
        EXACT_PHI_CACHE[k] = _phi_exact_table(D_max, k, ζ)
    end
    return EXACT_PHI_CACHE[k]::Tuple{Vector{nf_elem}, Vector{nf_elem}}
end

# ==============================================================================
# 1. MONOMIAL PROJECTOR
# ==============================================================================

"""
    project_exact(m::CyclotomicMonomial, k::Int)
Evaluates a single monomial exactly. Returns a Nemo nf_elem.
"""
function project_exact(m::CyclotomicMonomial, k::Int)
    m.sign == 0 && return zero(cyclotomic_field(2*(k+2), "ζ")[2])
    h = k + 2
    K, ζ = cyclotomic_field(2h, "ζ")
    V_exact, V_inv = get_phi_exact_table(m.max_d, k, ζ)
    
    val = one(ζ)
    @inbounds for (d, e) in m.phi_exps
        if d == h
            e > 0 && return zero(ζ)
            e < 0 && throw(DomainError(k, "Topological pole at level k"))
            continue
        end
        base = e > 0 ? V_exact[d] : V_inv[d]
        abs_e = abs(e)
        if abs_e == 1; val *= base
        elseif abs_e == 2; val *= (base * base)
        else; val *= base^abs_e end
    end
    
    # Apply q^P and sign
    res = val * ζ^m.q_pow
    return m.sign == 1 ? res : -res
end

# ==============================================================================
# 2. DCR PROJECTOR
# ==============================================================================

"""
    project_exact(dcr::DCR, k::Int)
Evaluates a full DCR series exactly. Returns a CycloExactResult.
"""
function project_exact(dcr::DCR, k::Int)
    h = k + 2
    K, ζ = cyclotomic_field(2h, "ζ")
    
    if dcr.radical.sign == 0 || dcr.base.sign == 0
        return CycloExactResult(k, ZERO_MONOMIAL, zero(ζ))
    end
    
    V_exact, V_inv = get_phi_exact_table(dcr.max_d, k, ζ)
    
    # --- HOT LOOP: Hypergeometric Sum ---
    sum_val = one(ζ)
    curr_term = one(ζ)
    
    for r in dcr.ratios
        # Internal optimized monomial projection to avoid field lookups
        r_val = _project_monomial_nemo_internal(r, V_exact, V_inv, ζ, h)
        curr_term *= r_val
        sum_val += curr_term
    end

    # Apply root prefactor and base term
    m_base_val = _project_monomial_nemo_internal(dcr.base, V_exact, V_inv, ζ, h)
    m_root_val = _project_monomial_nemo_internal(dcr.root, V_exact, V_inv, ζ, h)
    
    return CycloExactResult(k, dcr.radical, m_root_val * m_base_val * sum_val)
end

# Internal helper to keep the loops clean
@inline function _project_monomial_nemo_internal(m, V_exact, V_inv, ζ, h)
    m.sign == 0 && return zero(ζ)
    val = one(ζ)
    @inbounds for (d, e) in m.phi_exps
        (d == h) && (e > 0 ? (return zero(ζ)) : continue)
        base = e > 0 ? V_exact[d] : V_inv[d]
        abs_e = abs(e)
        if abs_e == 1; val *= base
        elseif abs_e == 2; val *= (base * base)
        else; val *= base^abs_e end
    end
    res = val * ζ^m.q_pow
    return m.sign == 1 ? res : -res
end

# ==============================================================================
# 3. ARITHMETIC FOR IDENTITY VERIFICATION
# ==============================================================================

# Allows adding two exact results if they share the same radical.
# Critical for Pentagons: Σ (F * F * F) - (F * F) == 0
function Base.:+(a::CycloExactResult, b::CycloExactResult)
    a.k != b.k && error("Level mismatch")
    if a.radical == b.radical
        return CycloExactResult(a.k, a.radical, a.sum_factor + b.sum_factor)
    elseif iszero(a.sum_factor)
        return b
    elseif iszero(b.sum_factor)
        return a
    else
        # If radicals differ, verification is done by checking if Sum(res_i) is zero.
        # This usually requires moving the radicals into the field (if they are perfect squares).
        error("Cannot sum exact results with different radicals algebraically.")
    end
end

Base.:-(a::CycloExactResult, b::CycloExactResult) = a + CycloExactResult(b.k, b.radical, -b.sum_factor)
Base.iszero(res::CycloExactResult) = iszero(res.sum_factor)

# Scalar multiplication
Base.:*(c::Number, res::CycloExactResult) = CycloExactResult(res.k, res.radical, c * res.sum_factor)