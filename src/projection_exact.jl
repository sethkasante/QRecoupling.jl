
# ------------------------------------------------
# Exact SU(2)k evaluator using Nemo.jl  
# (in cyclotomic fields Q(ζ)
# ------------------------------------------------


# -------------------------------------------------------------------------------
# !!! note  Performance note for high level k:
# 
# This exact engine maps a deferred cyclotomic representation (DCR) into Nemo.jl 
# cyclotomic number fields. It is optimized to use zero-division hot loops.
# However, Computer Algebra Systems (CAS) fundamentally consume exponential 
# memory as the cyclotomic degree grows. 
# 
# As a rule of thumb, this exact engine is incredibly fast and stable for k < 500. 
# Pushing past k = 1000 may cause RAM explosion or severe slowdowns due to 
# dense polynomial arithmetic. If you are working in ultra-high level regimes. 
# Consider using the `mode=:discrete` instead for exact evaluation.
# -------------------------------------------------------------------------------



# Stores (V_exact, V_inv) per level k. Rebuilds if max_d increases.
const EXACT_PHI_CACHE = LRU{Int, Any}(maxsize = 500)
const EXACT_PHI_LOCK = ReentrantLock()

"""
    _phi_exact_table(D_max::Int, k::Int, ζ::T)
Constructs the cyclotomic basis Φ_d(ζ²) and its inverses in Q(ζ).
"""
function _phi_exact_table(D_max::Int, k::Int, ζ::T) where T
    h = k + 2
    V_exact = Vector{T}(undef, D_max)
    V_inv   = Vector{T}(undef, D_max)
    
    # Initialize with (q^2n - 1) where q = ζ
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

@inline function get_phi_exact_table(D_max::Int, k::Int, ζ::T) where T
    lock(EXACT_PHI_LOCK) do
        if !haskey(EXACT_PHI_CACHE, k) || length((EXACT_PHI_CACHE[k])[1]) < D_max
            EXACT_PHI_CACHE[k] = _phi_exact_table(D_max, k, ζ)
        end
        return EXACT_PHI_CACHE[k]::Tuple{Vector{T}, Vector{T}}
    end
end


# Internal helper for the loops
@inline function _project_monomial_nemo_internal(m, V_exact, V_inv, ζ, h)
    m.sign == 0 && return zero(ζ)
    val = one(ζ)
    @inbounds for (d, e) in m.phi_exps
         if d == h
            e > 0 && return zero(ζ)
            e < 0 && throw(DomainError(k, "Topological pole at level k"))
            continue
        end

        base = e > 0 ? V_exact[d] : V_inv[d]
        abs_e = abs(e)
        if abs_e == 1 
            val *= base
        elseif abs_e == 2
            val *= (base * base)
        else
            val *= base^abs_e 
        end
    end
    res = val * ζ^m.q_pow
    return m.sign == 1 ? res : -res
end

# --- Project monomial and DCR results ----

"""
    project_exact(m::CyclotomicMonomial, k::Int)
Evaluates a single monomial exactly. Returns a Nemo `nf_elem`.
"""
function project_exact(m::CyclotomicMonomial, k::Int)
    h = k + 2
    K, ζ = cyclotomic_field(2h, "ζ")
    m.sign == 0 && return zero(ζ)
    V_exact, V_inv = get_phi_exact_table(m.max_d, k, ζ)

    return _project_monomial_nemo_internal(m, V_exact, V_inv, ζ, h)
end


# -- Project cyclotomic monomials and DCR to exact cyclotomic field ---

"""
    project_exact(dcr::DCR, k::Int)
Evaluates a DCR series exactly for discrete level `k`. Returns a CompositeExactResult.
"""
function project_exact(dcr::DCR, k::Int)
    h = k + 2
    K, ζ = cyclotomic_field(2h, "ζ")
    
    if dcr.radical.sign == 0 || dcr.base.sign == 0
        return CompositeExactResult(k, ZERO_MONOMIAL, zero(ζ))
    end
    
    V_exact, V_inv = get_phi_exact_table(dcr.max_d, k, ζ)
    
    # --- Loop: Hypergeometric Sum ---
    sum_val = one(ζ)
    curr_term = one(ζ)
    
    for r in dcr.ratios
        # internal monomial projection
        r_val = _project_monomial_nemo_internal(r, V_exact, V_inv, ζ, h)
        curr_term *= r_val
        sum_val += curr_term
    end

    # root prefactor and base term
    m_base_val = _project_monomial_nemo_internal(dcr.base, V_exact, V_inv, ζ, h)
    m_root_val = _project_monomial_nemo_internal(dcr.root, V_exact, V_inv, ζ, h)
    
    return CompositeExactResult(k, dcr.radical, m_root_val * m_base_val * sum_val)
end




"""
    evaluate_exact(res::CompositeExactResult, [T=ComplexF64])
Projects the deferred cyclotomic exact result into a complex/real number.
Just for consistency checks
"""
function evaluate_exact(res::CompositeExactResult, ::Type{T}=ComplexF64) where T
    h = res.k + 2
    target_z = cispi(one(BigFloat) / h)
    
    # 1. Horner evaluation for the Nemo polynomial (safe BigFloat casting)
    function _horner(poly, z)
        deg = degree(parent(poly))
        val = Complex{BigFloat}(0)
        for i in (deg-1):-1:0
            c = coeff(poly, i)
            c_bf = BigFloat(numerator(c)) / BigFloat(denominator(c))
            val = val * z + c_bf
        end
        return val
    end
    
    # Evaluate the exact polynomial sum
    B = _horner(res.sum_factor, target_z) 
    
    # evaluate the radical (CyclotomicMonomial) directly
    rad_val = project_discrete(res.radical, res.k )
    
    A = sqrt(max(zero(BigFloat), real(rad_val)))
    val = A * B
    
    return T <: Real ? T(real(val)) : T(val)
end