
# ------------------------------------------------
# Exact Nemo Evaluator (in cyclotomic fields Q(ζ)
# (Zero-Division Architecture)
# ------------------------------------------------


# -------------------------------------------------------------------------------
# !!! note  PERFORMANCE NOTE FOR HIGH LEVELS (k):
# 
# This exact engine maps our sparse cyclotomic arrays into rigorous Nemo.jl 
# cyclotomic number fields. While we have optimized this to use zero-division 
# hot loops, Computer Algebra Systems (CAS) fundamentally consume exponential 
# memory as the cyclotomic degree grows. 
# 
# As a rule of thumb, this exact engine is incredibly fast and stable for k < 700. 
# Pushing past k = 1000 may cause RAM explosion or severe slowdowns due to 
# dense polynomial arithmetic. If you are working in ultra-high level regimes, 
# strongly consider using the `:numeric` mode with BigFloats instead.
# -------------------------------------------------------------------------------


const EXACT_PHI_CACHE = LRU{Int, Any}(maxsize = 1000)
const EMPTY_MONOMIAL = CycloMonomial(1, 0, Pair{Int,Int}[])


"""
    _phi_exact_table(D_max::Int, k::Int, z::T) where {T}

Internal builder for the exact cyclotomic polynomial table.
Uses a multiplicative Möbius inversion to iteratively construct the exact 
algebraic values of Φ_d(z) in the cyclotomic field ζ. 

Crucially, this also precomputes and returns the *algebraic inverses* (`V_inv`), 
which completely eliminates the need for expensive matrix-inversion divisions 
during the hypergeometric summation loop. 
"""
function _phi_exact_table(D_max::Int, k::Int, z::T) where {T}
    h = k + 2
    V_exact = Vector{T}(undef, D_max)
    V_inv   = Vector{T}(undef, D_max)
    
    @inbounds for n in 1:D_max
        V_exact[n] = (n % h == 0) ? zero(z) : z^(2n) - one(z)
    end
    
    @inbounds for d in 1:D_max
        iszero(V_exact[d]) && continue
        for m in (2 * d):d:D_max
            if !iszero(V_exact[m])
                V_exact[m] = divexact(V_exact[m], V_exact[d])
            end
        end
    end
    
    # Precompute all inverses here to eliminate division in the hot loop
    @inbounds for d in 1:D_max
        if !iszero(V_exact[d])
            V_inv[d] = inv(V_exact[d])
        else
            V_inv[d] = zero(z)
        end
    end
    
    return V_exact, V_inv
end

"""
    get_phi_exact_table(D_max::Int, k::Int, z::T) where {T}

Retrieves the precomputed `(V_exact, V_inv)` table from the global LRU cache.
If the requested `D_max` exceeds the currently cached size, it safely expands 
and recompiles the cache for the current level `k`.
"""
@inline function get_phi_exact_table(D_max::Int, k::Int, z::T) where {T}
    # Corrected Type Assertion: We are caching a Tuple of two Vectors!
    if !haskey(EXACT_PHI_CACHE, k) || length((EXACT_PHI_CACHE[k]::Tuple{Vector{T}, Vector{T}})[1]) < D_max
        EXACT_PHI_CACHE[k] = _phi_exact_table(max(D_max, k + 2), k, z)
    end
    return EXACT_PHI_CACHE[k]::Tuple{Vector{T}, Vector{T}} 
end


"""
    _project_ratio_nemo(m::CycloMonomial, V_exact::Vector{T}, V_inv::Vector{T}, z::T, h::Int, z_zero::T, z_one::T) where {T}

Maps a sparse `CycloMonomial` into an exact Nemo number field element.
Small powers are explicitly unrolled to prevent generic allocation overhead.
"""
@inline function _project_ratio_nemo(m::CycloMonomial, V_exact::Vector{T}, V_inv::Vector{T}, z::T, h::Int, z_zero::T, z_one::T) where {T}
    # Fixed early exit return type
    m.sign == 0 && return z_zero
    
    val = z_one
    @inbounds for (d, e) in m.exps
        if d == h
            e > 0 && return z_zero
            e < 0 && throw(DomainError(h-2, "Topological pole: Level k=$(h-2)"))
            continue
        end
        
        # Multiply by V_exact for positive powers, V_inv for negative powers
        base = e > 0 ? V_exact[d] : V_inv[d]
        abs_e = abs(e)
        
        # Explicit unrolling for small powers to avoid generic `^` allocations
        if abs_e == 1
            val *= base
        elseif abs_e == 2
            val *= base * base
        else
            val *= base^abs_e
        end
    end
    
    val *= z^(m.z_pow)
    return m.sign == 1 ? val : -val
end



export evaluate_level_exact


"""
    evaluate_level_exact(res::CycloResult, k::Int)

The master exact evaluator for quantum recoupling symbols.
Takes a compiled, unevaluated `CycloResult` and projects it into the precise 
algebraic cyclotomic field for SU(2)_k.

# Returns:
An `ExactResult{T}` containing:
1. `pref_rad`: The strictly square-free remainder (kept symbolic for fast O(1) multiplication).
2. `sum_part`: The evaluated rational sum (a Nemo field element).

# Architecture:
Extracts the perfect algebraic squares during the construction phase, this 
evaluator bypasses CAS square-root factorization entirely. The hypergeometric 
loop itself executes using zero memory allocations and zero algebraic divisions.
"""
function evaluate_level_exact(res::CycloResult, k::Int)
    h = k + 2
    _, z = cyclotomic_field(2 * h, "ζ")
    z_zero, z_one = zero(z), one(z)
    
    if res.pref_rad.sign == 0 || res.m_min.sign == 0
        return ExactResult(k, EMPTY_MONOMIAL, z_zero)
    end
    
    V_exact, V_inv = get_phi_exact_table(res.max_d, k, z)
    
    # Notice: ZERO divisions in this hot loop! 
    sum_val, curr_term = z_one, z_one
    for r in res.ratios
        r_val = _project_ratio_nemo(r, V_exact, V_inv, z, h, z_zero, z_one)
        curr_term *= r_val
        sum_val += curr_term
    end

    c_mmin = _project_ratio_nemo(res.m_min, V_exact, V_inv, z, h, z_zero, z_one)
    exact_root = _project_ratio_nemo(res.pref_root, V_exact, V_inv, z, h, z_zero, z_one)
    
    return ExactResult(k, res.pref_rad, exact_root * c_mmin * sum_val)
end

export evaluate_level_exact




# ==============================================================================
# 5. Proof Scripts
# ==============================================================================
function exact_qdim(j::Spin, k::Int)
    h = k + 2
    K, z = cyclotomic_field(2 * h, "ζ")
    
    num = z^(Int(4j + 2)) - K(1)
    den = z^2 - K(1)
    dim_val = divexact(num, den) * z^(-Int(2j))
    
    # Quantum dimension has no prefactor root, so we pass EMPTY_MONOMIAL (i.e. '1')
    return ExactResult(k, EMPTY_MONOMIAL, dim_val)
end

export prove_orthogonality

function prove_orthogonality(j1::Spin, j2::Spin, j4::Spin, j5::Spin, j3::Spin, k::Int)
    println("======================================================")
    println(" Proving Exact Quantum Orthogonality (Bubble Move)")
    println("======================================================")
    
    h = k + 2
    K, z = cyclotomic_field(2 * h, "ζ")
    
    # Accumulate starting at 0 (empty monomial remainder)
    total_sum = ExactResult(k, EMPTY_MONOMIAL, K(0))
    
    x_min = max(abs(j1 - j5), abs(j2 - j4))
    x_max = min(j1 + j5, j2 + j4, k - max(j1+j5, j2+j4))
    
    println("Summing over internal spin x from $x_min to $x_max at level k=$k...")
    
    for x in x_min:1:x_max
        res_A = q6j_cyclo(j1, j2, j3, j4, j5, x)
        res_B = q6j_cyclo(j1, j2, j3, j4, j5, x)
        
        exact_A = evaluate_level_exact(res_A, k)
        exact_B = evaluate_level_exact(res_B, k)
        
        dim_x = exact_qdim(x, k)
        
        # When exact_A * exact_B happens, the remainder squares itself,
        # becomes perfect, and gets automatically sucked into the sum_part!
        term = dim_x * exact_A * exact_B 
        
        total_sum = total_sum + term
    end
    
    dim_j3 = exact_qdim(j3, k)
    theoretical_rhs = ExactResult(k, EMPTY_MONOMIAL, divexact(K(1), dim_j3.sum_part))
    
    println("\n LHS (Summed over x): ", total_sum)
    println(" RHS (1 / [2j3+1]): ", theoretical_rhs)
    
    is_proven = (total_sum == theoretical_rhs)
    println("\n Identity Proven? => ", is_proven ? "✅ YES" : "❌ NO")
end