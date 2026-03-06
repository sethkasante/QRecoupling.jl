# src/ExactAlgebra.jl


# Global Polynomial Cache (Independent of k)
# ============================================================

# Create the generic integer polynomial ring exactly once globally
const ZX_RING, ZX_VAR = polynomial_ring(ZZ, "x")

# Cache the abstract cyclotomic polynomials Φ_d(x)
const CYCLO_POLY_CACHE = LRU{Int, elem_type(ZX_RING)}(maxsize = 1024)


"""
    get_cyclo_poly(d::Int)

Retrieves the abstract integer cyclotomic polynomial Φ_d(x) from the global cache,
or computes and caches it if it doesn't exist.
"""
@inline function get_cyclo_poly(d::Int)
    return get!(CYCLO_POLY_CACHE, d) do
        cyclotomic(d, ZX_VAR)
    end
end


# ============================================================
# Exact Algebraic Model using (Nemo.jl)
# ============================================================

function ExactSU2kModel(k::Int)
    N = k + 2
    # z is a primitive 2N-th root of unity: z = exp(i * pi / N)
    K, z = cyclotomic_field(2N, "ζ") 
    
    # We cache evaluations of Φ_d at z^2 (which represents q)
    max_d = 2 * N 
    
    # Pre-allocate a dense vector for O(1) field evaluations during Racah sums
    Phi_eval = Vector{nf_elem}(undef, max_d)
    
    z_sq = z^2
    @inbounds for d in 1:max_d
        if d == 1
            Phi_eval[d] = z_sq - K(1) # Φ_1(q) = q - 1
        else
            # 1. Grab the abstract polynomial from the GLOBAL cache
            poly = get_cyclo_poly(d)
            # 2. Evaluate it in the LOCAL number field
            Phi_eval[d] = evaluate(poly, z_sq)  
        end
    end

    # Return the model with its local field evaluation cache
    return ExactSU2kModel(k, K, z, Phi_eval)
end

# function ExactSU2kModel(k::Int)
#     N = k + 2
#     # z is a primitive 2N-th root of unity: z = exp(i * pi / N)
#     K, z = cyclotomic_field(2N, "ζ") 
    
#     Zx, x = polynomial_ring(ZZ, "x") 
    
#     # Cache evaluations of Φ_d at z^2 (which represents q)
#     # Replaced Dict with a dense Vector for O(1) instant access
#     max_d = 2 * N 
#     Phi_cache = Vector{nf_elem}(undef, max_d)
    
#     z_sq = z^2
#     @inbounds for d in 1:max_d
#         if d == 1
#             Phi_cache[d] = z_sq - K(1) # Φ_1(q) = q - 1
#         else
#             poly = cyclotomic(d, x)
#             Phi_cache[d] = evaluate(poly, z_sq)  
#         end
#     end

#     return ExactSU2kModel(k, K, z, Phi_cache)
# end

# ============================================================
# Field Evaluation in Cyclotomic Fields 
# ============================================================

"""
Maps a purely symbolic CycloMonomial exactly into the Cyclotomic number field.
"""
function evaluate_cyclofield(M::CycloMonomial, model::ExactSU2kModel)
    K = model.K
    res = K(M.sign) 

    # Safe exponentiation for field generators
    if M.z_pow != 0
        base_z = M.z_pow > 0 ? model.z : inv(model.z)
        res *= base_z^abs(M.z_pow)
    end
    
    # Iterate over the Vector indices (which represent the cyclotomic base d)
    @inbounds for d in 1:length(M.exps)
        e = M.exps[d]
        if e != 0
            term = model.Phi_eval[d]
            # If e < 0 and term == 0, Nemo will throw a DivideError. 
            # However, for an admissible 6j symbol, denominators never exceed k+1, 
            # so we never divide by Φ_{2N}(z^2) == 0.
            res *= e > 0 ? term^e : inv(term)^abs(e)
        end
    end
    return res
end

# ============================================================
# Exact Cyclo Field Evaluations 
# ============================================================

# Evaluate the triangle coeffs 
function qtricoeff2_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    Δ2 = qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    return evaluate_cyclofield(Δ2, model)
end

# Evaluate the racah sum arithmetics 
function q6jseries_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    # Removed artificial `model.k` cap. The field arithmetic handles [k+2]_q = 0 exactly.
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3) 
    
    sum_cf = model.K(0)

    @inbounds for z in zrange
        val = q6jsummand_symb(z, α1, α2, α3, α4, β1, β2, β3) 
        sum_cf += evaluate_cyclofield(val, model)
    end
    return sum_cf
end

# Function for quantum6j symbols 
function qracah6j_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return ExactResult(model.k, model.K(0), model.K(0))
    end
    
    Tc2 = qtricoeff2_exact(model, j1, j2, j3, j4, j5, j6)
    Sum_cf = q6jseries_exact(model, j1, j2, j3, j4, j5, j6)
    
    return ExactResult(model.k, Tc2, Sum_cf)
end

# TODO Fulfilled: Standalone Wrapper
"""
    qracah6j_exact(j1, j2, j3, j4, j5, j6, k::Int) -> ExactResult

Creates a temporary ExactSU2kModel and evaluates the exact algebraic 6j symbol.
"""
function qracah6j_exact(j1::Number, j2, j3, j4, j5, j6, k::Int)
    model = ExactSU2kModel(k)
    return qracah6j_exact(model, j1, j2, j3, j4, j5, j6)
end

# ============================================================
# Float Projection (Horner's Method)
# ============================================================

"""
    evaluate_exact(res::ExactResult, prec::Int=256) -> BigFloat

Converts an exact algebraic result into a high-precision float.
"""
function evaluate_exact(res::ExactResult, prec::Int=256)
    # If admissibility failed, it's exactly 0
    if res.pref_sq == 0 && res.sum_cf == 0
        return zero(BigFloat)
    end
    
    setprecision(BigFloat, prec) do
        N = res.k + 2
        # Target z = exp(iπ/N)
        target_z = cispi(big"1.0" / N)
        
        # 1. Evaluate S (Racah Sum) via Horner
        sum_bf = horner_eval(res.sum_cf, target_z)
        
        # 2. Evaluate Δ² (Prefactor squared) via Horner
        pref2_bf = horner_eval(res.pref_sq, target_z)
        
        # 3. Final Assembly
        return real(sqrt(pref2_bf) * sum_bf)
    end
end

"""
    horner_eval(ev::nf_elem, root_val::Complex{BigFloat}) -> Complex{BigFloat}

Evaluates a Nemo number field element at a specific complex root using Horner's scheme.
"""
function horner_eval(ev::nf_elem, root_val::Complex{BigFloat})
    d = degree(parent(ev)) - 1
    d < 0 && return zero(Complex{BigFloat})
    
    res = Complex{BigFloat}(coeff(ev, d))
    for i in (d-1):-1:0
        c = coeff(ev, i)
        # Use Nemo.numerator and Nemo.denominator strictly to avoid Base conflicts
        val = BigFloat(Nemo.numerator(c)) / BigFloat(Nemo.denominator(c))
        res = res * root_val + val
    end
    return res
end