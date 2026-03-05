
#ExactAlgebra.jl

# ============================================================
# Exact Algebraic Model using (Nemo.jl)
# ============================================================

function ExactSU2kModel(k::Int)
    N = k + 2
    # z is a primitive 2N-th root of unity: z = exp(i * pi / N)
    K, z = cyclotomic_field(2N, "ζ") 
    
    Zx, x = polynomial_ring(ZZ, "x") 
    
    # We cache evaluations of Φ_d at z^2 (which represents q)
    # The maximum value inside a factorial is 2k (which is < 2N)
    max_d = 2 * N 
    Phi_cache = Dict{Int, nf_elem}()
    
    z_sq = z^2
    for d in 2:max_d
        poly = cyclotomic(d, x)
        Phi_cache[d] = evaluate(poly, z_sq)  
    end

    return ExactSU2kModel(k, K, z, Phi_cache)
end


# ========================================
# Field Evaluation in Cyclotomic Fields 
# ========================================

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
    for (d, e) in M.exps
        if e != 0
            term = model.Phi_cache[d]
            res *= e > 0 ? term^e : inv(term)^(-e)
        end
    end
    return res
end


#exact cyclo field evaluations 

#evaluate the triangle coeffs 
function qtricoeff2_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    Δ2 = qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    return evaluate_cyclofield(Δ2, model::ExactSU2kModel)
end


#evaluate the racah sum arithmetics 
#evaluate the racah sum arithmetics 
function q6jseries_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)

    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    zrange = max(α1, α2, α3, α4):min(min(β1, β2, β3),model.k) #take care of [k+2]! terms!! z has to be less than k  
    sum_cf = model.K(0)

    @inbounds for z in zrange
        val = q6jsummand_symb(z, α1, α2, α3, α4, β1, β2, β3) 
        sum_cf += evaluate_cyclofield(val, model)
    end
    return sum_cf
end

# function q6jseries_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
#     series = q6jseries_symb(j1, j2, j3, j4, j5, j6)
#     sum_cf = model.K(0)
#     for s in series
#         sum_cf += evaluate_cyclofield(s, model)
#     end
#     return sum_cf
# end


#function for quantum6j symbols 
function qracah6j_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return Exact6jResult(model.k,model.K(0), model.K(0))
    end
    
    Tc2 = qtricoeff2_exact(model, j1, j2, j3, j4, j5, j6)
    Sum_cf = q6jseries_exact(model, j1, j2, j3, j4, j5, j6)
    
    return Exact6jResult(model.k,Tc2, Sum_cf) # Exact results  
end

# src/ExactAlgebra.jl

"""
    evaluate_exact(res::Exact6jResult, prec::Int=256) -> BigFloat

Converts an exact algebraic result into a high-precision float.
"""
function evaluate_exact(res::Exact6jResult, prec::Int=256)
    setprecision(BigFloat, prec) do
        N = res.k + 2
        # Target z = exp(iπ/N). Note z² = q = exp(i2π/N)
        target_z = cispi(big"1.0" / N)
        
        # 1. Evaluate S (Racah Sum) via Horner
        sum_bf = horner_eval(res.sum_cf, target_z)
        
        # 2. Evaluate Δ² (Prefactor squared) via Horner
        pref2_bf = horner_eval(res.pref_sq, target_z)
        
        # 3. Final Assembly: sqrt(Δ²) * S
        # We take real() because the final SU(2)_k symbol is strictly real
        return real(sqrt(pref2_bf) * sum_bf)
    end
end

"""
    horner_eval(ev::nf_elem, root_val::Complex{BigFloat})

Evaluates a Nemo number field element at a specific complex root using Horner's scheme.
"""
function horner_eval(ev::nf_elem, root_val::Complex{BigFloat})
    d = degree(parent(ev)) - 1
    d < 0 && return zero(Complex{BigFloat})
    
    # Start with the highest coefficient
    res = Complex{BigFloat}(coeff(ev, d))
    for i in (d-1):-1:0
        c = coeff(ev, i)
        # Convert rational coefficient to BigFloat
        val = BigFloat(numerator(c)) / BigFloat(denominator(c))
        res = res * root_val + val
    end
    return res
end



# for when numeric is needed: not necessary here but we keep it anyways to compare to numeric later. 
"""
Safely evaluates an exact number field element to a Complex{BigFloat} using Horner's method.
"""
function evaluate_horner(ev::nf_elem, root_val)
    d = degree(parent(ev)) - 1
    d < 0 && return zero(Complex{BigFloat})
    
    res = BigFloat(coeff(ev, d))
    for i in (d-1):-2:0
        c = BigFloat(coeff(ev, i))
        res = res * root_val + c
    end
    return res
end

#TODO: Write a function qracah6j_exact(j1,j2,j3,j4,j5,j6,k::Int) 