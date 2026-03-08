#exactalgebra.jl

struct ExactSU2kModel
    k::Int
    K::AnticNumberField            
    z::nf_elem                     
    q_facts::Vector{nf_elem}       
end

# ============================================================
# Exact Algebraic Model using (Nemo.jl)
# ============================================================

function ExactSU2kModel(k::Int)
    N = k + 2
    K, z = cyclotomic_field(2N, "ζ") 
    
    q_facts = Vector{nf_elem}(undef, k + 3)
    q_facts[1] = K(1) 
    
    z_inv = inv(z)
    den_inv = inv(z - z_inv)
    z_n = z
    z_inv_n = z_inv
    
    @inbounds for n in 1:(k+2)
        q_int = (z_n - z_inv_n) * den_inv
        q_facts[n+1] = q_facts[n] * q_int
        z_n *= z
        z_inv_n *= z_inv
    end
    return ExactSU2kModel(k, K, z, q_facts)
end

# ============================================================
# Exact Cyclo Field Evaluations 
# ============================================================

"""
    qinteger(n, k; mode=:exact)
Returns the exact algebraic representation of [n]_q in Q(ζ).
"""
function qinteger(n::Int, k::Int; mode=:exact)
    model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
    #enforce [1] = 1 = [0] 
    if n == 0 || n==1 return model.K(1) end 
    return model.q_facts[n+1] * inv(model.q_facts[n])
end

@inline function q_delta2_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin)
    a = Int(j1 + j2 - j3); b = Int(j1 - j2 + j3)
    c = Int(-j1 + j2 + j3); d = Int(j1 + j2 + j3)
    
    num = model.q_facts[a+1] * model.q_facts[b+1] * model.q_facts[c+1]
    den = model.q_facts[d+2]
    return num * inv(den)
end

@inline function qtricoeff2_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    return q_delta2_exact(model, j1, j2, j3) * q_delta2_exact(model, j1, j5, j6) * q_delta2_exact(model, j2, j4, j6) * q_delta2_exact(model, j3, j4, j5)
end


function q6jseries_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3, model.k) 
    sum_cf = model.K(0)

    @inbounds for z in zrange
        num = model.q_facts[z+2]
        den = model.q_facts[z-α1+1] * model.q_facts[z-α2+1] * model.q_facts[z-α3+1] * model.q_facts[z-α4+1] * model.q_facts[β1-z+1] * model.q_facts[β2-z+1] * model.q_facts[β3-z+1]
        
        term = num * inv(den)
        sum_cf += iseven(z) ? term : -term
    end
    return sum_cf
end

function qracah6j_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    Tc2 = qtricoeff2_exact(model, j1, j2, j3, j4, j5, j6)
    Sum_cf = q6jseries_exact(model, j1, j2, j3, j4, j5, j6)
    return ExactResult(model.k, Tc2, Sum_cf)
end

# ============================================================
# Exact Quantum 3j Symbol
# ============================================================

@inline function q3j_pref_sq_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin)
    delta2 = q_delta2_exact(model, j1, j2, j3)
    facts = model.q_facts[Int(j1+m1)+1] * model.q_facts[Int(j1-m1)+1] * model.q_facts[Int(j2+m2)+1] * model.q_facts[Int(j2-m2)+1] * model.q_facts[Int(j3-m1-m2)+1] * model.q_facts[Int(j3+m1+m2)+1]
    return delta2 * facts
end


function q3jseries_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin)
    α1 = Int(j3 - j2 + m1) 
    α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 - j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    # Corrected z_min bound to prevent negative factorials
    zrange = max(-α1, -α2, 0):min(β1, β2, β3, model.k)
    sum_cf = model.K(0)
    
    # Pre-calculate the static part of the phase exponent
    phase_offset = α1 - α2

    @inbounds for z in zrange
        den = model.q_facts[z+1] * model.q_facts[α1+z+1] * model.q_facts[α2+z+1] * model.q_facts[β1-z+1] * model.q_facts[β2-z+1] * model.q_facts[β3-z+1]
        
        term = inv(den)
        # Phase is (-1)^{z + α1 - α2}
        sum_cf += isodd(z + phase_offset) ? -term : term
    end
    
    return sum_cf
end

function qracah3j_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin)
    pref_sq = q3j_pref_sq_exact(model, j1, j2, j3, m1, m2)
    Sum_cf = q3jseries_exact(model, j1, j2, j3, m1, m2)
    return ExactResult(model.k, pref_sq, Sum_cf)
end



# ==============================================================================
# Exact Field Evaluation Engine
# ==============================================================================

"""
    horner_eval(poly_elem, z::Complex{BigFloat})

Evaluates a Nemo cyclotomic polynomial element at a complex point `z` 
using a zero-allocation Horner's method.
"""
function horner_eval(poly_elem, z::Complex{BigFloat})
    # Fast path for pure rational numbers (QQFieldElem)
    if poly_elem isa Nemo.QQFieldElem
        num = BigFloat(numerator(poly_elem))
        den = BigFloat(denominator(poly_elem))
        return Complex{BigFloat}(num / den)
    end
    
    deg = Nemo.degree(parent(poly_elem))
    res = zero(Complex{BigFloat})
    
    # Horner's method: Evaluate without allocating an array
    for i in (deg - 1):-1:0
        c = Nemo.coeff(poly_elem, i)
        
        # Safely extract Nemo coefficients to BigFloat without Float64 truncation
        c_val = BigFloat(numerator(c)) / BigFloat(denominator(c))
        
        res = res * z + c_val
    end
    
    return res
end

"""
    evaluate_exact(ev::ExactValue, [T=Complex{BigFloat}]; prec=256)

Projects a self-contained ExactValue into a floating-point number.
"""
function evaluate_exact(ev::ExactValue, ::Type{T}=Complex{BigFloat}; prec=256) where {T}
    val = setprecision(BigFloat, prec) do
        target_z = cispi(big(1.0) / (ev.k + 2))
        return horner_eval(ev.val, target_z)
    end
    
    if T <: Real
        return T(real(val))
    else
        return T(val)
    end
end


"""
    evaluate_exact(res::ExactResult, ::Type{T}=Complex{BigFloat}; prec=256)

Projects a full ExactResult (6j or 3j symbol) into a floating-point number.
"""
function evaluate_exact(res::ExactResult, ::Type{T}=Complex{BigFloat}; prec=256) where {T}
    # Calculate in high precision
    val = setprecision(BigFloat, prec) do
        target_z = cispi(big(1.0) / (res.k + 2))
        
        val_sum = horner_eval(res.sum_cf, target_z)
        val_pref_sq = horner_eval(res.pref_sq, target_z)
        
        # In SU(2)_k, pref_sq is mathematically real and non-negative
        val_pref = sqrt(abs(val_pref_sq))
        
        return val_pref * val_sum
    end
    
    # Clean casting
    if T <: Real
        return T(real(val))
    else
        return T(val)
    end
end



# ============================================================
# Float Projection (Horner's Method)
# ============================================================

# """
#     evaluate_exact(res::ExactResult, [T=Complex{BigFloat}]; prec=256)
# Projects an exact algebraic result into a floating-point number.
# Defaults to Complex{BigFloat}. Can be cast to Float64, ComplexF64, etc.
# """
# function evaluate_exact(res::ExactResult, ::Type{T}=Complex{BigFloat}; prec=256) where {T}
#     # Calculate in high precision
#     val = setprecision(BigFloat, prec) do
#         target_z = cis(BIG_PI / (res.k + 2))
#         val_sum = horner_eval(res.sum_cf, target_z)
#         val_pref_sq = horner_eval(res.pref_sq, target_z)
        
#         val_pref = sqrt(abs(val_pref_sq))
#         return val_pref * val_sum
#     end
    
#     # Intelligently cast based on the requested type T
#     if T <: Real
#         # For Real types (e.g., Float64), we assume the imaginary part is numerical noise
#         return T(real(val))
#     else
#         # For Complex types, we return the full value
#         return T(val)
#     end
# end




