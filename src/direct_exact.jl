# ---------------------------------------------------------------------------------
# Exact Algebraic Model (:exact) uses `Nemo.jl` to handle cyclotomic fields Q(ζ)
# ---------------------------------------------------------------------------------
# WARNING: The :exact mode relies on abstract computer algebra in the cyclotomic 
# field Q(ζ). It is designed strictly for small levels (k ≤ 100) for rigorous 
# mathematical proofs. For large k, the polynomial arithmetic will trigger 
# extreme memory and time complexity limits. Use :cyclo or :numeric instead.
# ---------------------------------------------------------------------------------

# Define the exact result container
struct ExactResult
    k::Int
    pref_sq::nf_elem  # Square of the prefactor (Δ)
    alt_sum::nf_elem   # The exact hypergeometric sum
end

#----- Construct exact model ------- 

struct ExactSU2kModel
    k::Int
    K::AnticNumberField            
    z::nf_elem                     
    q_facts::Vector{nf_elem} 
    q_ints::Vector{nf_elem}       
end

# LRU Caches
# const Q6J_EXACT_CACHE   = LRU{Tuple{NTuple{6, Float64}, Int}, ExactResult}(maxsize=5000)
const EXACT_MODEL_CACHE = LRU{Int, ExactSU2kModel}(maxsize=50)

function ExactSU2kModel(k::Int)
    @warn "Building ExactSU2kModel for k=$k. This is highly memory-intensive for k > 100." 
    get!(EXACT_MODEL_CACHE, k) do
        N = k + 2
        # The cyclotomic field of order 2N
        K, z = cyclotomic_field(2N, "ζ") 
        
        q_facts = Vector{nf_elem}(undef, k + 3)
        q_ints  = Vector{nf_elem}(undef, k + 3)
        
        q_facts[1] = K(1) 
        q_ints[1]  = K(0) # [0]_q = 0
        
        z_inv = inv(z)
        den_inv = inv(z - z_inv)
        
        z_n = z
        z_inv_n = z_inv
        
        @inbounds for n in 1:(k+2)
            q_int = (z_n - z_inv_n) * den_inv
            q_ints[n+1] = q_int
            q_facts[n+1] = q_facts[n] * q_int
            
            z_n *= z
            z_inv_n *= z_inv
        end
        return ExactSU2kModel(k, K, z, q_facts, q_ints)
    end
end

# ------- Exact cyclotomic field computations -------

"""
    qinteger(n, k; mode=:exact)
Returns the exact algebraic representation of [n]_q in Q(ζ).
"""
function qinteger(n::Int, k::Int; mode=:exact)
    model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
    return model.q_ints[n+1]
end

# get num and den of triangle coefficients
@inline function qΔ2_exact_numden(model::ExactSU2kModel, j1::Real, j2::Real, j3::Real)
    a = Int(j1 + j2 - j3); b = Int(j1 - j2 + j3)
    c = Int(-j1 + j2 + j3); d = Int(j1 + j2 + j3)
    
    num = model.q_facts[a+1] * model.q_facts[b+1] * model.q_facts[c+1]
    den = model.q_facts[d+2]
    return num, den
end

@inline function qtricoeff2_exact(model::ExactSU2kModel, j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real)
    n1, d1 = qΔ2_exact_numden(model, j1, j2, j3)
    n2, d2 = qΔ2_exact_numden(model, j1, j5, j6)
    n3, d3 = qΔ2_exact_numden(model, j2, j4, j6)
    n4, d4 = qΔ2_exact_numden(model, j3, j4, j5)
    
    total_num = n1 * n2 * n3 * n4
    total_den = d1 * d2 * d3 * d4
    
    return total_num * inv(total_den) 
end

function q6jseries_exact(model::ExactSU2kModel, j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real)
    # bounds 
    α1 = Int(j1 + j2 + j3) 
    α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6) 
    α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5) 
    β2 = Int(j1 + j3 + j4 + j6)
    β3 = Int(j2 + j3 + j5 + j6)
    
    z_min = max(α1, α2, α3, α4)
    z_max = min(β1, β2, β3, model.k)
    
    z_min > z_max && return model.K(0)
    
    # Calculate the initial term in the sum explicitly
    num_term = model.q_facts[z_min+2]
    den_term = model.q_facts[z_min-α1+1] * model.q_facts[z_min-α2+1] * model.q_facts[z_min-α3+1] * model.q_facts[z_min-α4+1] * model.q_facts[β1-z_min+1] * model.q_facts[β2-z_min+1] * model.q_facts[β3-z_min+1]
    
    term = num_term * inv(den_term)
    if isodd(z_min)
        term = -term
    end
    
    sum_cf = term
    
    # Iteratively update the term using the exact algebraic ratio R_z
    @inbounds for z in z_min : (z_max - 1)
        #Direct lookup of quantum integers entirely bypasses GCD inverses!
        num_ratio = -(model.q_ints[z+3]) * (model.q_ints[β1-z+1]) * (model.q_ints[β2-z+1]) * (model.q_ints[β3-z+1])
        den_ratio = (model.q_ints[z-α1+2]) * (model.q_ints[z-α2+2]) * (model.q_ints[z-α3+2]) * (model.q_ints[z-α4+2])
                    
        term = term * (num_ratio * inv(den_ratio))
        sum_cf += term
    end
    
    return sum_cf
end

function q6j_exact(model::ExactSU2kModel, j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real)
    Tc2 = qtricoeff2_exact(model, j1, j2, j3, j4, j5, j6)
    Sum_cf = q6jseries_exact(model, j1, j2, j3, j4, j5, j6)
    return ExactResult(model.k, Tc2, Sum_cf)
end

export q6j_exact, q3j_exact

# ============================================================
# Exact Quantum 3j Symbol
# ============================================================

@inline function q3j_pref_sq_exact(model::ExactSU2kModel, j1::Real, j2::Real, j3::Real, m1::Real, m2::Real)
    num, den = qΔ2_exact_numden(model, j1, j2, j3) # FIXED OUTDATED NAME
    facts = model.q_facts[Int(j1+m1)+1] * model.q_facts[Int(j1-m1)+1] * model.q_facts[Int(j2+m2)+1] * model.q_facts[Int(j2-m2)+1] * model.q_facts[Int(j3-m1-m2)+1] * model.q_facts[Int(j3+m1+m2)+1]
    
    return (num * facts) * inv(den)
end

function q3jseries_exact(model::ExactSU2kModel, j1::Real, j2::Real, j3::Real, m1::Real, m2::Real)
    α1 = Int(j3 - j2 + m1) 
    α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 - j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    z_min = max(-α1, -α2, 0)
    z_max = min(β1, β2, β3, model.k)
    
    z_min > z_max && return model.K(0)
    
    term = inv(
        model.q_facts[z_min+1] * model.q_facts[α1+z_min+1] * model.q_facts[α2+z_min+1] * model.q_facts[β1-z_min+1] * model.q_facts[β2-z_min+1] * model.q_facts[β3-z_min+1]
    )
    
    phase_offset = α1 - α2
    if isodd(z_min + phase_offset)
        term = -term
    end
    
    sum_cf = term

    @inbounds for z in z_min : (z_max - 1)
        #Direct lookup of quantum integers
        num_ratio = -(model.q_ints[β1-z+1]) * (model.q_ints[β2-z+1]) * (model.q_ints[β3-z+1])
        den_ratio = (model.q_ints[z+2]) * (model.q_ints[α1+z+2]) * (model.q_ints[α2+z+2])
                    
        term = term * (num_ratio * inv(den_ratio))
        sum_cf += term
    end
    
    return sum_cf
end

function q3j_exact(model::ExactSU2kModel, j1::Real, j2::Real, j3::Real, m1::Real, m2::Real)
    pref_sq = q3j_pref_sq_exact(model, j1, j2, j3, m1, m2)
    Sum_cf = q3jseries_exact(model, j1, j2, j3, m1, m2)
    return ExactResult(model.k, pref_sq, Sum_cf)
end





# ==============================================================================
# Exact Field Evaluation Engine
# ==============================================================================

"""
    horner_eval(poly_elem, z::Complex{BigFloat})
Evaluates a Nemo cyclotomic polynomial element at a complex point `z` using a 
zero-allocation Horner's method.
"""
function horner_eval(poly_elem, z::Complex{BigFloat})
    # Fast path for pure rational numbers
    if poly_elem isa Nemo.QQFieldElem || poly_elem isa Rational
        num = BigFloat(numerator(poly_elem))
        den = BigFloat(denominator(poly_elem))
        return Complex{BigFloat}(num / den)
    end
    
    deg = Nemo.degree(parent(poly_elem))
    res = zero(Complex{BigFloat})
    
    # Horner's method safely extracting QQFieldElem coefficients
    for i in (deg - 1):-1:0
        c = Nemo.coeff(poly_elem, i)
        
        c_num = BigFloat(BigInt(numerator(c)))
        c_den = BigFloat(BigInt(denominator(c)))
        
        res = res * z + (c_num / c_den)
    end
    
    return res
end

"""
    evaluate_exact(res::ExactResult, [T=Complex{BigFloat}]; prec=256)
Projects a full ExactResult into a floating-point number.
"""
function evaluate_exact(res::ExactResult, ::Type{T}=Complex{BigFloat}; prec=256) where {T}
    val = setprecision(BigFloat, prec) do
        # cispi(x) computes exp(i * pi * x)
        target_z = cispi(one(BigFloat) / (res.k + 2))
        
        val_sum = horner_eval(res.alt_sum, target_z)
        val_pref_sq = horner_eval(res.pref_sq, target_z)
        
        val_pref = sqrt(abs(val_pref_sq))
        
        return val_pref * val_sum
    end
    
    return T <: Real ? T(real(val)) : T(val)
end