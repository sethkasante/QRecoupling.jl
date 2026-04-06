
# ---------------------------------------------------------------------------------
# Exact Algebraic Model (:exact) uses `Nemo.jl` to handle cyclotomic fields Q(ζ)
# Exact SU(2)k eager evaluation
# Performs immediate polynomial arithmetic in the Cyclotomic Field Q(ζ).
# ---------------------------------------------------------------------------------
# WARNING: The :exact mode relies on abstract computer algebra in the cyclotomic 
# field Q(ζ). It is designed strictly for small levels (k ≤ 100) for rigorous 
# mathematical proofs. For large k, the polynomial arithmetic will trigger 
# extreme memory and time complexity limits. 
# ---------------------------------------------------------------------------------





#  ---  Result Container & Printing ----

"""
    ExactResult
The exact topological invariant computed via eager dense polynomial arithmetic.
- `radical_sq`: The exact square of the prefactor (Δ²). 
- `factor_sum`: The evaluated hypergeometric alternating sum.
"""
struct ExactResult
    k::Int                
    radical_sq::nf_elem   # Exactly Δ²
    factor_sum::nf_elem   # The series sum
end

function _print_truncated_nemo(io::IO, x::nf_elem)
    s = string(x)
    if length(s) > 120
        print(io, s[1:60], " ... ", s[end-40:end])
    else
        print(io, s)
    end
end

function Base.show(io::IO, res::ExactResult)
    k_sub = to_subscript(res.k)
    println(io, "Exact SU(2)", k_sub, " Symbol:")
    K = parent(res.factor_sum)
    
    if res.radical_sq == K(1)
        print(io, "  Value: ")
        _print_truncated_nemo(io, res.factor_sum)
    else
        println(io, "  Value: √(A) * B")
        print(io, "  A (Δ²): ")
        _print_truncated_nemo(io, res.radical_sq)
        print(io, "\n  B (Sum): ")
        _print_truncated_nemo(io, res.factor_sum)
    end
end


# ---- Build exact model (Factorial Cache) ----

struct ExactSU2kModel
    k::Int
    K::AnticNumberField            
    z::nf_elem                     
    q_facts::Vector{nf_elem} 
    q_ints::Vector{nf_elem}       
end

const EXACT_MODEL_CACHE = LRU{Int, ExactSU2kModel}(maxsize=150)

function ExactSU2kModel(k::Int)
    # @warn "Building ExactSU2kModel for SU(2)_{$k}. Dense polynomial caching can trigger memory exhaustion for k > 100." maxlog=1
    get!(EXACT_MODEL_CACHE, k) do
        h = k + 2
        # Q(ζ) where ζ = exp(i2π/2h)
        K, z = cyclotomic_field(2h, "ζ") 
        
        q_facts = Vector{nf_elem}(undef, k + 4)
        q_ints  = Vector{nf_elem}(undef, k + 4)
        
        q_facts[1] = K(1) # [0]! = 1
        q_ints[1]  = K(0) # [0] = 0
        
        z_inv = inv(z)
        den_inv = inv(z - z_inv)
        
        @inbounds for n in 1:(k+3)
            q_int = (z^n - z_inv^n) * den_inv
            q_ints[n+1] = q_int
            q_facts[n+1] = q_facts[n] * q_int
        end
        return ExactSU2kModel(k, K, z, q_facts, q_ints)
    end
end



# ---- Algebraic Symbol Builders ----- 

# get num and den of triangle coefficients
@inline function _qΔ2_exact_numden(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin)
    a, b, c = Int(j1+j2-j3), Int(j1-j2+j3), Int(-j1+j2+j3)
    d = Int(j1+j2+j3)
    num = model.q_facts[a+1] * model.q_facts[b+1] * model.q_facts[c+1]
    den = model.q_facts[d+2]
    return num, den
end


"""
    q6j_exact(j1, j2, j3, j4, j5, j6, k)
Computes the exact SU(2)k Racah symbol using the memory-optimized iterative eager model.
"""
function q6j_exact(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int)
    !qδtet(j1, j2, j3, j4, j5, j6, k) && return ExactResult(k, ExactSU2kModel(k).K(1), ExactSU2kModel(k).K(0))
    
    model = ExactSU2kModel(k)
    
    # prefactor (radical squared)
    n1, d1 = _qΔ2_exact_numden(model, j1, j2, j3)
    n2, d2 = _qΔ2_exact_numden(model, j1, j5, j6)
    n3, d3 = _qΔ2_exact_numden(model, j2, j4, j6)
    n4, d4 = _qΔ2_exact_numden(model, j3, j4, j5)
    radical_sq = (n1 * n2 * n3 * n4) * inv(d1 * d2 * d3 * d4)
    

    α1 = Int(j1+j2+j3); α2 = Int(j1+j5+j6); α3 = Int(j2+j4+j6); α4 = Int(j3+j4+j5)
    β1 = Int(j1+j2+j4+j5); β2 = Int(j1+j3+j4+j6); β3 = Int(j2+j3+j5+j6)
    
    z_min = max(α1, α2, α3, α4)
    z_max = min(β1, β2, β3, model.k)
    
    z_min > z_max && return ExactResult(k, radical_sq, model.K(0))
    
    #intial term
    num_term = model.q_facts[z_min+2]
    den_term = model.q_facts[z_min-α1+1] * model.q_facts[z_min-α2+1] * model.q_facts[z_min-α3+1] * model.q_facts[z_min-α4+1] * model.q_facts[β1-z_min+1] * model.q_facts[β2-z_min+1] * model.q_facts[β3-z_min+1]
    
    term = num_term * inv(den_term)
    if isodd(z_min)
        term = -term
    end
    sum_val = term
    
    # iterative loop
    @inbounds for z in z_min : (z_max - 1)
        num_ratio = -(model.q_ints[z+3]) * (model.q_ints[β1-z+1]) * (model.q_ints[β2-z+1]) * (model.q_ints[β3-z+1])
        den_ratio = (model.q_ints[z-α1+2]) * (model.q_ints[z-α2+2]) * (model.q_ints[z-α3+2]) * (model.q_ints[z-α4+2])
                    
        term = term * (num_ratio * inv(den_ratio))
        sum_val += term
    end
    
    return ExactResult(k, radical_sq, sum_val)
end

"""
    q3j_exact(j1, j2, j3, m1, m2, m3, k)
Computes the exact SU(2)k Wigner 3j symbol using the memory-optimized iterative eager model.
"""
function q3j_exact(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin, k::Int)
    (!qδ(j1, j2, j3, k) || m1+m2+m3 != 0) && return ExactResult(k, ExactSU2kModel(k).K(1), ExactSU2kModel(k).K(0))
    
    model = ExactSU2kModel(k)
    
    # prefactor (radical squared)
    num, den = _qΔ2_exact_numden(model, j1, j2, j3)
    facts = model.q_facts[Int(j1+m1)+1] * model.q_facts[Int(j1-m1)+1] * model.q_facts[Int(j2+m2)+1] * model.q_facts[Int(j2-m2)+1] * model.q_facts[Int(j3+m3)+1] * model.q_facts[Int(j3-m3)+1]
    triangle_sq = (num * facts) * inv(den)
    
    α1 = Int(j3 - j2 + m1); α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 - j3); β2 = Int(j1 - m1); β3 = Int(j2 + m2)
    
    z_min = max(0, -α1, -α2)
    z_max = min(β1, β2, β3, model.k)
    
    z_min > z_max && return ExactResult(k, triangle_sq, model.K(0))
    
    # base term
    term = inv(model.q_facts[z_min+1] * model.q_facts[α1+z_min+1] * model.q_facts[α2+z_min+1] * model.q_facts[β1-z_min+1] * model.q_facts[β2-z_min+1] * model.q_facts[β3-z_min+1])
    
    phase_offset = α1 - α2
    if isodd(z_min + phase_offset)
        term = -term
    end
    sum_val = term

    # loop sum
    @inbounds for z in z_min : (z_max - 1)
        num_ratio = -(model.q_ints[β1-z+1]) * (model.q_ints[β2-z+1]) * (model.q_ints[β3-z+1])
        den_ratio = (model.q_ints[z+2]) * (model.q_ints[α1+z+2]) * (model.q_ints[α2+z+2])
                    
        term = term * (num_ratio * inv(den_ratio))
        sum_val += term
    end
    
    return ExactResult(k, triangle_sq, sum_val)
end


#  --- Fsymbol ----- 

function fsymbol_exact(j1, j2, j3, j4, j5, j6, k)
    res = q6j_exact(j1, j2, j3, j4, j5, j6, k)
    model = ExactSU2kModel(k)
    # scale by sqrt([2j3+1][2j6+1]) and phase
    dim_sq = model.q_ints[Int(2j3+1)+1] * model.q_ints[Int(2j6+1)+1]
    phase = iseven(Int(j1+j2+j4+j5)) ? 1 : -1
    return ExactResult(k, res.radical_sq * dim_sq, phase * res.factor_sum)
end



#  ----  numerical evaluation to Float/Complex  ---- 


"""
    evaluate_exact(res::ExactResult, [T=ComplexF64])
Projects the cyclotomic exact result into a complex/real number using Horner's method.
Inherits ambient `BigFloat` precision automatically for intermediate steps.
"""
function evaluate_exact(res::ExactResult, ::Type{T}=ComplexF64) where T
    h = res.k + 2
    target_z = cispi(one(BigFloat) / h)
    
    # Helper: Horner evaluation with SECURE BigFloat casting
    function _horner(poly, z)
        deg = degree(parent(poly))
        val = Complex{BigFloat}(0)
        for i in (deg-1):-1:0
            c = coeff(poly, i)
            # Bypass Float64 truncation by casting Ints natively
            c_bf = BigFloat(numerator(c)) / BigFloat(denominator(c))
            val = val * z + c_bf
        end
        return val
    end
    
    B = _horner(res.factor_sum, target_z)
    A_sq = _horner(res.radical_sq, target_z)
    
    A = sqrt(max(zero(BigFloat), real(A_sq)))
    val = A * B
    
    return T <: Real ? T(real(val)) : T(val)
end

