
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
    get!(EXACT_MODEL_CACHE, k) do
        h = k + 2
        # Q(ζ) where ζ = exp(iπ/h)
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

@inline function qΔ2_exact(model::ExactSU2kModel, j1, j2, j3)
    a, b, c = Int(j1+j2-j3), Int(j1-j2+j3), Int(-j1+j2+j3)
    d = Int(j1+j2+j3)
    # [a]![b]![c]! / [d+1]!
    return (model.q_facts[a+1] * model.q_facts[b+1] * model.q_facts[c+1]) * inv(model.q_facts[d+2])
end

"""
    q6j_exact(j1, j2, j3, j4, j5, j6, k)
Computes the exact SU(2)k Racah symbol using the eager model.
"""
function q6j_exact(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int)
    !qδtet(j1, j2, j3, j4, j5, j6,k) && return ExactResult(k, ExactSU2kModel(k).K(1), ExactSU2kModel(k).K(0))
    
    model = ExactSU2kModel(k)
    # Δ² = Δ123² * Δ156² * Δ246² * Δ345²
    radical_sq = qΔ2_exact(model, j1, j2, j3) * qΔ2_exact(model, j1, j5, j6) *
                 qΔ2_exact(model, j2, j4, j6) * qΔ2_exact(model, j3, j4, j5)
    
    α = (Int(j1+j2+j3), Int(j1+j5+j6), Int(j2+j4+j6), Int(j3+j4+j5))
    β = (Int(j1+j2+j4+j5), Int(j1+j3+j4+j6), Int(j2+j3+j5+j6))
    z_min, z_max = max(α...), min(β...)
    
    sum_val = model.K(0)
    for z in z_min:z_max
        num = model.q_facts[z+2]
        den = model.q_facts[z-α[1]+1] * model.q_facts[z-α[2]+1] * model.q_facts[z-α[3]+1] * 
                    model.q_facts[z-α[4]+1] * model.q_facts[β[1]-z+1] * model.q_facts[β[2]-z+1] * 
                        model.q_facts[β[3]-z+1]
        term = num * inv(den)
        sum_val = iseven(z) ? (sum_val + term) : (sum_val - term)
    end
    
    return ExactResult(k, radical_sq, sum_val)
end

"""
    q3j_exact(j1, j2, j3, m1, m2, m3, k)
Computes the exact SU(2)k Wigner 3j symbol.
"""
function q3j_exact(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin, k::Int)
    (!qδ(j1, j2, j3,k) || m1+m2+m3 != 0) && return ExactResult(k, ExactSU2kModel(k).K(1), ExactSU2kModel(k).K(0))
    
    model = ExactSU2kModel(k)
    # radical_sq = Δ123² * Π [j±m]!
    triangle_sq = qΔ2_exact(model, j1, j2, j3)
    m_facts = model.q_facts[Int(j1+m1)+1] * model.q_facts[Int(j1-m1)+1] *
              model.q_facts[Int(j2+m2)+1] * model.q_facts[Int(j2-m2)+1] *
              model.q_facts[Int(j3+m3)+1] * model.q_facts[Int(j3-m3)+1]
    
    α = (Int(j3-j2+m1), Int(j3-j1-m2))
    β = (Int(j1+j2-j3), Int(j1-m1), Int(j2+m2))
    z_min, z_max = max(0, -α[1], -α[2]), min(β[1], β[2], β[3])
    
    sum_val = model.K(0)
    for z in z_min:z_max
        den = model.q_facts[z+1] * model.q_facts[α[1]+z+1] * model.q_facts[α[2]+z+1] *
              model.q_facts[β[1]-z+1] * model.q_facts[β[2]-z+1] * model.q_facts[β[3]-z+1]
        term = inv(den)
        # Sign: (-1)^{z + α1 - α2}
        sum_val = iseven(z + α[1] - α[2]) ? (sum_val + term) : (sum_val - term)
    end
    
    return ExactResult(k, triangle_sq * m_facts, sum_val)
end


#  --- TQFT Kernels (F & G Symbols) ----- 

function fsymbol_exact(j1, j2, j3, j4, j5, j6, k)
    res = q6j_exact(j1, j2, j3, j4, j5, j6, k)
    model = ExactSU2kModel(k)
    # Unitary F-symbol scales by sqrt([2j3+1][2j6+1]) and phase
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

