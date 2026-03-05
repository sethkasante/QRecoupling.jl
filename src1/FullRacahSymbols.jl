module QRacah

using Nemo
using LRUCache

export CycloMonomial, ExactSU2kModel, NumericSU2kModel
export qfactorial_symb, q6jseries_symb, qtricoeff2_symb, evaluate_nf, evaluate_horner
export exact_qracah6j, numeric_qracah6j, qracah6j_symb, qracah6jsq_symb

# ============================================================
# 1. Symbolic Data Structures (Quantum Primes)
# ============================================================

"""
    CycloMonomial

Represents exactly a product of cyclotomic polynomials evaluated at q:
Value = sign * z^(z_pow) * Π (Φ_d(z^2))^(exps[d])

Here, z = q^{1/2} = e^{i π / (k+2)}.
"""
struct CycloMonomial
    sign::Int
    z_pow::Int           # Handles the integer powers of z (which are fractional powers of q)
    exps::Dict{Int, Int} # Handles Π Φ_d(z^2)^{exps[d]}
end

import Base: *, /

function *(a::CycloMonomial, b::CycloMonomial)
    small, large = length(a.exps) ≤ length(b.exps) ? (a,b) : (b,a)
    exps = copy(large.exps)
    for (d, e) in small.exps
        val = get(exps, d, 0) + e
        if val == 0
            delete!(exps, d)
        else
            exps[d] = val
        end
    end
    return CycloMonomial(a.sign * b.sign, a.z_pow + b.z_pow, exps)
end

function /(a::CycloMonomial, b::CycloMonomial)
    exps = copy(a.exps)
    for (d, e) in b.exps
        val = get(exps, d, 0) - e
        if val == 0
            delete!(exps, d)
        else
            exps[d] = val
        end
    end
    return CycloMonomial(a.sign * b.sign, a.z_pow - b.z_pow, exps)
end

# ============================================================
# 2. Exact Algebraic Model (Nemo.jl)
# ============================================================

struct ExactSU2kModel
    k::Int
    K::AnticNumberField            
    z::nf_elem                     # Primitive 2N-th root of unity (z = q^{1/2})
    Phi_cache::Dict{Int, nf_elem}  # Cached evaluations of Φ_d(z^2)
end

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

# ============================================================
# 3. Numeric Model (Log-Sum-Exp Fallback)
# ============================================================

# struct NumericSU2kModel{T<:AbstractFloat}
#     k::Int
#     logqnfact::Vector{T}    
# end

# const LOGQFACT_CACHE = LRU{Tuple{Int,Int}, Vector{BigFloat}}(maxsize = 1024)

# function NumericSU2kModel(k::Int, prec=256)::NumericSU2kModel
#     tab = get!(LOGQFACT_CACHE, (k, prec)) do
#         logqnfact_table(k, prec) 
#     end
#     return NumericSU2kModel{BigFloat}(k, tab)
# end

# function logqn_table(k::Int, prec::Int=256)::Vector{BigFloat}
#     N = k + 2
#     setprecision(prec) do
#         θ = big(pi) / BigFloat(N)
#         logsinθ = log(sin(θ))

#         tab = Vector{BigFloat}(undef, N)
#         half = N ÷ 2
#         tab[1] = zero(BigFloat) 

#         @inbounds for n in 1:half
#             val = log(sin(n*θ)) - logsinθ
#             tab[N + 1 - n] = tab[n+1] = val
#         end
#         return tab
#     end
# end

# function logqnfact_table(k::Int, prec::Int=256)::Vector{BigFloat}
#     logqn = logqn_table(k, prec)
#     setprecision(prec) do
#         tab = Vector{BigFloat}(undef, k+2)
#         tab[1] = logqn[1]  
#         @inbounds for n in 2:k+2
#             tab[n] = tab[n-1] + logqn[n]
#         end
#         return tab
#     end
# end

# ============================================================
# 4. Symbolic Math Engine
# ============================================================

function qfactorial_symb(n::Int)
    if n == 0
        return CycloMonomial(1, 0, Dict{Int, Int}())
    end

    exps = Dict{Int, Int}()
    for d in 2:n
        power = div(n, d) 
        if power > 0
            exps[d] = power
        end
    end
    
    # Power of z = q^{1/2}. q^{-n(n-1)/4} = z^{-n(n-1)/2}
    z_pow = - (n * (n - 1)) ÷ 2
    
    return CycloMonomial(1, z_pow, exps)
end

function qdelta2_symb(a, b, c)
    num = qfactorial_symb(Int(a+b-c)) * qfactorial_symb(Int(a-b+c)) * qfactorial_symb(Int(-a+b+c))
    den = qfactorial_symb(Int(a+b+c+1))
    return num / den
end

function qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    return qdelta2_symb(j1, j2, j3) * qdelta2_symb(j1, j5, j6) * qdelta2_symb(j2, j4, j6) * qdelta2_symb(j3, j4, j5)
end

function q6jsummand_symb(z, α1, α2, α3, α4, β1, β2, β3)
    num = qfactorial_symb(z+1)
    den = qfactorial_symb(z-α1) * qfactorial_symb(z-α2) * qfactorial_symb(z-α3) *
          qfactorial_symb(z-α4) * qfactorial_symb(β1-z) * qfactorial_symb(β2-z) *
          qfactorial_symb(β3-z) 
    
    res = num / den
    if isodd(z)
        return CycloMonomial(-res.sign, res.z_pow, res.exps)
    else
        return res
    end
end

function q6jseries_symb(j1, j2, j3, j4, j5, j6)::Vector{CycloMonomial}
    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    S_z = CycloMonomial[]
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3)
    @inbounds for z in zrange
        push!(S_z, q6jsummand_symb(z, α1, α2, α3, α4, β1, β2, β3)) 
    end
    return S_z
end

# ============================================================
# 5. Field Evaluation & BigFloat Bridge
# ============================================================

"""
Maps a purely symbolic CycloMonomial exactly into the Nemo number field.
"""
function evaluate_nf(M::CycloMonomial, model::ExactSU2kModel)
    K = model.K
    res = K(M.sign) * model.z^M.z_pow
    for (d, e) in M.exps
        res *= (model.Phi_cache[d])^e
    end
    return res
end


function qracah6j_symb(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return zero(BigFloat)
    end
    
    # 1. Exact Summation in the Field
    S_z = q6jseries_symb(j1, j2, j3, j4, j5, j6)
    sum_nf = model.K(0)
    for M in S_z
        sum_nf += evaluate_nf(M, model)
    end
    
    # 2. Exact Prefactor ^ 2 in the Field
    Pref2_symb = qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    Pref2_nf = evaluate_nf(Pref2_symb, model)
    return Pref2_nf, sum_nf
end

function qracah6jsq_symb(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    Pref2_nf, sum_nf = qracah6j_symb(model, j1, j2, j3, j4, j5, j6)
    return Pref2_nf * sum_nf^2
    #What about the sign? 
end


# function evaluate_nf(M::CycloMonomial, model::ExactSU2kModel)
#     K = model.K
#     res = K(M.sign)
    
#     # Handle the generator z precisely (avoiding Nemo negative power errors)
#     if M.z_pow != 0
#         base_z = M.z_pow > 0 ? model.z : inv(model.z)
#         res *= base_z^abs(M.z_pow)
#     end
    
#     # Handle the cyclotomic polynomials evaluated at z^2
#     for (d, e) in M.exps
#         if e != 0
#             term = model.Phi_cache[d]
#             if e > 0
#                 res *= term^e
#             else
#                 res *= inv(term)^abs(e)
#             end
#         end
#     end
#     return res
# end

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

# ============================================================
# 6. Final Evaluator API
# ============================================================

@inline ishalfInt(j)::Bool = isinteger(2*j) && j ≥ 0
@inline δ(j1, j2, j3)::Bool = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)
@inline qδ(j1, j2, j3, k)::Bool = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 

@inline function qδtet(j1, j2, j3, j4, j5, j6, k)
    return ishalfInt(j1) && ishalfInt(j2) && ishalfInt(j3) && 
           ishalfInt(j4) && ishalfInt(j5) && ishalfInt(j6) && 
           qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && 
           qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)
end

"""
    exact_qracah6j(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6; prec=256)

Computes the 6j symbol exactly using number fields to prevent cancellation, 
then safely projects to BigFloat.
"""
function exact_qracah6j(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6; prec=256)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return zero(BigFloat)
    end
    
    setprecision(BigFloat, prec) do
        # 1. Exact Summation in the Field
        S_z = q6jseries_symb(j1, j2, j3, j4, j5, j6)
        sum_nf = model.K(0)
        for M in S_z
            sum_nf += evaluate_nf(M, model)
        end
        
        # 2. Exact Prefactor ^ 2 in the Field
        Pref2_symb = qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
        Pref2_nf = evaluate_nf(Pref2_symb, model)
        
        # 3. Project exactly to float using Horner's Method
        N = model.k + 2
        target_z = cispi(big"1.0" / N) # Target complex root for z
        
        sum_bf = evaluate_horner(sum_nf, target_z)
        pref2_bf = evaluate_horner(Pref2_nf, target_z)
        
        # 4. Final Math
        final_res = sqrt(pref2_bf) * sum_bf
        
        # Ensure it maps to purely real (stripping machine epsilon imaginary parts)
        return real(final_res)
    end
end

end # module