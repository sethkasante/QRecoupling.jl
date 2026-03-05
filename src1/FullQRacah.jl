module FullQRacah


using Nemo
using LRUCache

export CycloMonomial, qfactorial_symb, ExactSU2kModel, NumericSU2kModel,
        qdelta2_symb, qtricoeff2_symb, q6jseries_symb, q6jsummand_symb, evaluate_nf

struct CycloMonomial
    sign::Int
    q_pow::Int     # Handles the n(n-1)/2 fractional powers of [n]!.  
    exps::Dict{Int, Int} # Handles Π (Φ_d(q^2))^(exps[d]) # Note the polynomials are evaluated at powers are q^2 where q = e^{i pi/ (k+2)}
end

# Exact Symbolic Multiplication
function Base.:*(a::CycloMonomial, b::CycloMonomial)
    # loop over smaller dict
    small, large = length(a.exps) ≤ length(b.exps) ? (a,b) : (b,a)
    exps = copy(large.exps)
    for (d,e) in small.exps
        val = get(exps, d, 0) + e
        if val == 0
            delete!(exps, d)
        else
            exps[d] = val
        end
    end
    return CycloMonomial(a.sign * b.sign, a.q_pow + b.q_pow, exps)
end

# Exact Symbolic Division
function Base.:/(a::CycloMonomial, b::CycloMonomial)
    exps = copy(a.exps)
    for (d, e) in b.exps
        val = get(exps, d, 0) - e
        if val == 0
            delete!(exps, d)
        else
            exps[d] = val
        end
    end
    CycloMonomial(a.sign * b.sign, a.q_pow - b.q_pow, exps)
end



#---- Models  ---- 

struct ExactSU2kModel
    k::Int
    K::AnticNumberField            
    z::nf_elem                     # Primitive root
    Phi_cache::Dict{Int, nf_elem}  # Cached evaluations of Φ_d(z)
    # qfacts::Vector{CycloMonomial}  # Precomputed symbolic factorials
end

# Cache for precomputed log-q-factorials: key (k, DataType) -> Vector{BigFloat}
const PHI_CACHE = LRU{Int, nf_elem}(maxsize = 1024)
# const QFACT_SYMB_CACHE = LRU{Int, CycloMonomial}(maxsize = 1024)

function ExactSU2kModel(k::Int)
    N = k + 2
    # q = exp(i * pi / N) is a primitive 2N-th root of unity
    K, z = cyclotomic_field(2N, "z") 
    
    Zx, x = polynomial_ring(ZZ, "x") # x is z^2: to be used at evaluation
    
    # Cache evaluations of Φ_d(z)
    max_d = 2 * N 
    Phi_cache = Dict{Int, nf_elem}()
    for d in 2:max_d
        poly = cyclotomic(d, x)
        Phi_cache[d] = get!(PHI_CACHE , d) do
            evaluate(poly, z^2)  # use z^2 at evaluation
        end
    end
    
    # Precompute symbolic factorials
    # max_n = 2 * k + 4 
    # qfacts = Dict{Int, nf_CycloMonomialelem}() #Vector{CycloMonomial}(undef, max_n + 1)
    # for n in 0:max_n
    #     qfacts[n+1] = get!(QFACT_SYMB_CACHE, n) do 
    #         qfactorial_symb(n)
    #     end
    # end
    
    return ExactSU2kModel(k, K, z, Phi_cache) #, qfacts)
end

# function ExactSU2kModel(k::Int)
#     N = k + 2
#     # q = exp(i * pi / N) is a primitive 2N-th root of unity
#     K, z = cyclotomic_field(2N, "z")
    
#     Zx, x = polynomial_ring(ZZ, "x") # x is z^2: to be used at evaluation
    
#     # Cache evaluations of Φ_d(z)
#     max_d = 2 * N 
#     Phi_cache = Dict{Int, nf_elem}()
#     for d in 2:max_d
#         poly = cyclotomic(d, x)
#         Phi_cache[d] = get!(PHI_CACHE , d) do
#             evaluate(poly, z^2)  # use z^2 at evaluation
#         end
#     end
    
#     # Precompute symbolic factorials
#     max_n = 2 * k + 4 
#     qfacts = Vector{CycloMonomial}(undef, max_n + 1)
#     for n in 0:max_n
#         qfacts[n+1] = get!(QFACT_SYMB_CACHE, n) do 
#             qfactorial_symb(n)
#         end
#     end
    
#     return ExactSU2kModel(k, K, z, Phi_cache, qfacts)
# end





struct NumericSU2kModel{T<:AbstractFloat}
    k::Int
    logqnfact::Vector{T}    #Fallback for large k (Float64 or BigFloat)
end

# Cache for precomputed log-q-factorials: key (k, DataType) -> Vector{BigFloat}
const LOGQFACT_CACHE = LRU{Tuple{Int,Int}, Vector{BigFloat}}(maxsize = 1024)

function NumericSU2kModel(k::Int,prec=256)::NumericSU2kModel
    # Use current precision
    tab = get!(LOGQFACT_CACHE, (k, prec)) do
        logqnfact_table(k,prec) 
    end
    return NumericSU2kModel{BigFloat}(k, tab)
end

# using Base.MPFR

function logqn_table(k::Int, prec::Int=256)::Vector{BigFloat}
    N = k + 2
    setprecision(prec) do
        θ = big(pi) / BigFloat(N)
        logsinθ = log(sin(θ))

        tab = Vector{BigFloat}(undef, N)

        tab[1] = -Inf  # log([0]_q) not used; define safely
        # Use symmetry: n ↔ k+2-n
        half = N ÷ 2
        tab[1] = zero(BigFloat) # log([0]) = 0

        @inbounds for n in 1:half
            val = log(sin(n*θ)) - logsinθ
            tab[N + 1 - n] = tab[n+1] = val
        end
        return tab
    end
end


function logqnfact_table(k::Int, prec::Int=256)::Vector{BigFloat}
    logqn = logqn_table(k, prec)

    setprecision(prec) do
        tab = Vector{BigFloat}(undef, k+1)
        tab[1] = logqn[1]  # log([0]_q!) = 0
        @inbounds for n in 2:k+2
            tab[n] = tab[n-1] + logqn[n]
        end
        return tab
    end
end

#-------------


function qfactorial_symb(n::Int)
    if n == 0
        return CycloMonomial(1, 0, Dict{Int, Int}())
    end

    exps = Dict{Int, Int}()
    for d in 2:n
        power = div(n, d) # Integer division gives floor(n/d)
        if power > 0
            exps[d] = power
        end
    end
    
    # Calculate the q^1/2 power
    q_pow = - (n * (n - 1)) ÷ 2
    
    return CycloMonomial(1, q_pow, exps)
end

function qdelta2_symb(a, b, c)
    qfactorial_symb(Int(a+b-c)) * qfactorial_symb(Int(a-b+c)) * qfactorial_symb(Int(-a+b+c)) / qfactorial_symb(Int(a+b+c+1))
end


function qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    qdelta2_symb(j1, j2, j3) * qdelta2_symb(j1, j5, j6) * 
    qdelta2_symb(j2, j4, j6) * qdelta2_symb(j3, j4, j5)
end

function q6jsummand_symb(z,α1,α2,α3,α4,β1,β2,β3)
    num = qfactorial_symb(z+1)
    den = qfactorial_symb(z-α1) * qfactorial_symb(z-α2) * qfactorial_symb(z-α3) *
        qfactorial_symb(z-α4) * qfactorial_symb(β1-z) * qfactorial_symb(β2-z) *
        qfactorial_symb(β3-z) 
    res = num / den
    if isodd(z)
        return CycloMonomial(-1, 0, Dict{Int,Int}()) * res
    else
        return res
    end
end


function q6jseries_symb(j1, j2, j3, j4, j5, j6)::Vector{CycloMonomial}
    α1 = Int(j1 + j2 + j3) 
    α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6) 
    α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5)
    β2 = Int(j1 + j3 + j4 + j6) 
    β3 = Int(j2 + j3 + j5 + j6)
    #can this be faster?, the order doesn't matter
    S_z = CycloMonomial[]
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3)
    @inbounds for z in zrange
        push!(S_z, q6jsummand_symb(z, α1, α2, α3, α4, β1, β2, β3)) 
    end
    return S_z
end

#evaluate in cyclotomic fields 
function eevaluate_nf(M::CycloMonomial, model::ExactSU2kModel)
    K = model.K
    res = K(M.sign) * model.z^M.z_pow
    for (d, e) in M.exps
        res *= (model.Phi_cache[d])^e
    end
    return res
end

function eevaluate_nf(M::CycloMonomial, model::ExactSU2kModel)
    K = model.K
    # Start with the sign
    res = K(M.sign)
    
    # Handle q_pow (which is power of q = z)
    if M.q_pow != 0
        z_pow = M.q_pow > 0 ? model.z : inv(model.z)
        res *= z_pow^abs(M.q_pow)
    end
    
    # Handle Cyclotomic Polynomials
    for (d, e) in M.exps
        if e != 0
            term = model.Phi_cache[d]
            if e > 0
                res *= term^e
            else
                res *= inv(term)^abs(e)
            end
        end
    end
    return res
end

# evaluate in real field or complex field

end #end module 