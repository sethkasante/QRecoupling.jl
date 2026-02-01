module FactoredQuantum

using Primes
using LRUCache

export CyclotomicBasis,
       CycloMonomial,CycloSum, _qRacahsum, 
       qinteger,
       qfactorial,
       qΔ2, qδ, δ, 
       DIVISOR_CACHE, QINTEGER_CACHE, QFACTORIAL_CACHE,
       q6jsummand,
       qRacahsum, factored_sum,
       _qracah6j,
       cyclotomic_eval, q_root

# using LinearAlgebra

# """
# CyclotomicBasis(N)

"""
A monomial of the cyclotomic polynomial basis.
Represents ∏ Φ_d(q)^{e_d}. These are our basis objects
"""
struct CycloMonomial
    powers::Dict{Int, Int}   # d => powers of Φ_d 
end

struct CycloExpr
    unit_exp::Rational{Int}              # exponent m such that unit = q^(m/2)
    cyclo::CycloMonomial
end

# struct CycloMonomial
#     exponent::Dict{Int, Int}   # d => powers of Φ_d 
# end


# struct CycloMonomial
#     powers::Dict{Int,Int}
#     function CycloTerm(powers::Dict{Int,Int})
#         # remove zero exponents automatically
#         new(Dict(filter((d,e)->e!=0, powers)))
#     end
# end

#TODO: can we set it up such that any element with value 0 is automatically removed?

CycloMonomial() = CycloMonomial(Dict{Int,Int}())

CycloMonomial(d::Int, e::Int=1) = CycloMonomial(Dict(d => e))


const DIVISOR_CACHE = LRU{Int, Vector{Int}}(maxsize=10_000)
const QINTEGER_CACHE = LRU{Int, CycloMonomial}(maxsize=10_000)
const QFACTORIAL_CACHE = LRU{Int, CycloMonomial}(maxsize=5_000)

function Base.:*(a::CycloMonomial, b::CycloMonomial)
    expo = copy(a.powers)
    for (d,e) in b.powers
        expo[d] = get(expo,d,0) + e
        # exp[d] == 0 && delete!(exp,d)
    end
    return CycloMonomial(expo)
end

Base.:/(a::CycloMonomial, b::CycloMonomial) =
    a * CycloMonomial(Dict(d => -e for (d,e) in b.powers))

Base.://(a::CycloMonomial, b::CycloMonomial) = a / b

function Base.:^(a::CycloMonomial, n::Int)
    CycloMonomial(Dict(d => n*e for (d,e) in a.powers))
end



@inline function divisors_cached(n::Int)
    n == 0 && return Int[]
    return get!(DIVISOR_CACHE, n) do
        return divisors(n) # use divisors from Prime 
    end
end


"""
qinteger(n)

Returns cyclotomic representation of integer [n]_q
"""
function qinteger(n::Int)
    get!(QINTEGER_CACHE, n) do
        divs = divisors_cached(n)[2:end] # omit d = 1 
        CycloMonomial(Dict(divs .=> 1))
    end
end

"""
qfactorial(n)

Returns cyclotomic representation of [n]_q!
"""
function qfactorial(n::Int)
    n <= 1 && return CycloMonomial()
    get!(QFACTORIAL_CACHE, n) do
        qfactorial(n-1) * qinteger(n)
    end
end

# function _qfactorial(n::Int)
#     n <= 1 && return CycloMonomial()
#     dict = Dict{Int,Int}()
#     for d in 2:n
#         dict[d] = floor(n/d)
#     end
#     return CycloMonomial(dict)
# end


# classical and quantum triangle conditions at level k: 
# checks admissible triple 
"""
    δ(j1, j2, j3,k) -> ::Bool
    qδ(j1, j2, j3,k) -> ::Bool

Checks the triangle conditions `j3 ≤ j1 + j2`, `j1 ≤ j2 + j3`, `j2 ≤ j3 + j1` and extra condition `j1 + j2 + j3 ≤ k` for qδ.
"""
δ(j1, j2, j3) = isinteger(j1+j2+j3) && (j3 <= j1 + j2) && (j1 <= j2 + j3) && (j2 <= j3 + j1) 

qδ(j1, j2, j3, k) = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 


δtet(j1, j2, j3, j4, j5, j6) = δ(j1, j2, j3) && δ(j1, j5, j6) && δ(j2, j4, j6) && δ(j3, j4, j5)

qδtet(j1, j2, j3, j4, j5, j6) = qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && qδ(j2, j4, j6, k), qδ(j3, j4, j5, k)


"""
quantum triangle coefficient in cyclotomic representation
"""
function qΔ2(j1,j2,j3)
    qfactorial(j1 + j2 - j3) * qfactorial(j1 - j2 + j3) * 
        qfactorial(-j1 + j2 + j3) / qfactorial(j1 + j2 + j3 + 1)
end



# function _qΔ2(j1,j2,j3)
#     δ(j1, j2, j3) || return CycloMonomial()
#     dict = Dict{Int,Int}()
#     for d in 2:j1+j2+j3+1
#         exponent = floor(Int,(j1+j2-j3)/d) + floor(Int,(j1-j2+j3)/d) + 
#                     floor(Int,(-j1+j2+j3)/d) - floor(Int,(j1+j2+j3+1)/d)
#         exponent != 0 ? dict[d] = exponent : nothing
#     end
#     return CycloMonomial(dict)
# end

function _qΔ2(j1, j2, j3)
    δ(j1, j2, j3) || return CycloExpr(0, CycloMonomial())

    dict = Dict{Int,Int}()
    n1 = j1 + j2 - j3
    n2 = j1 - j2 + j3
    n3 = -j1 + j2 + j3
    n4 = j1 + j2 + j3 + 1

    for d in 2:n4
        exponent =
            floor(n1/d) +
            floor(n2/d) +
            floor(n3/d) -
            floor(n4/d)

        exponent != 0 && (dict[d] = exponent)
    end

    # unit exponent from factorials

    unit_exp = -((n1*(n1-1) + n2*(n2-1) + n3*(n3-1) - n4*(n4-1)) // 4)

    return CycloExpr(unit_exp, CycloMonomial(dict))
end

function _q6jsummand(z,α1,α2,α3,α4,β1,β2,β3)
    dict = Dict{Int,Int}()
    for d in 2:max(z+1,β1-z,β2-z,β3-z)
        exponent = floor(Int,(z+1)/d) - floor(Int,(z-α1)/d) - floor(Int,(z-α2)/d) - 
                    floor(Int,(z-α3)/d) - floor(Int,(z-α4)/d) - floor(Int,(β1-z)/d) -
                        floor(Int,(β2-z)/d) - floor(Int,(β3-z)/d)
        exponent != 0 ? dict[d] = exponent : nothing
    end
    return CycloMonomial(dict)
end

function _qRacahsum(α1,α2,α3,α4,β1,β2,β3)
    #range of z values
    zrange = max(α1,α2,α3,α4):min(β1,β2,β3)
    # n = length(zrange)
    sgns = Int[] 
    terms = CycloMonomial[] # add sizehint = n 
    for z in zrange
        push!(terms, _q6jsummand(z,α1,α2,α3,α4,β1,β2,β3) )
        push!(sgns, iseven(z) ? 1 : -1 )
    end
    return CycloSum(terms,sgns)
end
# ----- dealing with summation -----

"""
Restructure sums into  Σ (-1)^z_i F_i =  C ⋅ Σ (-1)^z_i (R_i):
  C is common factor, R_i is residuals F_i/C and coeffs are signs (-1)^z_i
"""

struct CycloSum               # global common factors
    terms::Vector{CycloMonomial}  # reduced monomials → coeff
    coeffs::Vector{Int}
end


"""
individual elements in the Racah sum 
"""
function q6jsummand(z,α1,α2,α3,α4,β1,β2,β3)
    num = qfactorial(z+1)
    den = qfactorial(z-α1) * qfactorial(z-α2) * qfactorial(z-α3) * 
            qfactorial(z-α4) * qfactorial(β1-z) * qfactorial(β2-z) *
                qfactorial(β3-z)
    return num / den
end

# look for common factors in a list of CycloMonomials  -> CycloMonomial
function common_factors(terms)
    common = Dict{Int,Int}()
    all_d = Set{Int}()
    for m in terms
        for d in keys(m.powers)
            push!(all_d, d)
        end
    end
    for d in all_d
        cmin = minimum(get(m.powers,d,0) for m in terms)
        cmin != 0 && (common[d] = cmin)
    end
    return CycloMonomial(common)
end

function residual_terms(terms, common) #terms::Vector{CycloMonomial}
    return [m / common for m in terms]
end

function factored_sum(terms)
    cm = common_factors(terms) 
    residue = isempty(cm.powers) ? terms : residual_terms(terms,cm)
    return cm, residue 
end


"""
Racah sum as a sum over cyclotomic polynomial basis -> CycloSum
"""
# α1,α2,α3,α4,β1,β2,β3 are functions of spins js
function qRacahsum(α1,α2,α3,α4,β1,β2,β3)
    #range of z values
    zrange = max(α1,α2,α3,α4):min(β1,β2,β3)
    # n = length(zrange)
    sgns = Int[] 
    terms = CycloMonomial[] # add sizehint = n 
    for z in zrange
        push!(terms, q6jsummand(z,α1,α2,α3,α4,β1,β2,β3) )
        push!(sgns, iseven(z) ? 1 : -1 )
    end
    cm, res = factored_sum(terms)
    return cm, res, sgns
    # CycloSum(cm,res,sgns)
end


function _qracah6j(j1, j2, j3, j4, j5, j6)
    
    # check quantum triangle conditions for all faces in tetrahedron
    !(δtet(j1, j2, j3, j4, j5, j6)) || throw(ArgumentError("The input spins are inadmissible"))
    # if !(δtet(j1, j2, j3, j4, j5, j6))
    #     return CycloMonomial()
    # end
    
    #WignerSymbols reorders this for storage purposes 
    α1, α2, α3, α4 = j1 + j2 + j3, j1 + j5 + j6, j2 + j4 + j6, j3 + j4 + j5
    β1, β2, β3 = j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6

    #WignerSymbols reorders this for storage purposes 
    t1 = (j1, j2, j3)
    t2 = (j1, j5, j6)
    t3 = (j2, j4, j6)
    t4 = (j3, j4, j5)

    tri2_coeff = qΔ2(t1...) * qΔ2(t2...) *
                    qΔ2(t3...) * qΔ2(t4...)

    cm, res, sgns = qRacahsum(α1,α2,α3,α4,β1,β2,β3)
    #square root term squared 
    R = simplify(tri2_coeff * cm^2)
    res = [simplify(r) for r in res]
    S = CycloSum(res,sgns)
    return R, S
end

function simplify(ct::CycloMonomial)
    CycloMonomial(Dict(filter((d,e)->e!=0, ct.powers)))
end

function simplify(ct::CycloMonomial)
    for (d,e) in collect(ct.powers)
        if e == 0
            delete!(ct.powers, d)
        end
    end
    return ct
end

#-----------Now evaluation k,needed ------

#note q-admissible relation ensure no d > k+2 in the evaluation.




# function q6jsummand(z,α1,α2,α3,α4,β1,β2,β3)
#     num = qfactorial_factored(z+1)
#     den = qfactorial_factored(z-α1) * qfactorial_factored(z-α2) *
#             qfactorial_factored(z-α3) * qfactorial_factored(z-α4) *
#             qfactorial_factored(z-β1) * qfactorial_factored(z-β2) *
#             qfactorial_factored(z-β3)
#     return num / den
# end




# function q6jsummand(z,α1,α2,α3,α4,β1,β2,β3)
#     num = qfactorial_factored(z+1)
#     den = qfactorial_factored(z-α1) * qfactorial_factored(z-α2) *
#             qfactorial_factored(z-α3) * qfactorial_factored(z-α4) *
#             qfactorial_factored(z-β1) * qfactorial_factored(z-β2) *
#             qfactorial_factored(z-β3)
#     return num / den 
#     # sgn = iseven(z) ? 1 : -1
#     # return sgn, factors # sign and CycloMonomial
#     # return FactoredMonomial( iseven(z) ? 1 : -1 , factors) 
# end














function _divisors(n::Int)
    n ≤ 0 && return Int[]
    ds = Int[]
    for d in 1:isqrt(n)
        if n % d == 0
            push!(ds, d)
            d*d != n && push!(ds, div(n,d))
        end
    end
    return ds
end

function mobius(n::Int)
    n == 1 && return 1
    f = Dict{Int,Int}()
    m = n
    p = 2
    while p*p ≤ m
        while m % p == 0
            f[p] = get(f,p,0) + 1
            m ÷= p
        end
        p += 1
    end
    m > 1 && (f[m] = 1)

    for (_,e) in f
        e > 1 && return 0
    end
    return (-1)^length(f)
end

function q_root(k::Int; prec=256)
    setprecision(prec) do
        return cis(big(2pi)/(k+2))
    end
end

function cyclotomic_eval(n::Int, k::Int; prec=256)
    q = q_root(k; prec=prec)
    val = one(Complex{BigFloat})

    for d in divisors(n)
        μ = mobius(div(n,d))
        μ == 0 && continue
        val *= (q^d - 1)^μ
    end

    return val
end

function evaluate(term::CycloMonomial, k::Int; prec=256)
    val = Complex{BigFloat}(term.coeff)
    for (d,e) in term.factors
        φ = cyclotomic_eval(d, k; prec=prec)
        val *= φ^e
    end
    return val
end


# function evaluate(fe::FactoredElement, Q::QSU2kField)
#     q = qroot(Q)
#     result = one(Q.K)
#     for (i, d) in enumerate(fe.basis.divisors)
#         e = fe.powers[i]
#         e == 0 && continue
#         Φd = cyclotomic_polynomial(Q.K, d, q)
#         result *= Φd^e
#     end
#     return result
# end


end #end module


# export CycloMonomial,
#        FactoredElement,
#        qinteger_factored,
#        qfactorial_factored

# """
# Factored element in the cyclotomic polynomial basis.
# Represents ∏ Φ_d(q)^{e_d}.
# """
# struct CycloMonomial
#     exp::Dict{Int, Int}   # d ↦ exponent of Φ_d
# end

# CycloMonomial() = CycloMonomial(Dict{Int,Int}())

# CycloMonomial(d::Int, e::Int=1) =
#     CycloMonomial(Dict(d => e))


# function Base.:*(a::CycloMonomial, b::CycloMonomial)
#     exp = copy(a.powers)
#     for (d,e) in b.powers
#         exp[d] = get(exp,d,0) + e
#         exp[d] == 0 && delete!(exp,d)
#     end
#     return CycloMonomial(exp)
# end

# Base.:/(a::CycloMonomial, b::CycloMonomial) =
#     a * CycloMonomial(Dict(d => -e for (d,e) in b.powers))

# function Base.:^(a::CycloMonomial, n::Int)
#     CycloMonomial(Dict(d => n*e for (d,e) in a.powers))
# end

# function qinteger_factored(n::Int)
#     exp = Dict{Int,Int}()
#     for d in divisors(2n)
#         if d > 2
#             exp[d] = get(exp,d,0) + 1
#         end
#     end
#     CycloMonomial(exp)
# end

# function qfactorial_factored(n::Int)
#     f = CycloMonomial()
#     for m in 1:n
#         f *= qinteger_factored(m)
#     end
#     return f
# end

# end
