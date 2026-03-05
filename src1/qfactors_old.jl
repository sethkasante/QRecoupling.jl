module FactoredQuantum

using Primes
using LRUCache

export CyclotomicBasis,
       CycloMonomial,FactoredSum, 
       qinteger_factored,
       qfactorial_factored,
       delta2_factored,
       DIVISOR_CACHE, QINT_FACTOR_CACHE, QFACT_FACTOR_CACHE,
       cyclotomic_eval,
       q6jsummand,
       q_root,
       qRacahsum_factored, factored_sum,
       _qracah6j

# using LinearAlgebra

# """
# CyclotomicBasis(N)

# Stores divisors of N and an index map.
# Used to represent products of Φ_d(q).
# """
# struct CyclotomicBasis
#     k::Int
#     N::Int 
#     # divisors::Vector{Int}
#     # index::Dict{Int,Int}
# end


"""
Factored element in the cyclotomic polynomial basis.
Represents ∏ Φ_d(q)^{e_d}.
"""
struct CycloMonomial
    exp::Dict{Int, Int}   # d ↦ exponent of Φ_d 
end

# struct CycloMonomial2
#     d::Vector{Int}
#     e::Vector{Int}
# end

CycloMonomial() = CycloMonomial(Dict{Int,Int}())

CycloMonomial(d::Int, e::Int=1) =
    CycloMonomial(Dict(d => e))


const DIVISOR_CACHE = LRU{Int, Vector{Int}}(maxsize=10_000)
const QINT_FACTOR_CACHE = LRU{Int, CycloMonomial}(maxsize=10_000)
const QFACT_FACTOR_CACHE = LRU{Int, CycloMonomial}(maxsize=5_000)

function Base.:*(a::CycloMonomial, b::CycloMonomial)
    exp = copy(a.exp)
    for (d,e) in b.exp
        exp[d] = get(exp,d,0) + e
        # exp[d] == 0 && delete!(exp,d)
    end
    return CycloMonomial(exp)
end

Base.:/(a::CycloMonomial, b::CycloMonomial) =
    a * CycloMonomial(Dict(d => -e for (d,e) in b.exp))

Base.://(a::CycloMonomial, b::CycloMonomial) = a / b

function Base.:^(a::CycloMonomial, n::Int)
    CycloMonomial(Dict(d => n*e for (d,e) in a.exp))
end

# function Base.show(io::IO, c::CycloMonomial)
#     isempty(c.exp) && return print(io,"1")
#     terms = sort(collect(c.exp))
#     strs = String[]
#     for (d,e) in terms
#         push!(strs, e == 1 ? "Φ_$d" : "Φ_$d^$e")
#     end
#     print(io, join(strs," "))
# end

@inline function divisors_cached(n::Int)
    n == 0 && return Int[]
    n = abs(n)

    return get!(DIVISOR_CACHE, n) do
        divs = Int[]
        r = isqrt(n)
        @inbounds for i in 1:r
            if n % i == 0
                push!(divs, i)
                j = div(n, i)
                if j != i
                    push!(divs, j)
                end
            end
        end
        sort!(divs)   # optional but useful for deterministic ordering
        return divs
    end
end


# function _divisors(n::Int)
#     if n == 0
#         return Int[]
#     end
#     n = abs(n) # work with the absolute value for positive divisors
#     divs = Int[]
#     for i = 1:isqrt(n)
#         if n % i == 0
#             push!(divs, i)
#             if i * i != n
#                 push!(divs, div(n, i))
#             end
#         end
#     end
#     return divs
# end



# function CyclotomicBasis(N::Int)
#     divs = divisors(N)
#     # idx  = Dict(d => i for (i,d) in enumerate(divs))
#     return CyclotomicBasis(N, divs) #, idx)
# end


"""
qinteger_factored(n, basis)

Returns factored representation of [n]_q
"""
# function qinteger_factored(n::Int)
#     exp = Dict{Int,Int}()
#     for d in _divisors(2n)
#         if d > 2
#             exp[d] = get(exp,d,0) + 1
#         end
#     end
#     CycloMonomial(exp)
# end

function qinteger_factored(n::Int)
    get!(QINT_FACTOR_CACHE, n) do
        exp = Dict{Int,Int}()
        for d in divisors_cached(n)
            d > 1 && (exp[d] = 1)
        end
        CycloMonomial(exp)
    end
end

function qfactorial_factored(n::Int)
    get!(QFACT_FACTOR_CACHE, n) do
        f = CycloMonomial()
        for m in 1:n
            f *= qinteger_factored(m)
        end
        f
    end
end

function qfactorial_factored3(n::Int)
    n <= 1 && return CycloMonomial()
    get!(QFACT_FACTOR_CACHE, n) do
        qfactorial_factored(n-1) * qinteger_factored(n)
    end
end

function qfactorial_factored_floor(n::Int)
    get!(QFACT_FACTOR_CACHE, n) do
        factors = Dict{Int,Int}()
        for d in 1:n
            d > 1 && (factors[d] = floor(n/d))
        end
        f
    end
end
"""
qfactorial_factored(n, basis)

Returns factored representation of [n]_q!
"""
# function qfactorial_factored(n::Int)
#     f = CycloMonomial()
#     for m in 1:n
#         f *= qinteger_factored(m)
#     end
#     return f
# end

# q triangle condition at level k: checks admissible triple 
"""
    qδ(j1, j2, j3,k) -> ::Bool

Checks the triangle conditions `j3 ≤ j1 + j2`, `j1 ≤ j2 + j3`, `j2 ≤ j3 + j1` and j1+j2+j3 ≤ k.
"""
qδ(j1, j2, j3, k) = isinteger(j1+j2+j3) && (j3 <= j1 + j2) && (j1 <= j2 + j3) && (j2 <= j3 + j1) && (j1+j2+j3) <= k 

dt(j1, j2, j3) = isinteger(j1+j2+j3) && (j3 <= j1 + j2) && (j1 <= j2 + j3) && (j2 <= j3 + j1) 

function dtall(j1, j2, j3, j4, j5, j6)
    t1 = (j1, j2, j3)
    t2 = (j1, j5, j6)
    t3 = (j2, j4, j6)
    t4 = (j3, j4, j5)
    dt(t1...), dt(t2...), dt(t3...), dt(t4...)
end

function _check(j1, j2, j3, j4, j5, j6)
    return j1 + j2 + j3, j1 + j5 + j6, j2 + j4 + j6, j3 + j4 + j5,  j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6
end


function delta2_factored(a,b,c)
    qfactorial_factored(a+b-c) *
    qfactorial_factored(a-b+c) *
    qfactorial_factored(-a+b+c) /
    qfactorial_factored(a+b+c+1)
end


function _qracah6j(j1, j2, j3, j4, j5, j6, k)
    
    # check quantum triangle conditions
    if !(qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k))
        return CycloMonomial()
    end
    
    #
    α1, α2, α3, α4 = j1 + j2 + j3, j1 + j5 + j6, j2 + j4 + j6, j3 + j4 + j5
    β1, β2, β3 = j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6

    #WignerSymbols reorders this for storage purposes 
    t1 = (j1, j2, j3)
    t2 = (j1, j5, j6)
    t3 = (j2, j4, j6)
    t4 = (j3, j4, j5)

    tri2_coeff = delta2_factored(t1...) * delta2_factored(t2...) *
                    delta2_factored(t3...) * delta2_factored(t4...)

    sumseries = qRacahsum_factored(α1,α2,α3,α4,β1,β2,β3)

    tri2_coeff, sumseries
end





function q6jsummand(z,α1,α2,α3,α4,β1,β2,β3)
    num = qfactorial_factored(z+1)
    den = qfactorial_factored(z-α1) * qfactorial_factored(z-α2) *
            qfactorial_factored(z-α3) * qfactorial_factored(z-α4) *
            qfactorial_factored(z-β1) * qfactorial_factored(z-β2) *
            qfactorial_factored(z-β3)
    return num / den
end



"""
Restructure sums into  Σ (-1)^z_i F_i =  C ⋅ Σ (-1)^z_i (R_i):
  C is common factor, R_i is residuals F_i/C and coeffs are signs (-1)^z_i
"""

struct FactoredSum
    common::CycloMonomial                # global common factors
    residuals::Vector{CycloMonomial}  # reduced monomials → coeff
    coeffs::Vector{Int}
end

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


#TODO: write a function to filter and remove terms that are 0, but this need k  
function common_factors(terms)
    common = Dict{Int,Int}()
    alld = Set{Int}()
    for m in terms
        for d in keys(m.exp)
            push!(alld, d)
        end
    end
    for d in alld
        cmin = minimum(get(m.exp,d,0) for m in terms)
        cmin > 0 && (common[d] = cmin)
    end
    return CycloMonomial(common)
end


function residual_terms(terms, common)
    return [m / common for m in terms]
end

function factored_sum(terms)
    cm = common_factors(terms) 
    residue = isempty(cm.exp) ? terms : residual_terms(terms,cm)
    return cm, residue 
end


# α1,α2,α3,α4,β1,β2,β3 are functions of spins js
function qRacahsum_factored(α1,α2,α3,α4,β1,β2,β3)
    zrange = max(α1,α2,α3,α4):min(β1,β2,β3)
    # n = length(zrange)
    sgns = Int[] #ones(Int,n)  #Vector{Int}(undef,n)
    terms = CycloMonomial[] # add sizehint = n 
    for z in zrange
        push!(terms, q6jsummand(z,α1,α2,α3,α4,β1,β2,β3) )
        push!(sgns, iseven(z) ? 1 : -1 )
    end
    cm, res = factored_sum(terms)
    FactoredSum(cm,res,sgns)
end









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
#         e = fe.exp[i]
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
#     exp = copy(a.exp)
#     for (d,e) in b.exp
#         exp[d] = get(exp,d,0) + e
#         exp[d] == 0 && delete!(exp,d)
#     end
#     return CycloMonomial(exp)
# end

# Base.:/(a::CycloMonomial, b::CycloMonomial) =
#     a * CycloMonomial(Dict(d => -e for (d,e) in b.exp))

# function Base.:^(a::CycloMonomial, n::Int)
#     CycloMonomial(Dict(d => n*e for (d,e) in a.exp))
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
