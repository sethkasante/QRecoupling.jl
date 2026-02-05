module FactoredQuantum


export CyclotomicBasis,
       CycloMonomial, CycloExpr, CycloSum, 
       _qΔ2, _q6jsummand, 
       _qRacahsum, _qracah6j,
       qδ, δ, δtet, qδtet




# classical and quantum triangle conditions at level k: 
# checks admissible triple 
"""
    δ(j1, j2, j3,k) -> ::Bool
    qδ(j1, j2, j3,k) -> ::Bool

Checks the triangle conditions `j3 ≤ j1 + j2`, `j1 ≤ j2 + j3`, `j2 ≤ j3 + j1` and extra condition `j1 + j2 + j3 ≤ k` for qδ.
"""
δ(j1, j2, j3) = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)

# (j3<= j1+j2) && (j1 <= j2 + j3) && (j2 <= j3 + j1) 

qδ(j1, j2, j3, k) = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 


δtet(j1, j2, j3, j4, j5, j6) = δ(j1, j2, j3) && δ(j1, j5, j6) && δ(j2, j4, j6) && δ(j3, j4, j5)

qδtet(j1, j2, j3, j4, j5, j6) = qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)



"""
A monomial of the cyclotomic polynomial basis.
Represents ∏ Φ_d(q)^{e_d}. These are our basis objects
"""
struct CycloMonomial
    exponent::Dict{Int, Int}   # d => e : Φ_d^e 
end


struct CycloExpr
    unit_exp::Int     # exponent m such that unit = q^(m/2)
    cyclo::CycloMonomial
end

#set q = q0^2: so that units are in q0 and monomials are in q

struct CycloSum               # global common factors
    unit_exp::Vector{Int}
    coeffs::Vector{Int}
    terms::Vector{CycloMonomial}  # reduced monomials → coeff
end

CycloMonomial() = CycloMonomial(Dict{Int,Int}())


Base.length(m::CycloMonomial) = length(m.exponent)

function Base.:*(a::CycloMonomial, b::CycloMonomial)
    # loop over smaller dict
    small, large = length(a) ≤ length(b) ? (a,b) : (b,a)
    dict = copy(large.exponent)
    for (d,e) in small.exponent
        dict[d] = get(dict,d,0) + e
        dict[d] == 0 && delete!(dict,d)
    end
    return CycloMonomial(dict)
end

function Base.:*(a::CycloExpr, b::CycloExpr)
    CycloExpr(a.unit_exp + b.unit_exp,
        a.cyclo * b.cyclo)
end

"""
quantum triangle coefficient squared in cyclotomic form
"""
function _qΔ2(j1, j2, j3)
    δ(j1, j2, j3) || return CycloExpr(0, CycloMonomial())

    dict = Dict{Int,Int}()
    n1 = j1 + j2 - j3
    n2 = j1 - j2 + j3
    n3 = -j1 + j2 + j3
    n4 = j1 + j2 + j3 + 1
    sizehint!(dict,n4-1)

    for d in 2:n4
        exponent =
            floor(Int, n1/d) + floor(Int, n2/d) + 
            floor(Int, n3/d) - floor(Int, n4/d)
        exponent != 0 && (dict[d] = exponent)
    end

    # unit exponent from factorials

    unit_exp = -((n1*(n1-1) + n2*(n2-1) + n3*(n3-1) - n4*(n4-1)) ÷ 2 )

    return CycloExpr(unit_exp, CycloMonomial(dict))
end

"""
cyclo expression for each term in Racah sum
"""
function _q6jsummand(z,α1,α2,α3,α4,β1,β2,β3)
    dict = Dict{Int,Int}()

    #factorials 
    n1 = (z+1)
    a1, a2, a3, a4 = z-α1, z-α2, z-α3, z-α4
    b1, b2, b3 = β1-z, β2-z, β3-z
    N = max(z+1,β1-z,β2-z,β3-z)
    sizehint!(dict,N - 1)
    for d in 2:N
        exponent = floor(Int,n1/d) - 
                    floor(Int,a1/d) - floor(Int,a2/d) - floor(Int,a3/d) - 
                    floor(Int,a4/d) - floor(Int,b1/d) - floor(Int,b2/d) -
                    floor(Int,b3/d)
        exponent != 0 && (dict[d] = exponent)
    end
    # monomial unit exponents

    unit_expr = -((n1*(n1-1) - a1*(a1-1) - a2*(a2-1) - a3*(a3-1) -
                a4*(a4-1) - b1*(b1-1) - b2*(b2-1) - b3*(b3-1)) ÷ 2)

    return CycloExpr(unit_expr, CycloMonomial(dict))
end

"""
Sums of cyclotomic expressions
"""
function _qRacahsum(α1,α2,α3,α4,β1,β2,β3)
    #range of z values
    zrange = max(α1,α2,α3,α4):min(β1,β2,β3)
    N = length(zrange)
    # n = length(zrange)
    sgns = Int[] 
    sizehint!(sgns,N)
    terms = CycloMonomial[] # add sizehint = n 
    sizehint!(terms,N)
    unit_expr = Rational{Int}[] # add sizehint = n 
    sizehint!(unit_expr,N)
    for z in zrange
        cyclo_expr = _q6jsummand(z,α1,α2,α3,α4,β1,β2,β3)
        push!(terms, cyclo_expr.cyclo )
        push!(unit_expr, cyclo_expr.unit_exp )
        push!(sgns, iseven(z) ? 1 : -1 )
    end

    return CycloSum(unit_expr,sgns,terms)
end



function _qracah6j(j1, j2, j3, j4, j5, j6)
    
   
    #WignerSymbols reorders this for storage purposes 
    t1 = (j1, j2, j3)
    t2 = (j1, j5, j6)
    t3 = (j2, j4, j6)
    t4 = (j3, j4, j5)

    T2 = _qΔ2(t1...) * _qΔ2(t2...) * _qΔ2(t3...) * _qΔ2(t4...)
    
    #WignerSymbols reorders this for storage purposes 
    α1, α2, α3, α4 = j1 + j2 + j3, j1 + j5 + j6, j2 + j4 + j6, j3 + j4 + j5
    β1, β2, β3 = j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6

    S = _qRacahsum(α1,α2,α3,α4,β1,β2,β3)
    return T2, S
end


using Nemo

export val_unit, evaluate, val_unit, valuation, is_forced_zero, is_forced_poles, evaluate_6j

#note that the units are monomials in Z[q^1/2] and cyclotomic polynomials are Φ_d(q)
#set q0^2 = q and compute in terms of q0 

const _, x = polynomial_ring(ZZ, "x")

q0(k) = 

function evaluate(m::CycloMonomial,k)
    K, ζ = cyclotomic_field(2*(k+2), "ζ") #2(k+2) from q0
    val = one(ζ)
    for (d, e) in m.exponent
        val *= K(cyclotomic(d, x^2))^e   # fast numerical Φ_d(q)^e
    end
    return val
end


@inline function val_unit(exp::Int, k)
    #val = x^exp
    K, ζ = cyclotomic_field(2*(k+2), "ζ")
    #K(val)
    ζ^exp
end

function is_forced_zero(S::CycloExpr, k)
    valuation(S, k+2) > 0 
end

function is_forced_poles(S::CycloExpr, k)
    valuation(S, k+2) < 0 
end

# @inline 
function evaluate(e::CycloExpr, k)
    is_forced_zero(e,k+2) && return 0
    val = valuation(e, k+2)
    val < 0 && return throw(ArgumentError("The expression has poles of order $val"))
    val_unit(e.unit_exp, k) * evaluate(e.cyclo, k)
end

function evaluate(S::CycloSum,k)
    if is_forced_zero(S, k)
        return 0
    end

    total = 0
    @inbounds for i in eachindex(S.terms)
        term = evaluate(S.terms[i], k)
        total += S.coeffs[i] * val_unit(S.unit_exp[i], k) * term
    end
    return total
end


# ------ evaluation -------- 

"""
    eval(m::CycloMonomial, q)

Evaluate ∏_d Φ_d(q)^{e_d} numerically.
Assumes q is already specialized.
"""
# function eval2(m::CycloMonomial, q)
#     val = one(q)
#     for (d, e) in m.exponent
#         Φd = cyclotomic(d, q)   # fast numerical Φ_d(q)
#         val *= Φd^e
#     end
#     return val
# end

# @inline eval_unit(exp::Rational{Int}, q) = q^exp

"""
    eval(e::CycloExpr, q)

Evaluate q^(unit_exp) * cyclotomic monomial.
"""
# @inline function eval2(e::CycloExpr, q)
#     eval_unit(e.unit_exp, q) * eval(e.cyclo, q)
# end

"""
    eval(S::CycloSum, q; k=nothing)

Evaluate the Racah sum numerically.
If k is provided, valuation at d=k+2 is checked.
"""
function evalua(S::CycloSum, q; k=nothing)
    if k !== nothing && is_forced_zero(S, k)
        return zero(q)
    end

    total = zero(q)
    @inbounds for i in eachindex(S.terms)
        term = evaluate(S.terms[i], q)
        total += S.coeffs[i] * eval_unit(S.unit_exp[i], q) * term
    end
    return total
end

"""
    evaluate_6j(j1,...,j6; q, k=nothing)

Compute the quantum or classical 6j symbol.
"""
function evaluate_6j(j1,j2,j3,j4,j5,j6, k)
    T2, S = _qracah6j(j1,j2,j3,j4,j5,j6)

    # evaluate triangle product
    tri_val = evaluate(T2, k)

    # evaluate Racah sum
    sum_sq = evaluate(S, k)^2

    return (tri_val) * sum_sq
end

# function evaluate_6j_q(js...; k)
#     q = cis(2π/(k+2))
#     evaluate_6j(js...; q=q, k=k)
# end

"""
    is_forced_zero(S::CycloSum, k)

Returns true if valuation at d = k+2 is positive.
"""
function is_forced_zero(S::CycloSum, k)
    valuation(S, k+2) > 0
end

# ------ valuation ----- 

valuation(m::CycloMonomial, d::Int) = get(m.exponent, d, 0)

valuation(e::CycloExpr, d::Int) = valuation(e.cyclo, d)

function valuation(S::CycloSum, d::Int)
    vmin = typemax(Int)
    for m in S.terms
        v = valuation(m, d)
        v < vmin && (vmin = v)
    end
    return vmin
end

function support(S::CycloSum)
    ds = Set{Int}()
    for m in S.terms
        union!(ds, keys(m.exponent))
    end
    return ds
end

function valuation_q6j(j1, j2, j3, j4, j5, j6, d)
    T2, S = _qracah6j(j1, j2, j3, j4, j5, j6)
    return (valuation(T2, d) + valuation(S, d) )
end

# function valuation_q6j(j1, j2, j3, j4, j5, j6) 
#     T2, S = _qracah6j(j1, j2, j3, j4, j5, j6)  
#     vals = Tuple{Int64, Int64}[]
#     for d in support(S)
#         push!(vals, (d, valuation(T2, d) + valuation(S, d)) )
#     end
#     return vals
# end







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

# function evaluate(term::CycloMonomial, k::Int; prec=256)
#     val = Complex{BigFloat}(term.coeff)
#     for (d,e) in term.factors
#         φ = cyclotomic_eval(d, k; prec=prec)
#         val *= φ^e
#     end
#     return val
# end


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
