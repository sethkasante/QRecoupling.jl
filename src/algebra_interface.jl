# Stable public wrappers. Their polynomial payload is supplied by the optional Nemo extension.
"An expanded generic expression; load Nemo to construct and manipulate it."
struct XValue{P}
    rad::Vector{Int}
    num::P
    den::P
end
"An exact level value P(x)*sqrt(R(x)); the optional backend supplies its polynomials."
struct ExactX{P}
    k::Int
    p::P
    sqclass::Vector{Int}
    r::P
    rule::Union{Nothing,FactorialSum}
end
"A sum of exact level values, grouped by radical class."
struct ExactXSum{P}
    k::Int
    terms::Dict{Vector{Int},P}
end
"A nested real radical expression returned by the optional algebra backend."
struct RadExpr
    rat::Rational{BigInt}
    terms::Vector{Tuple{Rational{BigInt},RadExpr}}
end
"The reason a compact square-root expression was not returned."
struct NoRadical
    k::Int
    kind::Symbol
    degree::Int
    field::Int
    limit::Int
    approx::Union{Nothing,Float64}
end

_missing_algebra(name) = throw(ArgumentError("$name requires the optional Nemo backend. Install Nemo and load it with `using Nemo`. Numerical values, `Exact()` and `Exact(k)` do not require it."))
_qfact_exactx(args...; kwargs...) = _missing_algebra(:Exact)
exact_x(args...; kwargs...) = _missing_algebra(:exact_x)
generic_value(args...; kwargs...) = _missing_algebra(:generic_value)
generic_sixj(args...; kwargs...) = _missing_algebra(:generic_sixj)
radical(args...; kwargs...) = _missing_algebra(:radical)
radical_form(args...; kwargs...) = _missing_algebra(:radical_form)
radical_levels(args...; kwargs...) = _missing_algebra(:radical_levels)
has_radical_form(args...; kwargs...) = _missing_algebra(:has_radical_form)
prove_identity(args...; kwargs...) = _missing_algebra(:prove_identity)
xpolynomial(args...; kwargs...) = _missing_algebra(:xpolynomial)
radicand(args...; kwargs...) = _missing_algebra(:radicand)
phi_form(args...; kwargs...) = _missing_algebra(:phi_form)
splits_completely(args...; kwargs...) = _missing_algebra(:splits_completely)
xvalue(rad, num) = _missing_algebra(:xvalue)
xvalue(rad, num, den) = _missing_algebra(:xvalue)
ExactXSum(k) = _missing_algebra(:ExactXSum)
function x_form end
function numeric_value end
Base.zero(::Type{ExactX}, k) = _missing_algebra(:Exact)

"""
    numeric_value(v::ClassicalResult; bits=precision(BigFloat)) -> BigFloat

The classical exact value σ√r, rounded to `bits` bits.
"""
function numeric_value(v::ClassicalResult; bits::Integer=precision(BigFloat))
    bits > 0 || throw(ArgumentError("bits must be positive"))
    target=Int(bits)
    result=setprecision(() -> v.sign * sqrt(BigFloat(v.sq_val)), BigFloat, target+32)
    return setprecision(()->BigFloat(result),BigFloat,target)
end


"""
    psi_exponents(pairs) -> Dict{Int,Int}

Multiplicity of each `ψ_e` in `∏ₙ [n]!^{cₙ}`, from `pairs` a collection of `n => c`. Since
`[n]! = ∏_{m≤n} [m]` and `[m] = ∏_{e | 2m, e ≥ 3} ψ_e`, the multiplicity of `ψ_e` is a divisor count:
no factoring is involved.
"""
function psi_exponents(pairs)
    isempty(pairs) && return Dict{Int,Int}()
    N = maximum(Int(p.first) for p in pairs)
    N <= 0 && return Dict{Int,Int}()
    E = Dict{Int,Int}()
    # ψ_e divides [m] exactly when e divides 2m. Count those multiples in n! directly.
    for e in 3:2N
        step = iseven(e) ? e ÷ 2 : e
        v = sum(Int(c) * (Int(n) ÷ step) for (n,c) in pairs; init=0)
        iszero(v) || (E[e] = v)
    end
    return E
end

"""
Halve a ψ exponent vector into square-part exponents and odd-multiplicity indices, still symbolic, so the
square part can merge with the first term before anything is expanded.
"""
function _halve_exponents(E::Dict{Int,Int})
    half = Dict{Int,Int}(); rad = Int[]
    for e in sort!(collect(keys(E)))
        v = E[e]
        r = mod(v, 2)
        r == 1 && push!(rad, e)
        f = (v - r) ÷ 2
        iszero(f) || (half[e] = f)
    end
    return half, rad
end
