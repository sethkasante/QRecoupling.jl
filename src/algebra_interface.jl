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
"""
    LevelXValue

A value at a level written out in x = 2cos(π/(k+2)) as `P(x)√R(x)` over one integer denominator, from
[`x_form`](@ref) on an `Exact(k)` form. The polynomials are reduced below the degree of the number field, so they hold only
at this level. `QRecoupling.xpolynomial` and `QRecoupling.radicand` give their coefficients, lowest degree first.
"""
struct LevelXValue{E}
    value::E
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

const _NEMO = Base.PkgId(Base.UUID("2edaba10-b0f1-5616-af89-8c11ac63239a"), "Nemo")

"Why the backend is missing, which differs once Nemo itself has been loaded: then only the extension can be absent."
function _algebra_hint(name, nemo_loaded::Bool, top_level::Bool)
    nemo_loaded || return "$name requires the optional Nemo backend. Install Nemo and load it with `using Nemo`. Numerical values, `Exact()` and `Exact(k)` do not require it."
    top_level || return "$name requires the Nemo extension, and Nemo is loaded, but a QRecoupling made by `include` cannot have package extensions. Restart Julia and load it with `using QRecoupling, Nemo`."
    return "$name requires the Nemo extension, and Nemo is loaded, but the extension is not. Restart Julia and load both with `using QRecoupling, Nemo`. If it persists, check `Base.get_extension(QRecoupling, :QRecouplingNemoExt)` and that Nemo meets the version in QRecoupling's compat."
end
_missing_algebra(name) = throw(ArgumentError(_algebra_hint(name, Base.root_module_exists(_NEMO), parentmodule(@__MODULE__) === @__MODULE__)))
_qfact_exactx(args...; kwargs...) = _missing_algebra(:Exact)
exact_x(args...; kwargs...) = _missing_algebra(:exact_x)
generic_value(args...; kwargs...) = _missing_algebra(:x_form)
level_x_value(args...; kwargs...) = _missing_algebra(:x_form)
generic_sixj(args...; kwargs...) = _missing_algebra(:generic_sixj)
"""
    radical(v; kwargs...)

A value at a level in nested square roots, where they exist. It needs Nemo: load it with `using Nemo`.
"""
radical(args...; kwargs...) = _missing_algebra(:radical)
radical_form(args...; kwargs...) = _missing_algebra(:radical_form)
radical_levels(args...; kwargs...) = _missing_algebra(:radical_levels)
has_radical_form(args...; kwargs...) = _missing_algebra(:has_radical_form)
"""
    prove_identity(lhs, rhs; kwargs...)

A proof that two expressions in x = q + 1/q agree at every q. It needs Nemo: load it with `using Nemo`.
"""
prove_identity(args...; kwargs...) = _missing_algebra(:prove_identity)
xpolynomial(args...; kwargs...) = _missing_algebra(:xpolynomial)
radicand(args...; kwargs...) = _missing_algebra(:radicand)
"""
    phi_form(f; kwargs...)

A formula expanded in q and factored over the cyclotomic polynomials. It needs Nemo: load it with `using Nemo`.
"""
phi_form(args...; kwargs...) = _missing_algebra(:phi_form)
splits_completely(args...; kwargs...) = _missing_algebra(:splits_completely)
xvalue(rad, num) = _missing_algebra(:xvalue)
xvalue(rad, num, den) = _missing_algebra(:xvalue)
ExactXSum(k) = _missing_algebra(:ExactXSum)
function x_form end

"Largest generic polynomial the display of an `XValue` writes out, and its length."
const GENERIC_MAX_DEGREE = 16
const GENERIC_MAX_CHARS = 120
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
