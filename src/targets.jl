# --------------------------------------------
#  Evaluation targets
#
#  `q6j(Level(10), j...)` instead of `q6j(j...; k = 10)`. The keyword form keeps working; targets exist so
#  that a valid request is a single value rather than a combination of keywords, and so that a target can
#  be stored, passed around and swept over.
# --------------------------------------------

"""
    EvalTarget

Where a symbol should be evaluated: [`Level`](@ref), [`Exact`](@ref), [`At`](@ref), [`Classical`](@ref), or [`Symbolic`](@ref).
Any of them can be passed as the first argument of a symbol function in place of the `k`, `q`, `exact` and
`T` keywords.
"""
abstract type EvalTarget end

"""
    Level(k; T = Float64)

Floating-point values at the root of unity q = exp(iπ/(k+2)). `k` may be a range, which sweeps levels.

```julia
q6j(Level(10), 1, 1, 1, 1, 1, 1)
q6j(Level(2:50), 1, 1, 1, 1, 1, 1)          # one value per level
q6j(Level(10; T = BigFloat), 1, 1, 1, 1, 1, 1)
```
"""
struct Level{T,K} <: EvalTarget
    k::K
end
Level(k; T::Type = Float64) = Level{T,typeof(k)}(k)

"""
    Exact(k)

Exact values in the cyclotomic field ℚ(ζ_{2(k+2)}). Symbols that vanish are recognised by the structural
and modular zero tests, so an exact zero costs microseconds rather than milliseconds.
"""
struct Exact{K} <: EvalTarget
    k::K
end

"Internal target selected by Exact(k; backend=:native)."
struct NativeExact <: EvalTarget
    k::Int
end

"""
    Exact(k::Integer; backend=:nemo)

Choose the existing Nemo projection or the experimental `:native` direct
factorial-rule backend. Native results use BigInt polynomial arithmetic and a
positive radical; see [`exact_level`](@ref) for supported rules and limitations.
The package still loads Nemo for its existing symbolic/exact APIs.
"""
function Exact(k::Integer;backend::Symbol=:nemo)
    k>=0 || throw(DomainError(k,"level must be nonnegative"))
    backend===:nemo && return Exact{typeof(k)}(k)
    backend===:native && return NativeExact(Int(k))
    throw(ArgumentError("exact backend must be :nemo or :native"))
end
Base.show(io::IO,t::NativeExact)=print(io,"Exact(",t.k,"; backend=:native)")
target_kwargs(::NativeExact)=throw(ArgumentError("native exact target requires a supported scalar symbol or FactorialSum"))

qeval(t::NativeExact,s::FactorialSum;workspace=nothing)=exact_level(s,t.k;workspace=workspace)
function q6j(t::NativeExact,js::Vararg{Spin,6};workspace=nothing)
    J=doubled(js...)
    s=_qδtet(J...,t.k) ? sixj_sum(J...) : EMPTY_FACTORIAL_SUM
    exact_level(s,t.k;workspace=workspace)
end
function q3j(t::NativeExact,j1::Spin,j2::Spin,j3::Spin,m1::Spin,m2::Spin,m3::Spin=-m1-m2;workspace=nothing)
    J=doubled(j1,j2,j3,m1,m2,m3)
    s=_qδ(J[1],J[2],J[3],t.k) ? threej_sum(J...) : EMPTY_FACTORIAL_SUM
    exact_level(s,t.k;workspace=workspace)
end
for (f,rule) in ((:fsymbol,:fsymbol_sum),(:gsymbol,:gsymbol_sum),(:tetrahedron,:tetrahedron_sum))
    @eval function $f(t::NativeExact,js::Vararg{Spin,6};workspace=nothing)
        J=doubled(js...)
        s=_qδtet(J...,t.k) ? $rule(J...) : EMPTY_FACTORIAL_SUM
        exact_level(s,t.k;workspace=workspace)
    end
end
function evaluate_exact(r::NativeLevelResult,::Type{T}=ComplexF64) where T
    v=NativeLevelArithmetic.approximate(r)
    T<:Real ? T(real(v)) : T(v)
end

"""
    At(q)

Evaluation at a given q, symbolic or numeric, without assuming a root of unity.
"""
struct At{Q} <: EvalTarget
    q::Q
end

"""
    Classical(; exact = false)

The q → 1 limit: ordinary SU(2) recoupling. With `exact = true` the result is a rational/radical value
rather than floating point.
"""
struct Classical <: EvalTarget
    exact::Bool
end
Classical(; exact::Bool = false) = Classical(exact)

target_kwargs(t::Level{T}) where {T} = (; k = t.k, T = T)
target_kwargs(t::Exact) = (; k = t.k, exact = true)
target_kwargs(t::At) = (; q = t.q)
target_kwargs(t::Classical) = (; q = 1, exact = t.exact)

Base.show(io::IO, t::Level{T}) where {T} = print(io, "Level(", t.k, T === Float64 ? "" : "; T = $T", ")")
Base.show(io::IO, t::Exact) = print(io, "Exact(", t.k, ")")
Base.show(io::IO, t::At) = print(io, "At(", t.q, ")")
Base.show(io::IO, t::Classical) = print(io, "Classical(", t.exact ? "; exact = true" : "", ")")

for f in (:q6j, :q3j, :fsymbol, :gsymbol, :rmatrix, :tetrahedron, :theta_value, :qdim, :qeval)
    @eval function $f(t::EvalTarget, args...; kw...)
        t isa Symbolic && throw(ArgumentError("Symbolic() does not accept evaluation keywords"))
        fixed = target_kwargs(t)
        any(key -> key in (:k,:q) || haskey(fixed,key),keys(kw)) &&
            throw(ArgumentError("evaluation target conflicts with explicit evaluation keywords"))
        return $f(args...; fixed...,kw...)
    end
end

"""
    Symbolic()

Request a parameter-independent representation instead of evaluation. Recoupling symbols
return a `DCR` built from their factorial rule; `qdim` and `theta_value` return a
`CyclotomicMonomial`, and `rmatrix` returns a `QPhase`. No evaluation keywords are accepted.
Use `qeval(Symbolic(), rule)` to lower a `FactorialSum` to a DCR.
"""
struct Symbolic <: EvalTarget end
Base.show(io::IO,::Symbolic) = print(io,"Symbolic()")

q6j(::Symbolic,js::Vararg{Spin,6}) = _factorial_dcr(sixj_sum(doubled(js...)...))
q3j(::Symbolic,j1::Spin,j2::Spin,j3::Spin,m1::Spin,m2::Spin,m3::Spin=-m1-m2) =
    _factorial_dcr(threej_sum(doubled(j1,j2,j3,m1,m2,m3)...))
fsymbol(::Symbolic,js::Vararg{Spin,6}) = _factorial_dcr(fsymbol_sum(doubled(js...)...))
gsymbol(::Symbolic,js::Vararg{Spin,6}) = _factorial_dcr(gsymbol_sum(doubled(js...)...))
tetrahedron(::Symbolic,js::Vararg{Spin,6}) = _factorial_dcr(tetrahedron_sum(doubled(js...)...))
qdim(::Symbolic,j::Spin) = qdim_mono(doubled(j))
theta_value(::Symbolic,js::Vararg{Spin,3}) = theta_mono(doubled(js...)...)
function rmatrix(::Symbolic,js::Vararg{Spin,3})
    Js = doubled(js...)
    return _δ(Js...) ? rmatrix_mono(Js...) : zero(QPhase)
end
qeval(::Symbolic,s::FactorialSum) = _factorial_dcr(_validate_rule(s))
qeval(::Symbolic,s::Union{DCR,CyclotomicMonomial,QPhase}) = s

for (f,nlab) in ((:q6j,:((6,))),(:q3j,:((5,6))),(:fsymbol,:((6,))),(:gsymbol,:((6,))))
    @eval function $f(t::Symbolic,labels::Union{AbstractVector,Base.Generator,Base.Iterators.Filter,Tuple{Any,Vararg{Any}}})
        L = _normalize_labels(labels,$nlab,$(string(f)))
        return DCR[$f(t,l...) for l in L]
    end
end
