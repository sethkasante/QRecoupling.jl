# --------------------------------------------
#  Evaluation targets
#
#  `q6j(Level(10), j...)` instead of `q6j(j...; k = 10)`. The keyword form keeps working; targets exist so
#  that a valid request is a single value rather than a combination of keywords, and so that a target can
#  be stored, passed around and swept over.
# --------------------------------------------

"""
    EvalTarget

Where a symbol should be evaluated: [`Level`](@ref), [`Exact`](@ref), [`At`](@ref) or [`Classical`](@ref).
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
    @eval $f(t::EvalTarget, args...; kw...) = $f(args...; target_kwargs(t)..., kw...)
end
