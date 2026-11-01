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
    Exact()
    Exact(k)

One number, exactly. `Exact()` is the value at q = 1, the same as `Classical(exact = true)`: a sign and a rational under a
root, with its own arithmetic. `Exact(k)` is the value at level `k` as a [`QForm`](@ref) in quantum integers, folded and
simplified, with its zeros found and proved. It needs no backend: `iszero`, `==`, `float`, `*`, `+` and `^` are exact, and
[`radical`](@ref), which needs Nemo, rewrites it in nested square roots where they exist.

Numbers come from [`Level`](@ref) and [`At`](@ref), and the formula at generic q from [`Symbolic`](@ref). Braiding phases
stay a `QPhase`.

```julia
q6j(Exact(), 1, 1, 1, 1, 1, 1)                      # 1//6
q6j(Exact(3), 1, 1, 1, 1, 1, 1)                     # −1/[2]²
iszero(q6j(Exact(10), 1, 3//2, 3//2, 2, 3//2, 3//2))   # true, proved
```
"""
struct Exact{K} <: EvalTarget
    k::K
    function Exact(k::K) where {K<:Union{Integer,AbstractVector{<:Integer}}}
        all(>=(0), k) || throw(DomainError(k, "level must be nonnegative"))
        return new{K}(k)
    end
end
Exact() = Classical(exact = true)

"""
    At(q)

Numerical evaluation at a given real or complex q, without assuming a root of unity.
Use `Symbolic()` for the formula at generic q.

Away from |q| = 1 values shrink quickly with the spins. A value below `floatmin(Float64)` rounds to `0.0`, which is not an
exact zero: at q = 0.8 that is about half of random labels with spins up to 300. `T = BigFloat` gives the value.
"""
struct At{Q} <: EvalTarget
    q::Q
end

"""
    Classical(; exact = false)

The q → 1 limit: ordinary SU(2) recoupling. With `exact = true` the result is a rational/radical value
rather than floating point. `Exact()` is the short name for that.
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

for f in (:q6j, :q3j_factorial, :fsymbol, :gsymbol, :rmatrix, :tetrahedron, :theta_value, :qdim, :qeval,
          :twist, :qint, :qfact, :qbinomial)
    @eval function $f(t::EvalTarget, args...; kw...)
        t isa Symbolic && throw(ArgumentError("Symbolic() does not accept evaluation keywords"))
        if t isa Exact
            isempty(kw) || throw(ArgumentError("`Exact($(t.k))` already sets the target and takes no keywords; remove " *
                                               join(("`$key`" for key in keys(kw)), ", ")))
            return _exact_target($f,t.k,args...)
        end
        fixed = target_kwargs(t)
        any(key -> key in (:k,:q) || haskey(fixed,key),keys(kw)) &&
            throw(ArgumentError("evaluation target conflicts with explicit evaluation keywords"))
        return $f(args...; fixed...,kw...)
    end
    # `Level{T}`: a type passed as a keyword is a `DataType` field of a named tuple, so the call it reaches
    # is dispatched at run time. The default type is therefore left to the callee's own default.
    @eval function $f(t::Level{T}, args::Vararg{Any,N}; kw...) where {T,N}
        if isempty(kw)
            T === Float64 && return $f(args...; k = t.k)     # an explicit `T` keyword is not seen statically
            return $f(args...; k = t.k, T = T)
        end
        any(key -> key === :k || key === :q || key === :T, keys(kw)) &&
            throw(ArgumentError("evaluation target conflicts with explicit evaluation keywords"))
        return $f(args...; k = t.k, T = T, kw...)
    end
end

"""
    Symbolic()

The formula at generic q, as a [`QForm`](@ref) in quantum integers: a function of q, with no value of its own. `f(q)`
evaluates it at any q that [`At`](@ref) takes, and `f(1)` gives the exact value at q = 1. `==` and `iszero` mean at every q,
and are decided exactly without expanding anything. `float` and arithmetic need a value of q, and say so.

`rmatrix` and `twist` keep `QPhase` objects, and `qeval(Symbolic(), rule)` writes a `FactorialSum` the same way. With Nemo,
`x_form(f)` and `phi_form(f)` expand a form in x = q + 1/q or in q. No evaluation keywords are accepted.

```julia
f = q6j(Symbolic(), 1, 1, 1, 1, 1, 1)       # ([5] − 1)/[4]!
f(0.8), f(1)                                # (0.1704…, 1//6)
qcg(Symbolic(), 1//2, 1//2, 1//2, -1//2, 1)  # q^(−1/2) √(1/[2])
```
"""
struct Symbolic <: EvalTarget end
Base.show(io::IO,::Symbolic) = print(io,"Symbolic()")

q6j(::Symbolic,js::Vararg{Spin,6}) = qform(sixj_sum(doubled(js...)...))
q3j_factorial(::Symbolic,j1::Spin,j2::Spin,j3::Spin,m1::Spin,m2::Spin,m3::Spin=-m1-m2) =
    qform(threej_sum(doubled(j1,j2,j3,m1,m2,m3)...))
fsymbol(::Symbolic,js::Vararg{Spin,6}) = qform(fsymbol_sum(doubled(js...)...))
gsymbol(::Symbolic,js::Vararg{Spin,6}) = qform(gsymbol_sum(doubled(js...)...))
tetrahedron(::Symbolic,js::Vararg{Spin,6}) = qform(tetrahedron_sum(doubled(js...)...))
function qdim(::Symbolic,j::Spin)
    J=doubled(j)
    J>=0 || throw(DomainError(j,"spin must be nonnegative"))
    return qform(FactorialSum((J+1=>1,J=>-1),false,(),false,Int8(1),0,0))
end
function twist(::Symbolic,j::Spin)
    J=doubled(j)
    J>=0 || throw(DomainError(j,"spin must be nonnegative"))
    return QPhase(Int8(1),(J*(J+2))//2)
end
function theta_value(::Symbolic,js::Vararg{Spin,3})
    A,B,C=doubled(js...)
    _δ(A,B,C) || return qform(EMPTY_FACTORIAL_SUM)
    t=(A+B+C)÷2
    pre=(t+1=>1,(A+B-C)÷2=>1,(A-B+C)÷2=>1,(-A+B+C)÷2=>1,A=>-1,B=>-1,C=>-1)
    return qform(FactorialSum(pre,false,(),false,Int8(1),0,0))
end
function rmatrix(::Symbolic,js::Vararg{Spin,3})
    Js = doubled(js...)
    return _δ(Js...) ? rmatrix_mono(Js...) : zero(QPhase)
end
qint(::Symbolic,n::Integer,p::Integer=1) = _product_symbolic(_qint_pairs(Int(n),Int(p)))
qfact(::Symbolic,n::Integer,p::Integer=1) = _product_symbolic(_qfact_pairs(Int(n),Int(p)))
qbinomial(::Symbolic,n::Integer,m::Integer) = _product_symbolic(_qbinomial_pairs(Int(n),Int(m)))
qeval(::Symbolic,s::FactorialSum) = qform(_validate_rule(s))
qeval(::Symbolic,s::Union{DCR,CyclotomicMonomial,QPhase}) = s

for (f,nlab) in ((:q6j,:((6,))),(:q3j_factorial,:((5,6))),(:fsymbol,:((6,))),(:gsymbol,:((6,))))
    @eval function $f(t::Symbolic,labels::Union{AbstractVector,Base.Generator,Base.Iterators.Filter,Tuple{Any,Vararg{Any}}})
        L = _normalize_labels(labels,$nlab,$(string(f)))
        return QForm[$f(t,l...) for l in L]
    end
end

# The exact target at a level. It reuses the symbol interface for the rule and its admissibility, so a symbol added to
# `symbols.jl` reaches this form with no work here.
_exact_target(f,k::AbstractVector,args...) = [_exact_target(f,kk,args...) for kk in k]

"The form of a rule at a level, or zero where the labels are not admissible."
_exact_form(rule,k,admissible) = admissible(k) ? qform(rule;k=k) : _qform_zero(k)

function _exact_target(f,k::Integer,args...)
    kk = Int(k)
    f === qdim && return _exact_form(qdim(Symbolic(),only(args)).rule,kk,l -> doubled(only(args)) <= l)
    f === theta_value && return _exact_form(theta_value(Symbolic(),args...).rule,kk,l -> _qδ(doubled(args...)...,l))
    f === tetrahedron && return _exact_form(tetrahedron_sum(doubled(args...)...),kk,l -> _qδtet(doubled(args...)...,l))
    # phases: a `QPhase`, not a field value
    f === rmatrix && return rmatrix(args...;k=kk,exact=true)
    f === twist && return twist(only(args);k=kk,exact=true)
    # the bare products carry a rule but no labels
    f === qint && return qform(qint(Symbolic(),args...).rule;k=kk)
    f === qfact && return qform(qfact(Symbolic(),args...).rule;k=kk)
    f === qbinomial && return qform(qbinomial(Symbolic(),args...).rule;k=kk)
    if f === qeval
        length(args)==1 || throw(ArgumentError("qeval requires one expression"))
        s=only(args)
        s isa QForm && s.k === nothing && return qform(something(s.rule,EMPTY_FACTORIAL_SUM);k=kk,weight=s.weight)
        s isa FactorialSum && return qform(_validate_rule(s);k=kk)
        throw(ArgumentError("`Exact(k)` evaluates a factorial rule; got a $(typeof(s))"))
    end
    sym=symbol_of(f)
    sym === nothing && throw(ArgumentError(
        "`Exact` has no rule for $(f). The symbols it knows are q6j, q3j, qcg, QRecoupling.q3j_factorial, fsymbol, gsymbol, " *
        "tetrahedron, qdim, theta_value and rmatrix."))
    if length(args)==1 && !(only(args) isa Spin)
        L=_normalize_labels(only(args),(nlabels(sym),),string(f))
        return [_exact_target(f,kk,l...) for l in L]
    end
    return _exact_form(symbol_rule(sym,args...),kk,l -> level_admissible(sym,l,args...))
end
