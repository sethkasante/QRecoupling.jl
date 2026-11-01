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

The value written exactly, as a [`QForm`](@ref) in quantum integers. No value of q is needed: `Exact()` holds at every
generic q, and used as a number it means its value at q = 1, which [`classical`](@ref) gives exactly. `Exact(k)` is the
same value at the level, folded and simplified, with its zeros found and proved.

Numbers come from [`Level`](@ref) and [`At`](@ref), and exact questions from `iszero` and `==`. At a level, `*`, `+` and
`^` are exact too, and [`radical`](@ref), which needs Nemo, rewrites a value in nested square roots where they exist.
Braiding phases stay a `QPhase`.

```julia
q6j(Exact(), 1, 1, 1, 1, 1, 1)             # ([5] − 1)/[4]!, which is 1/6 at q = 1
q6j(Exact(), 1, 1, 1, 1, 1, 1) == 1//6     # true
q6j(Exact(3), 1, 1, 1, 1, 1, 1)            # −1/[2]²
qcg(Exact(), 1//2, 1//2, 1//2, -1//2, 1)   # q^(−1/2) √(1/[2])
```
"""
struct Exact{K} <: EvalTarget
    k::K
    function Exact(k::K) where {K<:Union{Nothing,Integer,AbstractVector{<:Integer}}}
        all(>=(0), something(k, 0)) || throw(DomainError(k, "level must be nonnegative"))
        return new{K}(k)
    end
end
Exact() = Exact(nothing)

"""
    At(q)

Numerical evaluation at a given real or complex q, without assuming a root of unity.
Use `Symbolic()` instead for a generic expression.
"""
struct At{Q} <: EvalTarget
    q::Q
end

"""
    Classical(; exact = false)

The q → 1 limit: ordinary SU(2) recoupling. With `exact = true` the result is a rational/radical value
rather than floating point, as [`classical`](@ref) gives from an `Exact()` value.
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
Base.show(io::IO, t::Exact) = print(io, "Exact(", t.k === nothing ? "" : t.k, ")")
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

Request a parameter-independent representation instead of evaluation. Recoupling symbols return
a [`SymbolicValue`](@ref) retaining the factorial rule without carrying out the sum.
Display shows the deferred factorial rule in x for every label size; it never carries out the sum.
Use `x_form(value)` to request full expansion, or `phi_form(value)` for cyclotomic factorization.
`qdim` and `theta_value` use the same x-form interface. `rmatrix` and `twist` retain `QPhase`
objects, displayed as algebraic phases over x with their q branch. No evaluation keywords are accepted.
Use `qeval(Symbolic(), rule)` to lower a `FactorialSum` the same way.
The weighted coefficients `q3j` and `qcg` do not yet support this target.
"""
struct Symbolic <: EvalTarget end
Base.show(io::IO,::Symbolic) = print(io,"Symbolic()")

q6j(::Symbolic,js::Vararg{Spin,6}) = SymbolicValue(sixj_sum(doubled(js...)...))
q3j_factorial(::Symbolic,j1::Spin,j2::Spin,j3::Spin,m1::Spin,m2::Spin,m3::Spin=-m1-m2) =
    SymbolicValue(threej_sum(doubled(j1,j2,j3,m1,m2,m3)...))
fsymbol(::Symbolic,js::Vararg{Spin,6}) = SymbolicValue(fsymbol_sum(doubled(js...)...))
gsymbol(::Symbolic,js::Vararg{Spin,6}) = SymbolicValue(gsymbol_sum(doubled(js...)...))
tetrahedron(::Symbolic,js::Vararg{Spin,6}) = SymbolicValue(tetrahedron_sum(doubled(js...)...))
function qdim(::Symbolic,j::Spin)
    J=doubled(j)
    J>=0 || throw(DomainError(j,"spin must be nonnegative"))
    return SymbolicValue(FactorialSum((J+1=>1,J=>-1),false,(),false,Int8(1),0,0))
end
twist(::Symbolic,j::Spin) = twist(j;exact=true)
function theta_value(::Symbolic,js::Vararg{Spin,3})
    A,B,C=doubled(js...)
    _δ(A,B,C) || return SymbolicValue(EMPTY_FACTORIAL_SUM)
    t=(A+B+C)÷2
    pre=(t+1=>1,(A+B-C)÷2=>1,(A-B+C)÷2=>1,(-A+B+C)÷2=>1,A=>-1,B=>-1,C=>-1)
    return SymbolicValue(FactorialSum(pre,false,(),false,Int8(1),0,0))
end
function rmatrix(::Symbolic,js::Vararg{Spin,3})
    Js = doubled(js...)
    return _δ(Js...) ? rmatrix_mono(Js...) : zero(QPhase)
end
qint(::Symbolic,n::Integer,p::Integer=1) = _product_symbolic(_qint_pairs(Int(n),Int(p)))
qfact(::Symbolic,n::Integer,p::Integer=1) = _product_symbolic(_qfact_pairs(Int(n),Int(p)))
qbinomial(::Symbolic,n::Integer,m::Integer) = _product_symbolic(_qbinomial_pairs(Int(n),Int(m)))
qeval(::Symbolic,s::FactorialSum) = SymbolicValue(_validate_rule(s))
qeval(::Symbolic,s::Union{SymbolicValue,DCR,CyclotomicMonomial,QPhase}) = s

for (f,nlab) in ((:q6j,:((6,))),(:q3j_factorial,:((5,6))),(:fsymbol,:((6,))),(:gsymbol,:((6,))))
    @eval function $f(t::Symbolic,labels::Union{AbstractVector,Base.Generator,Base.Iterators.Filter,Tuple{Any,Vararg{Any}}})
        L = _normalize_labels(labels,$nlab,$(string(f)))
        return SymbolicValue[$f(t,l...) for l in L]
    end
end






# The exact target. It reuses the symbol interface for the rule and its admissibility, so a symbol added to
# `symbols.jl` reaches this form with no work here.
_exact_target(f,k::AbstractVector,args...) = [_exact_target(f,kk,args...) for kk in k]

"The form of a rule, or zero at a level where the labels are not admissible."
_exact_form(rule,k,admissible) = k !== nothing && !admissible(k) ? _qform_zero(k) : qform(rule;k=k)

function _exact_target(f,k,args...)
    k === nothing || k isa Integer || throw(ArgumentError("exact level must be an integer"))
    kk = k === nothing ? nothing : Int(k)
    kk === nothing || kk >= 0 || throw(DomainError(kk,"level must be nonnegative"))
    f === qdim && return _exact_form(qdim(Symbolic(),only(args)).rule,kk,l -> doubled(only(args)) <= l)
    f === theta_value && return _exact_form(theta_value(Symbolic(),args...).rule,kk,l -> _qδ(doubled(args...)...,l))
    f === tetrahedron && return _exact_form(tetrahedron_sum(doubled(args...)...),kk,l -> _qδtet(doubled(args...)...,l))
    # phases: a `QPhase`, not a field value
    f === rmatrix && return kk === nothing ? rmatrix(Symbolic(),args...) : rmatrix(args...;k=kk,exact=true)
    f === twist && return kk === nothing ? twist(Symbolic(),only(args)) : twist(only(args);k=kk,exact=true)
    # the bare products carry a rule but no labels
    f === qint && return qform(qint(Symbolic(),args...).rule;k=kk)
    f === qfact && return qform(qfact(Symbolic(),args...).rule;k=kk)
    f === qbinomial && return qform(qbinomial(Symbolic(),args...).rule;k=kk)
    if f === qeval
        length(args)==1 || throw(ArgumentError("qeval requires one expression"))
        s=only(args)
        s isa SymbolicValue && return qform(s.rule;k=kk)
        s isa FactorialSum && return qform(_validate_rule(s);k=kk)
        kk === nothing && return qeval(s;q=1,exact=true)
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
