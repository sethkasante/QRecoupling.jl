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
    Exact(k; form=:x)

Exact level values.

* `form=:x` (**the default**) returns an [`ExactX`](@ref): the value in the real basis
  `x = q + q⁻¹ = 2cos(π/(k+2))`, at half the degree of the cyclotomic field, printed as the polynomial
  in `x` it is stored as. [`radical`](@ref) rewrites it in nested square roots where they exist. It is
  closed under `*`, `inv` and `^`, and sums of it are [`ExactXSum`](@ref).
* `form=:canonical` returns canonical cyclotomic coefficients in ℚ(ζ₂ₕ). This was the default before the
  real basis existed, and it is deprecated: it will be removed with the rest of the cyclotomic layer in
  v0.5.

Both use deterministic arithmetic; modular screening alone never establishes an exact zero. Exact
braiding phases retain their `QPhase` form in either.

```julia
q6j(Exact(3), 1, 1, 1, 1, 1, 1)                    # x − 2 at x = 2cos(π/5)
radical(q6j(Exact(3), 1, 1, 1, 1, 1, 1))           # (√5 − 3)/2, the Fibonacci level
q6j(Exact(3; form = :canonical), 1, 1, 1, 1, 1, 1) # the same value in ℚ(ζ₁₀)
```
"""
struct Exact{K} <: EvalTarget
    k::K
    form::Symbol
    function Exact(k::K; form::Symbol=:x) where K
        form in (:canonical,:x) ||
            throw(ArgumentError("exact form must be :x or :canonical"))
        new{K}(k,form)
    end
end

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
Base.show(io::IO, t::Exact) = print(io, "Exact(", t.k, t.form === :x ? "" : "; form=:$(t.form)", ")")
Base.show(io::IO, t::At) = print(io, "At(", t.q, ")")
Base.show(io::IO, t::Classical) = print(io, "Classical(", t.exact ? "; exact = true" : "", ")")

for f in (:q6j, :q3j, :fsymbol, :gsymbol, :rmatrix, :tetrahedron, :theta_value, :qdim, :qeval,
          :twist, :qint, :qfact, :qbinomial)
    @eval function $f(t::EvalTarget, args...; kw...)
        t isa Symbolic && throw(ArgumentError("Symbolic() does not accept evaluation keywords"))
        if t isa Exact && t.form === :x
            isempty(kw) || throw(ArgumentError("`form = :x` targets do not accept evaluation keywords"))
            return _exact_x_target($f,t.k,args...)
        end
        if t isa Exact && t.form === :canonical
            isempty(kw) || throw(ArgumentError("canonical exact targets do not accept evaluation keywords"))
            _deprecated_cyclotomic()
            return _canonical_target($f,t.k,args...)
        end
        fixed = target_kwargs(t)
        any(key -> key in (:k,:q) || haskey(fixed,key),keys(kw)) &&
            throw(ArgumentError("evaluation target conflicts with explicit evaluation keywords"))
        return $f(args...; fixed...,kw...)
    end
end

"""
    Symbolic()

Request a parameter-independent representation instead of evaluation. Recoupling symbols return
a [`SymbolicValue`](@ref) retaining the factorial rule without constructing a DCR or carrying out the
sum. Display shows the deferred factorial rule in x for every label size; it never carries out the sum.
Use `x_form(value)` to request full expansion, or `value.dcr` for a lazily constructed compatibility DCR.
`qdim` and `theta_value` use the same x-form interface. `rmatrix` and `twist` retain `QPhase`
objects, displayed as algebraic phases over x with their q branch. No evaluation keywords are accepted.
Use `qeval(Symbolic(), rule)` to lower a `FactorialSum` the same way.
"""
struct Symbolic <: EvalTarget end
Base.show(io::IO,::Symbolic) = print(io,"Symbolic()")

q6j(::Symbolic,js::Vararg{Spin,6}) = SymbolicValue(sixj_sum(doubled(js...)...))
q3j(::Symbolic,j1::Spin,j2::Spin,j3::Spin,m1::Spin,m2::Spin,m3::Spin=-m1-m2) =
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

for (f,nlab) in ((:q6j,:((6,))),(:q3j,:((5,6))),(:fsymbol,:((6,))),(:gsymbol,:((6,))))
    @eval function $f(t::Symbolic,labels::Union{AbstractVector,Base.Generator,Base.Iterators.Filter,Tuple{Any,Vararg{Any}}})
        L = _normalize_labels(labels,$nlab,$(string(f)))
        return SymbolicValue[$f(t,l...) for l in L]
    end
end






"""
Warn once. The cyclotomic carrier ℚ(ζ₂ₕ) is no longer what the package means by an exact level value —
`Exact(k)` and `exact = true` both give the real basis ℚ[x]/Ψ_h now. `form = :canonical` and
still works, and still returns exactly what it did, so that anything reading
`CompositeExactResult` keeps working while it migrates; it is scheduled for
removal in v0.5, on the same footing as `eager = true`.
"""
_deprecated_cyclotomic() = @warn(
    "the cyclotomic exact form `Exact(k; form = :canonical)` is deprecated; " *
    "`Exact(k)` and `exact = true` return the real basis ℚ[x]/Ψ_h. They will be removed in v0.5.",
    maxlog = 1)

# The cyclotomic target, kept reachable under its own name while it is deprecated. It goes straight to
# `project_exact` rather than through `exact = true`, which now means the real basis.
_canonical_target(f,k::AbstractVector,args...) = [_canonical_target(f,kk,args...) for kk in k]

function _canonical_target(f,k,args...)
    k isa Integer || throw(ArgumentError("exact level must be an integer"))
    kk=Int(k); kk>=0 || throw(DomainError(kk,"level must be nonnegative"))
    f === rmatrix && return rmatrix(args...;k=kk,exact=true)
    f === twist && return twist(only(args);k=kk,exact=true)
    if f === qint || f === qfact || f === qbinomial
        # the deprecated carrier keeps its monomial route; the real basis is what `Exact(k)` uses
        mono = f === qint ? qint_mono(Int.(args)...) :
               f === qfact ? qfact_mono(Int.(args)...) : qbinomial_mono(Int.(args)...)
        return qeval(mono;k=kk,exact=true)
    end
    if f === qdim || f === theta_value
        mono = f === qdim ? qdim_mono(doubled(only(args))) :
               (Js = doubled(args...); _qδ(Js...,kk) ? theta_mono(Js...) : ZERO_MONOMIAL)
        return qeval(mono;k=kk,exact=true)
    end
    if f === qeval
        length(args)==1 || throw(ArgumentError("qeval requires one expression"))
        s=only(args)
        s isa SymbolicValue && return project_exact(s.dcr,kk;rule=s.rule)
        s isa FactorialSum && (_validate_rule(s); return project_exact(_factorial_dcr(s),kk;rule=s))
        s isa DCR && return project_exact(s,kk)
        s isa CyclotomicMonomial && return project_exact(s,kk)
        throw(ArgumentError("unsupported canonical exact expression"))
    end
    if length(args)==1 && !(only(args) isa Spin)
        L=_normalize_labels(only(args),f === q3j ? (5,6) : (6,),string(f))
        return [_canonical_target(f,kk,l...) for l in L]
    end
    s = f === tetrahedron ? tetrahedron_sum(doubled(args...)...) : _rule_for(f,args...)
    admissible = f === tetrahedron ? _qδtet(doubled(args...)...,kk) : _admissible_at(f,kk,args...)
    admissible || return _cyclotomic_zero(kk)
    return project_exact(_factorial_dcr(s),kk;rule=s)
end

function _cyclotomic_zero(k::Int)
    _, ζ = cyclotomic_field(2*(k+2),"ζ")
    return CompositeExactResult(k, ZERO_MONOMIAL, zero(ζ))
end

# The real-basis target. It reuses the symbol interface for the rule and its admissibility, so a symbol
# added to `symbols.jl` reaches this form with no work here, and it mirrors the canonical route's handling
# of level sweeps and label collections.
_exact_x_target(f,k::AbstractVector,args...) = [_exact_x_target(f,kk,args...) for kk in k]

function _exact_x_target(f,k,args...)
    k isa Integer || throw(ArgumentError("exact level must be an integer"))
    kk=Int(k); kk>=0 || throw(DomainError(kk,"level must be nonnegative"))
    f === qdim && return qdim(only(args);k=kk,exact=true)
    f === theta_value && return theta_value(args...;k=kk,exact=true)
    f === tetrahedron && return (_qδtet(doubled(args...)...,kk) ?
        exact_x(tetrahedron_sum(doubled(args...)...),kk) : zero(ExactX,kk))
    f === rmatrix && return rmatrix(args...;k=kk,exact=true)   # a phase: a `QPhase`, not a field value
    f === twist && return twist(only(args);k=kk,exact=true)
    # the bare products: they carry a rule but no labels, so they answer for themselves
    f === qint && return qint(args...;k=kk,exact=true)
    f === qfact && return qfact(args...;k=kk,exact=true)
    f === qbinomial && return qbinomial(args...;k=kk,exact=true)
    if f === qeval
        length(args)==1 || throw(ArgumentError("qeval requires one expression"))
        s=only(args)
        s isa SymbolicValue && return exact_x(s.rule,kk)
        s isa FactorialSum && (_validate_rule(s); return exact_x(s,kk))
        throw(ArgumentError("`form = :x` evaluates a factorial rule; got a $(typeof(s))"))
    end
    sym=symbol_of(f)
    sym === nothing && throw(ArgumentError(
        "`form = :x` has no rule for $(f). The symbols it knows are q6j, q3j, fsymbol, gsymbol, " *
        "tetrahedron, qdim, theta_value and rmatrix."))
    if length(args)==1 && !(only(args) isa Spin)
        L=_normalize_labels(only(args),(nlabels(sym),),string(f))
        return [_exact_x_target(f,kk,l...) for l in L]
    end
    level_admissible(sym,kk,args...) || return zero(ExactX,kk)
    return exact_x(symbol_rule(sym,args...),kk)
end
