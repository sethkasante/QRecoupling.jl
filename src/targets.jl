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
  `x = q + q⁻¹ = 2cos(π/(k+2))`, at half the degree of the cyclotomic field, printed through the display
  ladder — rational, a surd, nested radicals, or a polynomial in `x`, whichever is exact and short. It
  is closed under `*`, `inv` and `^`, and sums of it are [`ExactXSum`](@ref).
* `form=:canonical` returns canonical cyclotomic coefficients in ℚ(ζ₂ₕ). This was the default before the
  real basis existed; it is what the `exact = true` keyword still returns, and it carries
  `ExactLevelFraction`, the deferred form and the braiding phases.
* `form=:deferred` retains numerator/denominator pairs where the cleared kernel applies, avoiding field
  division until `canonicalize_exact` or numerical conversion.

All three use deterministic arithmetic; modular screening alone never establishes an exact zero.

`form=:deferred` supports scalar/collection recoupling symbols, level sweeps, and `qeval` of factorial
rules, DCRs and cyclotomic monomials. Exact braiding phases retain their existing `QPhase` form.

```julia
q6j(Exact(3), 1, 1, 1, 1, 1, 1)                    # (√5 − 3)/2, the Fibonacci level
q6j(Exact(3; form = :canonical), 1, 1, 1, 1, 1, 1) # the same value in ℚ(ζ₁₀)
```
"""
struct Exact{K} <: EvalTarget
    k::K
    form::Symbol
    function Exact(k::K; form::Symbol=:x) where K
        form in (:canonical,:deferred,:x) ||
            throw(ArgumentError("exact form must be :x, :canonical or :deferred"))
        new{K}(k,form)
    end
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
Base.show(io::IO, t::Exact) = print(io, "Exact(", t.k, t.form === :x ? "" : "; form=:$(t.form)", ")")
Base.show(io::IO, t::At) = print(io, "At(", t.q, ")")
Base.show(io::IO, t::Classical) = print(io, "Classical(", t.exact ? "; exact = true" : "", ")")

for f in (:q6j, :q3j, :fsymbol, :gsymbol, :rmatrix, :tetrahedron, :theta_value, :qdim, :qeval, :twist)
    @eval function $f(t::EvalTarget, args...; kw...)
        t isa Symbolic && throw(ArgumentError("Symbolic() does not accept evaluation keywords"))
        if t isa Exact && t.form === :deferred
            isempty(kw) || throw(ArgumentError("deferred exact targets do not accept evaluation keywords"))
            _deprecated_cyclotomic()
            return _deferred_target($f,t.k,args...)
        end
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
sum. Display uses a simplified x-form for small rules and a deferred factorial sum in x for larger ones.
Use `xvalue(value)` to request full expansion, or `value.dcr` for a lazily constructed compatibility DCR.
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
qeval(::Symbolic,s::FactorialSum) = SymbolicValue(_validate_rule(s))
qeval(::Symbolic,s::Union{SymbolicValue,DCR,CyclotomicMonomial,QPhase}) = s

for (f,nlab) in ((:q6j,:((6,))),(:q3j,:((5,6))),(:fsymbol,:((6,))),(:gsymbol,:((6,))))
    @eval function $f(t::Symbolic,labels::Union{AbstractVector,Base.Generator,Base.Iterators.Filter,Tuple{Any,Vararg{Any}}})
        L = _normalize_labels(labels,$nlab,$(string(f)))
        return SymbolicValue[$f(t,l...) for l in L]
    end
end


# Deferred targets share rule/admissibility definitions with the structural queries.
function _deferred_target(f,k::AbstractVector,args...)
    return [_deferred_target(f,kk,args...) for kk in k]
end
function _deferred_target(f,k,args...)
    k isa Integer || throw(ArgumentError("exact level must be an integer"))
    kk=Int(k); kk>=0 || throw(DomainError(kk,"level must be nonnegative"))
    if f === rmatrix
        return rmatrix(args...;k=kk,exact=true)
    elseif f === twist
        return twist(args...;k=kk,exact=true)          # a phase, exact at any level
    elseif f === qeval
        length(args)==1 || throw(ArgumentError("qeval requires one expression"))
        s=only(args)
        s isa SymbolicValue && return _project_deferred(s.dcr,kk;rule=s.rule)
        if s isa FactorialSum
            _validate_rule(s)
            return _project_deferred(_factorial_dcr(s),kk;rule=s)
        elseif s isa DCR
            return _project_deferred(s,kk)
        elseif s isa CyclotomicMonomial
            return ExactLevelFraction(project_exact(s,kk))
        end
        throw(ArgumentError("unsupported deferred exact expression"))
    elseif f === qdim || f === theta_value
        # the deferred form is a cyclotomic fraction, so it takes the monomial route rather than
        # `exact = true`, which now means the real basis
        mono = f === qdim ? qdim_mono(doubled(only(args))) :
               (Js = doubled(args...); _qδ(Js...,kk) ? theta_mono(Js...) : ZERO_MONOMIAL)
        return ExactLevelFraction(qeval(mono;k=kk,exact=true))
    end
    # Match the existing collection convention, without sharing mutable arithmetic between calls.
    if length(args)==1
        L=_normalize_labels(only(args),f === q3j ? (5,6) : (6,),string(f))
        return [_deferred_target(f,kk,l...) for l in L]
    end
    s=f === tetrahedron ? tetrahedron_sum(doubled(args...)...) : _rule_for(f,args...)
    admissible=f === tetrahedron ? _qδtet(doubled(args...)...,kk) : _admissible_at(f,kk,args...)
    if !admissible
        _,z=cyclotomic_field(2(kk+2),"ζ")
        return CompositeExactResult(kk,typeof(ExactLevelFraction(zero(z))))
    end
    return _project_deferred(_factorial_dcr(s),kk;rule=s)
end


"""
    verify_biedenharn_elliott(Exact(k), labels; workspace = nothing)

Whether the Biedenharn–Elliott identity holds exactly at level `k` for nine spin labels.

Verified in the **real basis** ℚ[x]/Ψ_h, at half the degree of the cyclotomic field this used to run in;
the difference of the two sides comes back as an `ExactXSum` keyed by square class, and cancels to an
empty term list, which is a proof that needs neither a norm nor a numeric fallback.

`workspace` is accepted for compatibility and its level is still checked, but the real-basis route owns
no per-level scratch and ignores it. It will be removed with the rest of the cyclotomic layer in v0.5.
"""
function verify_biedenharn_elliott(t::Exact,labels;workspace=nothing)
    t.k isa Integer || throw(ArgumentError("identity verification requires a scalar integer level"))
    k=Int(t.k);k>=0 || throw(DomainError(k,"level must be nonnegative"))
    length(labels)==9 || throw(ArgumentError("expected nine spin labels"))
    if workspace !== nothing
        workspace isa ExactLevelWorkspace || throw(ArgumentError("expected ExactLevelWorkspace"))
        workspace.k == k || throw(ArgumentError("exact workspace level mismatch"))
    end
    return iszero(_verify_be_x(k,doubled(labels...)))
end


"""
Warn once. The cyclotomic carrier ℚ(ζ₂ₕ) is no longer what the package means by an exact level value —
`Exact(k)` and `exact = true` both give the real basis ℚ[x]/Ψ_h now. `form = :canonical` and
`form = :deferred` still work, and still return exactly what they did, so that anything reading
`CompositeExactResult` or `ExactLevelFraction` keeps working while it migrates; they are scheduled for
removal in v0.5, on the same footing as `eager = true`.
"""
_deprecated_cyclotomic() = @warn(
    "the cyclotomic exact forms (`Exact(k; form = :canonical)` and `form = :deferred`) are deprecated; " *
    "`Exact(k)` and `exact = true` return the real basis ℚ[x]/Ψ_h. They will be removed in v0.5.",
    maxlog = 1)

# The cyclotomic target, kept reachable under its own name while it is deprecated. It goes straight to
# `project_exact` rather than through `exact = true`, which now means the real basis.
_canonical_target(f,k::AbstractVector,args...) = [_canonical_target(f,kk,args...) for kk in k]

function _canonical_target(f,k,args...)
    k isa Integer || throw(ArgumentError("exact level must be an integer"))
    kk=Int(k); kk>=0 || throw(DomainError(kk,"level must be nonnegative"))
    f === rmatrix && return rmatrix(args...;k=kk,exact=true)
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
# added to `symbols.jl` reaches this form with no work here, and it mirrors `_deferred_target`'s handling
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
