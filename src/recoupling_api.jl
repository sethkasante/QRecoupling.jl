# --------------------------------------------
# Unified Public API for SU(2) TQFT Kernels
# ---------------------------------------------


"Exact zero in ℚ(ζ_{2(k+2)}), in the same form as `project_exact` returns for a vanishing symbol."
function _exact_zero(k::Int)
    _, ζ = cyclotomic_field(2 * (k + 2), "ζ")
    return CompositeExactResult(k, ZERO_MONOMIAL, zero(ζ))
end

"""
Level-k evaluation of a recoupling symbol from its factorial rule `s`. Inadmissible labels give zero.
Numeric values use the factorial rule directly; exact values short-circuit exact zeros (decided by
valuations and a modular test) and otherwise project the DCR built by `dcr()`.
"""
function _level_value(s::FactorialSum, admissible::Bool, dcr, k::Int, exact::Bool, ::Type{T};
                     labels = nothing, family=nothing, workspace=nothing) where {T}
    k >= 0 || throw(DomainError(k, "level must be nonnegative"))
    admissible || return exact ? _exact_zero(k) : zero(T)
    if exact
        is_zero_at_level(s, k) && return _exact_zero(k)
        return project_exact(dcr(), k)
    end
    return value_at_level(s, k, T; fallback = () -> project_discrete(dcr(), k, T), labels=labels, family=family, workspace=workspace)
end


"Warn once: eager is now a compatibility alias for the standard evaluator."
_deprecated_eager() = @warn("`eager=true` is deprecated and now uses the standard factorial-rule evaluator; remove the keyword. It will be removed in v0.5.", maxlog=1)

"Evaluate one symbol rule, constructing its DCR only for projections that need it."
@inline function _symbol_value(s::FactorialSum, admissible::A, k, q, exact, ::Type{T};
                               labels=nothing, family=nothing, workspace=nothing) where {A,T}
    if !isnothing(k)
        return _level_value(s, admissible(), () -> _factorial_dcr(s), Int(k), exact, T;
                            labels=labels, family=family, workspace=workspace)
    elseif _is_classical(q)
        return exact ? classical_exact(s) :
                       classical_value(s,T; labels=labels,workspace=workspace)
    end
    return analytic_value(s,q;workspace=workspace)
end

"""
    q6j(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

Evaluate the q6j symbol; the default is its classical value.
Use `k` for a root-of-unity level or `q` for an analytic parameter (mutually exclusive).
`exact=true` requests an exact classical or level value. `q6j(Symbolic(), ...)`
constructs a DCR from the factorial rule. Numerical classical/level calls evaluate the rule directly.
`eager=true` is deprecated and uses the same evaluator.
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin;
                  k=nothing, q=nothing, exact::Bool=false, T::Type{TT}=Float64,
                  workspace=nothing, eager::Bool=false) where {TT}
    q = _evaluation_q(k,q,exact)
    eager && _deprecated_eager()
    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps do not accept caller-owned scratch"))
    k isa AbstractVector && return _sweep_levels(q6j, (j1, j2, j3, j4, j5, j6), k;
                                                q=q,exact=exact,T=T,threads=nothing)
    Js = doubled(j1, j2, j3, j4, j5, j6)
    s = sixj_sum(Js...)
    return _symbol_value(s, () -> _qδtet(Js...,Int(k)), k,q,exact,T;
                         family=Val(:sixj),workspace=workspace,labels=Js)
end

"""
    q3j(j1, j2, j3, m1, m2, m3; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

Evaluate the Wigner 3j symbol; the default is its classical value.
Use `k` for a root-of-unity level or `q` for an analytic parameter (mutually exclusive).
`exact=true` requests an exact classical or level value. `q3j(Symbolic(), ...)`
constructs a DCR from the factorial rule. Numerical classical/level calls evaluate the rule directly.
`eager=true` is deprecated and uses the same evaluator.
"""
function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin=-m1-m2;
                  k=nothing, q=nothing, exact::Bool=false, T::Type{TT}=Float64,
                  workspace=nothing, eager::Bool=false) where {TT}
    q = _evaluation_q(k,q,exact)
    eager && _deprecated_eager()
    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps do not accept caller-owned scratch"))
    k isa AbstractVector && return _sweep_levels(q3j, (j1, j2, j3, m1, m2, m3), k;
                                                q=q,exact=exact,T=T,threads=nothing)
    Js = doubled(j1, j2, j3, m1, m2, m3)
    s = threej_sum(Js...)
    return _symbol_value(s, () -> _qδ(Js[1],Js[2],Js[3],Int(k)), k,q,exact,T;
                         family=Val(:threej),workspace=workspace)
end

"""
    fsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

Evaluate (−1)^(j1+j2+j4+j5) √([2j3+1][2j6+1]) {6j}, classically by default.
Use `k` for a root-of-unity level or `q` for an analytic parameter (mutually exclusive).
`exact=true` requests an exact classical or level value. `fsymbol(Symbolic(), ...)`
constructs a DCR from the factorial rule. Numerical classical/level calls evaluate the rule directly.
"""
function fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin;
                  k=nothing, q=nothing, exact::Bool=false, T::Type{TT}=Float64,
                  workspace=nothing) where {TT}
    q = _evaluation_q(k,q,exact)
    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps do not accept caller-owned scratch"))
    k isa AbstractVector && return _sweep_levels(fsymbol, (j1, j2, j3, j4, j5, j6), k;
                                                q=q,exact=exact,T=T,threads=nothing)
    Js = doubled(j1, j2, j3, j4, j5, j6)
    s = fsymbol_sum(Js...)
    return _symbol_value(s, () -> _qδtet(Js...,Int(k)), k,q,exact,T;
                         family=Val(:f),workspace=workspace)
end

"""
    gsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

Evaluate √(Πᵢ[2ji+1]) {6j}, classically by default.
Use `k` for a root-of-unity level or `q` for an analytic parameter (mutually exclusive).
`exact=true` requests an exact classical or level value. `gsymbol(Symbolic(), ...)`
constructs a DCR from the factorial rule. Numerical classical/level calls evaluate the rule directly.
"""
function gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin;
                  k=nothing, q=nothing, exact::Bool=false, T::Type{TT}=Float64,
                  workspace=nothing) where {TT}
    q = _evaluation_q(k,q,exact)
    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps do not accept caller-owned scratch"))
    k isa AbstractVector && return _sweep_levels(gsymbol, (j1, j2, j3, j4, j5, j6), k;
                                                q=q,exact=exact,T=T,threads=nothing)
    Js = doubled(j1, j2, j3, j4, j5, j6)
    s = gsymbol_sum(Js...)
    return _symbol_value(s, () -> _qδtet(Js...,Int(k)), k,q,exact,T;
                         family=Val(:g),workspace=workspace)
end

"""
    tetrahedron(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

Evaluate the closed tetrahedron in the package’s existing normalization, classically by default.
Use `k` for a root-of-unity level or `q` for an analytic parameter (mutually exclusive).
`exact=true` requests an exact classical or level value. `tetrahedron(Symbolic(), ...)`
constructs a DCR from the factorial rule. Numerical classical/level calls evaluate the rule directly.
"""
function tetrahedron(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin;
                  k=nothing, q=nothing, exact::Bool=false, T::Type{TT}=Float64,
                  workspace=nothing) where {TT}
    q = _evaluation_q(k,q,exact)
    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps do not accept caller-owned scratch"))
    k isa AbstractVector && return _sweep_levels(tetrahedron, (j1, j2, j3, j4, j5, j6), k;
                                                q=q,exact=exact,T=T,threads=nothing)
    Js = doubled(j1, j2, j3, j4, j5, j6)
    s = tetrahedron_sum(Js...)
    return _symbol_value(s, () -> _qδtet(Js...,Int(k)), k,q,exact,T;
                         family=nothing,workspace=workspace)
end

"""
    theta_value(j1, j2, j3; k=nothing, q=nothing, exact=false, T=Float64)

Theta-graph value, classical by default. Use `Symbolic()` for the monomial representation.
"""
function theta_value(j1::Spin,j2::Spin,j3::Spin; k=nothing,q=nothing,exact::Bool=false,T::Type=Float64)
    q = _evaluation_q(k,q,exact)
    Js = doubled(j1,j2,j3)
    mono = !isnothing(k) && !_qδ(Js...,Int(k)) ? ZERO_MONOMIAL : theta_mono(Js...)
    return qeval(mono;k=k,q=q,exact=exact,T=T)
end

"""
    qdim(j; k=nothing, q=nothing, exact=false, T=Float64)

Quantum dimension [2j+1], classical by default. Use `qdim(Symbolic(), j)` for a monomial.
"""
function qdim(j::Spin; k=nothing,q=nothing,exact::Bool=false,T::Type=Float64)
    q = _evaluation_q(k,q,exact)
    return qeval(qdim_mono(doubled(j));k=k,q=q,exact=exact,T=T)
end

#---- clear caches ---

clear_numeric_caches!() = (@lock _NUMERIC_MODEL_LOCK empty!(_NUMERIC_MODEL_CACHE); nothing)

function clear_exact_caches!()
    @lock EXACT_PHI_LOCK empty!(EXACT_PHI_CACHE)
    @lock EXACT_MODEL_LOCK empty!(EXACT_MODEL_CACHE)
    return nothing
end

function clear_sieve_caches!()
    @lock ROU_TABLE_LOCK empty!(ROU_TABLE_CACHE)
    @lock QINT_TABLES_LOCK empty!(QINT_TABLES)
    empty!(QINT_F64_TABLES)
    empty!(LEVEL_ZERO_TABLES)
    empty!(CLASSICAL_F64_TABLES)
    empty!(CLASSICAL_MOD_TABLES)
    empty!(CLASSICAL_EXACT_PRIMES)
    empty!(MW3_TABLES)
    empty!(MW4_TABLES)
    @lock CLASSICAL_TABLES_LOCK empty!(CLASSICAL_TABLES)
    return nothing
end

"""
    empty_caches!()

Clear the cached level-k tables (numeric, root-of-unity, q-integer, modular and exact cyclotomic).
Useful for freeing memory in long sessions or before benchmarking. The small classical prime-power
sieve is kept, because it is read without a lock.
"""
function empty_caches!()
    clear_numeric_caches!()
    clear_exact_caches!()
    clear_sieve_caches!()
    return nothing
end
