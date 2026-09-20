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


"""
Level-k (`q === nothing`) or classical (`q == 1`) value of the symbol with factorial rule `s`.
`admissible()` is the level admissibility test and `dcr()` builds the DCR for the exact and fallback routes;
both run only when needed. The public functions take `T::Type{TT}` so that this call is compiled for the
number type: with a plain `T::Type` keyword it was a dynamic dispatch costing ~100 ns and an allocation.
"""
@inline function _rule_value(s::FactorialSum, admissible::A, dcr::D, k, q, exact::Bool,
                             ::Type{T}; labels=nothing, family=nothing, workspace=nothing) where {A,D,T}
    isnothing(q) || return classical_value(s,T; labels=labels, workspace=workspace)
    return _level_value(s, admissible(), dcr, k,exact,T; labels=labels, family=family, workspace=workspace)
end

"Warn once that the eager route is on its way out."
_deprecated_eager() = @warn("`eager = true` is deprecated: the default level path is faster and more " *
                            "accurate, and the keyword will be removed in v0.5.", maxlog = 1)

"""
    q6j(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, eager::Bool=false, T=Float64)

Returns the Quantum 6j-symbol.
- If `k` and `q` are omitted, returns the abstract `DCR` object.
- If `k` is provided, projects to the root of unity (Float64 by default, Cyclotomic if `exact=true`).
  Zeros and poles are decided from the symbol's factorial rule before any arithmetic, and exact zeros
  are returned as zero.
- If `q=1`, computes the exact classical Ponzano-Regge limit.
- If `eager=true`, bypasses DCR construction for raw speed (only valid for root of unity `k`).
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin;
             k=nothing, q=nothing, exact::Bool=false, eager::Bool=false, T::Type{TT}=Float64, workspace=nothing) where {TT}

    # --- a range of levels is a sweep ---
    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps manage their own scratch"))
    k isa AbstractVector &&
        return _sweep_levels(q6j, (j1, j2, j3, j4, j5, j6), k; q=q, exact=exact, T=T, threads=nothing)

    # --- level k and q = 1: factorial rule (symmetric by construction, no canonical form needed) ---
    if !eager && ((!isnothing(k) && isnothing(q)) || (_is_classical(q) && !exact))
        Jd = doubled(j1, j2, j3, j4, j5, j6)
        return _rule_value(sixj_sum(Jd...), () -> _qδtet(Jd..., k),
                           () -> q6j_dcr(canonical_spins(j1, j2, j3, j4, j5, j6)...), k, q, exact, T;
                           labels=Jd, family=Val(:sixj), workspace=workspace)
    end

    js = canonical_spins(j1, j2, j3, j4, j5, j6)

    # --- eager evaluation circuit ---
    if eager && !isnothing(k) && isnothing(q)
        _deprecated_eager()
        return exact ? q6j_exact(js..., k) : q6j_direct(js..., k, T)
    end

    # --- DCR construction ---
    dcr = q6j_dcr(js...)

    # --- return raw graph if no target is specified ---
    isnothing(q) && return dcr

    # --- DCR Projections ---
    return qeval(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    q3j(j1, j2, j3, m1, m2, m3; k=nothing, q=nothing, exact::Bool=false, eager::Bool=false, T=Float64)
Returns the quantum Wigner 3j-symbol.
"""
function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin=-m1-m2;
             k=nothing, q=nothing, exact::Bool=false, eager::Bool=false, T::Type{TT}=Float64, workspace=nothing) where {TT}

    Js = doubled(j1, j2, j3, m1, m2, m3)

    # --- a range of levels is a sweep ---
    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps manage their own scratch"))
    k isa AbstractVector &&
        return _sweep_levels(q3j, (j1, j2, j3, m1, m2, m3), k; q=q, exact=exact, T=T, threads=nothing)

    # --- eager evaluation circuit ---
    if eager && !isnothing(k) && isnothing(q)
        _deprecated_eager()
        return exact ? q3j_exact(Js..., k) : q3j_direct(Js..., k, T)
    end

    # --- level k and q = 1: factorial rule ---
    if (!isnothing(k) && isnothing(q)) || (_is_classical(q) && !exact)
        return _rule_value(threej_sum(Js...), () -> _qδ(Js[1], Js[2], Js[3], k), () -> q3j_dcr(Js...),
                           k,q,exact,T; family=Val(:threej),workspace=workspace)
    end

    # --- DCR Construction ---
    dcr = q3j_dcr(Js...)

    # --- Return raw graph if no target is specified ---
    isnothing(q) && return dcr

    # --- DCR Projections ---
    return qeval(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    fsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Unitary crossing matrix element: √([d3][d6]) * {6j}.
"""
function fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin;
                 k=nothing, q=nothing, exact::Bool=false, T::Type{TT}=Float64, workspace=nothing) where {TT}

    # F-symbol not fully symmetric! Just double
    Js = doubled(j1, j2, j3, j4, j5, j6)

    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps manage their own scratch"))
    k isa AbstractVector &&
        return _sweep_levels(fsymbol, (j1, j2, j3, j4, j5, j6), k; q=q, exact=exact, T=T, threads=nothing)

    if (!isnothing(k) && isnothing(q)) || (_is_classical(q) && !exact)
        return _rule_value(fsymbol_sum(Js...), () -> _qδtet(Js..., k), () -> fsymbol_dcr(Js...),k,q,exact,T; family=Val(:f),workspace=workspace)
    end

    dcr = fsymbol_dcr(Js...)
    isnothing(q) && return dcr

    return qeval(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    gsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Tetrahedrally symmetric invariant: √(Π[di]) * {6j}.
"""
function gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin;
                 k=nothing, q=nothing, exact::Bool=false, T::Type{TT}=Float64, workspace=nothing) where {TT}

    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps manage their own scratch"))
    k isa AbstractVector &&
        return _sweep_levels(gsymbol, (j1, j2, j3, j4, j5, j6), k; q=q, exact=exact, T=T, threads=nothing)

    # --- level k and q = 1: factorial rule (symmetric by construction, no canonical form needed) ---
    if (!isnothing(k) && isnothing(q)) || (_is_classical(q) && !exact)
        Jd = doubled(j1, j2, j3, j4, j5, j6)
        return _rule_value(gsymbol_sum(Jd...), () -> _qδtet(Jd..., k),
                           () -> gsymbol_dcr(canonical_spins(j1, j2, j3, j4, j5, j6)...),k,q,exact,T;
                           family=Val(:g),workspace=workspace)
    end

    # G-symbol is fully symmetric.
    Js = canonical_spins(j1, j2, j3, j4, j5, j6)

    dcr = gsymbol_dcr(Js...)
    isnothing(q) && return dcr

    return qeval(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    tetrahedron(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Evaluates the standard closed tetrahedron network.
"""
function tetrahedron(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin;
                     k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)

    # Tetrahedron is symmetric.
    Js = canonical_spins(j1, j2, j3, j4, j5, j6)

    if !isnothing(k) && !_qδtet(Js..., k)
        dcr = ZERO_DCR
    else
        dcr = tetrahedron_dcr(Js...)
    end

    if isnothing(k) && isnothing(q)
        return dcr
    end

    return qeval(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    theta_value(j1, j2, j3; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Value of the Theta-graph.
"""
function theta_value(j1::Spin, j2::Spin, j3::Spin;
                     k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)

    Js = doubled(j1, j2, j3)

    if !isnothing(k) && !_qδ(Js..., k)
        mono = ZERO_MONOMIAL
    else
        mono = theta_mono(Js...)
    end

    # classical limit
    if !isnothing(q) && (q == 1 || q == 1.0)
        return exact ? project_classical_exact(mono) : project_classical(mono, T)
    end

    !isnothing(q) && return project_analytic(mono, q)
    isnothing(k) && return mono

    return exact ? project_exact(mono, k) : project_discrete(mono, k, T)
end


"""
    qdim(j; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Quantum dimension [2j+1]_q.
"""
function qdim(j::Spin; k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)
    J = doubled(j)
    mono = qdim_mono(J)

    if !isnothing(q) && (q == 1 || q == 1.0)
        return exact ? project_classical_exact(mono) : project_classical(mono, T)
    end

    !isnothing(q) && return project_analytic(mono, q)
    isnothing(k) && return mono

    return exact ? project_exact(mono, k) : project_discrete(mono, k, T)
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
