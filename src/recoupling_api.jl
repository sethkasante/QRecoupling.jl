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
Numeric and exact values use the factorial rule directly; `dcr()` supplies the numerical fallback.
A modular zero candidate is not sufficient to short-circuit an exact result.
"""
function _level_value(s::FactorialSum, admissible::Bool, dcr, k::Int, exact::Bool, ::Type{T};
                     labels = nothing, family=nothing, workspace=nothing) where {T}
    k >= 0 || throw(DomainError(k, "level must be nonnegative"))
    admissible || return exact ? zero(ExactX, k) : zero(T)
    # `exact = true` is the real basis, ℚ[x]/Ψ_h — the same value `Exact(k)` returns. The cyclotomic
    # carrier is reachable as `Exact(k; form = :canonical)` and is deprecated.
    exact && return exact_x(s, k)
    return value_at_level(s, k, T; fallback = () -> project_discrete(dcr(), k, T), labels=labels, family=family, workspace=workspace)
end


"Warn once: eager is now a compatibility alias for the standard evaluator."
_deprecated_eager() = @warn("`eager=true` is deprecated and now uses the standard factorial-rule evaluator; remove the keyword. It will be removed in later version", maxlog=1)

"Evaluate one symbol rule, constructing its DCR only for projections that need it."
@inline function _symbol_value(s::FactorialSum, admissible::A, k, q, exact, ::Type{T};
                               labels=nothing, family=nothing, workspace=nothing) where {A,T}
    if !isnothing(k)
        return _level_value(s, admissible(), () -> _factorial_dcr(s), Int(k), exact, T;
                            labels=labels, family=family, workspace=workspace)
    elseif _is_classical(q)
        return exact ? classical_exact(s) :
                       classical_value(s,T; labels=labels,family=family,workspace=workspace)
    end
    # the analytic evaluator reads `labels` as those of a bare 6j symbol
    return analytic_value(s,q,T;workspace=workspace,labels=family === Val(:sixj) ? labels : nothing)
end

"""
    q6j(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

Evaluate the q6j symbol; the default is its classical value.
Use `k` for a root-of-unity level or `q` for an analytic parameter (mutually exclusive).
`exact=true` requests an exact classical or level value. `q6j(Symbolic(), ...)`
retains a factorial rule with bounded x-form display; `x_form` expands it explicitly.
Numerical classical/level calls evaluate the rule directly.
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
    QRecoupling.q3j_factorial(j1, j2, j3, m1, m2, m3; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

The classical Racah formula for the 3j symbol with every factorial replaced by a symmetric quantum
factorial, and no q-power weights. It equals the Wigner 3j symbol at q = 1, and it is real at a level, but
for q ≠ 1 it is *not* the U_q(sl₂) coupling coefficient: the quantum 3j symbol is [`q3j`](@ref). This is
the function that was called `q3j` up to v0.4.0. It is not exported: call it as `QRecoupling.q3j_factorial`.
Use `k` for a root-of-unity level or `q` for an analytic parameter (mutually exclusive).
`exact=true` requests an exact classical or level value. `q3j_factorial(Symbolic(), ...)`
retains a factorial rule with bounded x-form display; `x_form` expands it explicitly.
Numerical classical/level calls evaluate the rule directly.
`eager=true` is deprecated and uses the same evaluator.
"""
function q3j_factorial(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin=-m1-m2;
                  k=nothing, q=nothing, exact::Bool=false, T::Type{TT}=Float64,
                  workspace=nothing, eager::Bool=false) where {TT}
    q = _evaluation_q(k,q,exact)
    eager && _deprecated_eager()
    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps do not accept caller-owned scratch"))
    k isa AbstractVector && return _sweep_levels(q3j_factorial, (j1, j2, j3, m1, m2, m3), k;
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
retains a factorial rule with bounded x-form display; `x_form` expands it explicitly.
Numerical classical/level calls evaluate the rule directly.
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
                         family=Val(:f),workspace=workspace,labels=Js)
end

"""
    gsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

Evaluate √(Πᵢ[2ji+1]) {6j}, classically by default.
Use `k` for a root-of-unity level or `q` for an analytic parameter (mutually exclusive).
`exact=true` requests an exact classical or level value. `gsymbol(Symbolic(), ...)`
retains a factorial rule with bounded x-form display; `x_form` expands it explicitly.
Numerical classical/level calls evaluate the rule directly.
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
                         family=Val(:g),workspace=workspace,labels=Js)
end

"""
    tetrahedron(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

Evaluate the closed tetrahedron in the package’s existing normalization, classically by default.
Use `k` for a root-of-unity level or `q` for an analytic parameter (mutually exclusive).
`exact=true` requests an exact classical or level value. `tetrahedron(Symbolic(), ...)`
retains a factorial rule with bounded x-form display; `x_form` expands it explicitly.
Numerical classical/level calls evaluate the rule directly.
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

Theta-graph value, classical by default. Use `Symbolic()` for a rule-backed x-form.
"""
function theta_value(j1::Spin,j2::Spin,j3::Spin; k=nothing,q=nothing,exact::Bool=false,T::Type=Float64)
    q = _evaluation_q(k,q,exact)
    Js = doubled(j1,j2,j3)
    if exact && !isnothing(k)
        kk = Int(k)
        _qδ(Js...,kk) || return zero(ExactX,kk)
        return _qfact_exactx(_theta_qfacts(Js...),kk)
    end
    mono = !isnothing(k) && !_qδ(Js...,Int(k)) ? ZERO_MONOMIAL : theta_mono(Js...)
    return qeval(mono;k=k,q=q,exact=exact,T=T)
end

"""
The q-factorial content of [`theta_mono`](@ref), as `n => c` pairs, so that the exact route builds the
same value the monomial route does rather than a second opinion about it.
"""
function _theta_qfacts(A::Int,B::Int,C::Int)
    m1 = (A+B-C)÷2; m2 = (A-B+C)÷2; m3 = (-A+B+C)÷2; t = (A+B+C)÷2+1
    return [t => 1, m1 => 1, m2 => 1, m3 => 1, A => -1, B => -1, C => -1]
end

"""
    qdim(j; k=nothing, q=nothing, exact=false, T=Float64)

Quantum dimension [2j+1], classical by default. Use `qdim(Symbolic(), j)` for a rule-backed x-form.
"""
function qdim end

# The docstring sits on the bare `function qdim end` above: a docstring written directly on a
# `Base.@constprop` definition is not attached by Julia 1.10, which broke the documentation build.
# Constant propagation drops the exact branch for ordinary calls; inlining keeps the numerical
# result unboxed even though generic real and complex q share this wrapper.
Base.@constprop :aggressive @inline function qdim(j::Spin; k=nothing,q=nothing,exact::Bool=false,
                                          T::Type{TT}=Float64) where {TT}
    q = _evaluation_q(k,q,exact)
    J = doubled(j)
    J >= 0 || throw(DomainError(j,"spin must be nonnegative"))
    n = J+1
    if !exact && T === Float64 && !(k isa AbstractVector) &&
       (q === nothing || real(typeof(float(q))) !== BigFloat)
        v = _qnumber_float(:int,n,1,k,q)
        v === nothing || return v
        return _product_value(_qint_pairs(n,1),k,q,exact,T)::Union{Float64,ComplexF64}
    end
    return _product_value(_qint_pairs(n,1),k,q,exact,T)
end

#---- clear caches ---

function clear_exact_caches!()
    @lock EXACT_PHI_LOCK empty!(EXACT_PHI_CACHE)
    return nothing
end

function clear_sieve_caches!()
    @lock _SIN_LOCK empty!(_SIN_TABLE)
    @lock ROU_TABLE_LOCK empty!(ROU_TABLE_CACHE)
    @lock QINT_TABLES_LOCK empty!(QINT_TABLES)
    empty!(QINT_F64_TABLES)
    empty!(LEVEL_HALF_PHASES)   # high and low phase parts in qcg.jl
    clear_cg_caches!()          # qcg_columns.jl
    empty!(LEVEL_ZERO_TABLES)
    empty!(LEVEL_PHI)
    empty!(CLASSICAL_F64_TABLES)
    empty!(CLASSICAL_MOD_TABLES)
    empty!(GENERIC_MOD_TABLES)
    empty!(CLASSICAL_EXACT_PRIMES)
    empty!(MW3_TABLES)
    empty!(MW4_TABLES)
    @lock CLASSICAL_TABLES_LOCK empty!(CLASSICAL_TABLES)
    return nothing
end

"""
    empty_caches!()

Clear the cached level-k tables (numeric, root-of-unity, q-integer, modular and exact cyclotomic).
Useful for freeing numerical-table memory in long sessions or before benchmarking. Reusable symbolic
polynomials and number fields are retained, as is the small classical prime-power sieve.
"""
function empty_caches!()
    clear_analytic_caches!()
    clear_exact_caches!()
    clear_sieve_caches!()
    return nothing
end
