# ---------------------------------------------------------------------------------
#  Quantum Clebsch–Gordan coefficients and the quantum 3j symbol
#
#  For U_q(sl₂) with K|m⟩ = q^m|m⟩, E|m⟩ = √([j−m][j+m+1]) |m+1⟩ and coproduct Δ(E) = E⊗K + K⁻¹⊗E,
#  Δ(K) = K⊗K, the coupled vectors |j m⟩ = Σ ⟨j₁m₁; j₂m₂|j m⟩_q |j₁m₁⟩⊗|j₂m₂⟩ have coefficients
#
#    ⟨j₁m₁; j₂m₂|j m⟩_q = δ_{m₁+m₂,m} q^{E₀} Δ(j₁j₂j) √([j₁+m₁]![j₁−m₁]![j₂+m₂]![j₂−m₂]![j+m]![j−m]![2j+1])
#        × Σ_z (−1)^z q^{−z(j₁+j₂+j+1)} / ([z]! [j₁+j₂−j−z]! [j₁−m₁−z]! [j₂+m₂−z]! [j−j₂+m₁+z]! [j−j₁−m₂+z]!),
#
#  E₀ = ½(j₁+j₂−j)(j₁+j₂+j+1) + j₁m₂ − j₂m₁, Δ² = [j₁+j₂−j]![j₁−j₂+j]![−j₁+j₂+j]!/[j₁+j₂+j+1]!.
#  (Kirillov–Reshetikhin 1989; written in the package's q, whose [n] = (qⁿ − q⁻ⁿ)/(q − q⁻¹).)
#
#  The factorials are those of the 3j rule; the q-power weights are what a symmetric rule cannot carry (the
#  substituted formula, `q3j_factorial`, is 0.76 from orthogonal at q = 0.8). They are kept as two integers,
#  `w` per summation step and the doubled overall power `e2`, and applied inside the scaled representation.
#  At q = 1 the weights are 1 (classical values bit for bit); at a level the values are complex. The 6j
#  needs no weights because it is invariant under q ↔ q⁻¹.
# ---------------------------------------------------------------------------------

"The q-power weights of the Clebsch–Gordan sum, doubled labels: `(w, e2)`, term z × q^{wz}, value × q^{e2/2}."
@inline function _cg_weights(J1::Int, J2::Int, J::Int, M1::Int, M2::Int)
    w = -((J1 + J2 + J) ÷ 2 + 1)
    e2 = ((J1 + J2 - J) ÷ 2) * ((J1 + J2 + J) ÷ 2 + 1) + (J1 * M2 - J2 * M1) ÷ 2
    return (w, e2)
end

"""
    qcg_rule(J1, M1, J2, M2, J, M) -> (rule, (w, e2))

The quantum Clebsch–Gordan coefficient as a factorial rule with its q-power weights, from doubled labels.
The rule is the 3j rule with the dimension [2j+1] under the root and without the Wigner phase; inadmissible
labels, including `M ≠ M1 + M2`, give an empty rule.
"""
function qcg_rule(J1::Int, M1::Int, J2::Int, M2::Int, J::Int, M::Int)
    s = with_dimensions(threej_sum(J1, J2, J, M1, M2, -M), J)
    # `threej_sum` carries (−1)^{j₁−j₂+m}; the coefficient has no phase in front of its sum.
    isodd((J1 - J2 + M) ÷ 2) && (s = negated(s))
    return s, _cg_weights(J1, J2, J, M1, M2)
end

"""
    q3j_rule(J1, J2, J3, M1, M2, M3) -> (rule, (w, e2))

The quantum 3j symbol, (−1)^{j₁−j₂−m₃} [2j₃+1]^{−1/2} ⟨j₁m₁; j₂m₂|j₃, −m₃⟩_q, as the 3j rule (which already
carries the phase) with the Clebsch–Gordan weights.
"""
q3j_rule(J1::Int, J2::Int, J3::Int, M1::Int, M2::Int, M3::Int) =
    (threej_sum(J1, J2, J3, M1, M2, M3), _cg_weights(J1, J2, J3, M1, M2))

"q = e^{iπ/(k+2)} in the precision of `T`."
_level_q(::Type{T}, k::Int) where {T} = cispi(one(real(float(T))) / (k + 2))

_weighted_exact_error(f) = ArgumentError(
    "$(f) at a level is complex, with values in ℚ(ζ_{4h}) up to a square root, and has no exact form yet. " *
    "Use `Level(k)` for numerical values, or `QRecoupling.q3j_factorial` for the real factorial-substitution symbol.")

_weighted_symbolic_error(f) = ArgumentError(
    "$(f) carries q-power weights that the factorial-rule display cannot hold yet. " *
    "Use `At(q)` or `Level(k)` for numerical values, or `QRecoupling.q3j_factorial(Symbolic(), ...)`.")

"""
    _level_half_phases(k) -> Vector{ComplexF64}

e^{iπr/(2h)} for r = 0, …, 4h − 1, h = k + 2, each correctly rounded: integer powers of q = e^{iπ/h} at
even r, and q^{e2/2} at any r. Built once per level, from 128-bit values.
"""
const LEVEL_HALF_PHASES = LevelCache{Vector{ComplexF64}}()
_level_half_phases(k::Int) = get_level!(LEVEL_HALF_PHASES, k) do
    h = k + 2
    setprecision(BigFloat, 128) do
        [ComplexF64(cispi(BigFloat(r) / (2h))) for r in 0:4h-1]
    end
end

"e^{iπr/(2h)} in the precision of `R`, rounded once."
_half_phase(::Type{Float64}, k::Int, r::Int) = @inbounds _level_half_phases(k)[mod(r, 4(k + 2)) + 1]
_half_phase(::Type{R}, k::Int, r::Int) where {R} = cispi(R(mod(r, 4(k + 2))) / R(2(k + 2)))

"Term `z` of a rule at a level, Π [a z + b]!^c, from the split factorial tables: mantissa and exponent."
@inline function _level_term(s::FactorialSum, z::Int, tab::QIntTables{Float64})
    mh, ml, e = 1.0, 0.0, 0
    for f in s.fac
        mh, ml, e = _split_mul_dw(mh, ml, e, tab, _arg(f, z), Int(f.c))
        mh < 0x1p-500 && ((mh, ml, e) = _renorm_dw(mh, ml, e))
    end
    return mh + ml, e
end
@inline function _level_term(s::FactorialSum, z::Int, tab::QIntTables{R}) where {R}
    m, e = one(R), 0
    for f in s.fac
        m, e = _split_mul(m, e, tab, _arg(f, z), Int(f.c))
        m, e = _renorm(m, e)
    end
    return m, e
end

"""
    _weighted_level_pass(s, w, e2, k, R) -> (value, relbound, κ)

One pass over a weighted rule at q = e^{iπ/h}, h = k + 2, for labels inside the level's fusion rule. Every
number is taken at the root itself (a rounded `cispi(1/h)` is amplified by q^{E₀}, |E₀| ~ j²): real level
q-integers, exact phases ζ^r with r = w·z mod 2h, and each term built directly from split factorials in a
double word, so the bound is a fixed multiple of κ = Σ|t_z|/|Σ| (about 4u·κ, charged at 6u·κ).
"""
function _weighted_level_pass(s::FactorialSum, w::Int, e2::Int, k::Int, ::Type{R}) where {R}
    tab = qint_tables(R, k)
    acc = _ascaled(zero(Complex{R})); comp = acc
    mass = AnalyticScaled(zero(R), 0)
    wide = R === BigFloat
    r = 2 * mod(w * s.zlo, 2(k + 2))           # on the half-step grid of `_half_phase`
    for z in s.zlo:s.zhi
        m, e = _level_term(s, z, tab)
        s.alternating && isodd(z) && (m = -m)
        y = _ascaled(m * _half_phase(R, k, r), e)
        if wide
            acc = _aadd(acc, y)
        else
            y = _asub(y, comp); next = _aadd(acc, y)
            comp = _asub(_asub(next, acc), y); acc = next
        end
        mass = _aadd(mass, _ascaled(abs(m), e))
        r += 2w
    end
    if R === Float64
        ph, pl, pe = 1.0, 0.0, 0
        for (n, c) in s.pre
            ph, pl, pe = _split_mul_dw(ph, pl, pe, tab, Int(n), Int(c))
            ph, pl, pe = _renorm_dw(ph, pl, pe)
        end
        if s.sqrt_pre
            isodd(pe) && ((ph, pl, pe) = (2ph, 2pl, pe - 1))
            ph, pl = _dw_sqrt(ph, pl); pe ÷= 2
        end
        pre = _ascaled(ph + pl, pe)
    else
        pm, pe = one(R), 0
        for (n, c) in s.pre
            pm, pe = _split_mul(pm, pe, tab, Int(n), Int(c))
            pm, pe = _renorm(pm, pe)
        end
        pre = AnalyticScaled(pm, pe)
        s.sqrt_pre && (pre = _asqrt(pre))
    end
    v = _amul(_amul(AnalyticScaled(_half_phase(R, k, e2), 0), AnalyticScaled(complex(pre.m), pre.e)), acc)
    s.sign0 < 0 && (v = _aneg(v))
    iszero(acc.m) && return v, R(Inf), R(Inf)
    κ = _avalue(_adiv(mass, _aabs(acc)))
    return v, (6κ + 16) * eps(R), κ
end

"""
Prove that the weighted sum vanishes at a level, after two inconclusive precision passes: a nonzero residue
modulo a prime excludes a zero; a zero residue is confirmed in the cyclotomic field.
"""
function _weighted_level_zero(s::FactorialSum, w::Int, k::Int)
    tab = level_zero_table(k)
    m = tab.m; order = 2(k + 2)
    r = to_mont(m, root_of_unity(m.p, order))
    step = mont_pow(m, r, mod(w, order))
    phase = mont_pow(m, r, mod(w * s.zlo, order))
    acc = UInt64(0)
    N = 0
    for z in s.zlo:s.zhi
        term = phase
        for f in s.fac
            n = _arg(f, z)
            0 <= n <= k + 1 || return false
            N = max(N, n)
            b = f.c > 0 ? tab.fact[n+1, 1] : tab.invfact[n+1, 1]
            for _ in 1:abs(f.c); term = mont_mul(m, term, b); end
        end
        acc = s.alternating && isodd(z) ? mont_sub(m, acc, term) : mont_add(m, acc, term)
        phase = mont_mul(m, phase, step)
    end
    iszero(acc) || return false
    # A modular zero alone is not a proof. This rare path uses exact arithmetic;
    # the square-root prefactor is not needed for deciding whether the sum is zero.
    _, ζ = cyclotomic_field(order, "ζ")
    facts = Vector{typeof(ζ)}(undef, N + 1); facts[1] = one(ζ)
    qi = inv(ζ); den = ζ - qi
    for n in 1:N
        facts[n+1] = facts[n] * divexact(ζ^n - qi^n, den)
    end
    step_exact = ζ^mod(w, order)
    phase_exact = ζ^mod(w * s.zlo, order)
    total = zero(ζ)
    for z in s.zlo:s.zhi
        term = phase_exact
        for f in s.fac; term *= facts[_arg(f, z)+1]^Int(f.c); end
        total += s.alternating && isodd(z) ? -term : term
        phase_exact *= step_exact
    end
    return iszero(total)
end

"""
Level value of a weighted rule: the Float64 pass when its bound certifies it,
then arbitrary precision at doubling widths. Numerical cancellation triggers
an algebraic zero test; agreement near zero at two widths is not a proof.
"""
function _weighted_level_value(s::FactorialSum, w::Int, e2::Int, k::Int, ::Type{T}; near_edge = nothing) where {T}
    R = real(float(T))
    if R !== BigFloat
        v, relb, _ = _weighted_level_pass(s, w, e2, k, R)
        isfinite(relb) && relb <= RTOL_PLAIN && return Complex{R}(_avalue(v))
        # near-edge tier: the coefficient from its Casimir column, O(column) and free of the sum's cancellation
        if near_edge !== nothing
            vr = near_edge()
            vr === nothing || return Complex{R}(vr)
        end
    end
    tol = R === BigFloat ? eps(BigFloat) * 64 : RTOL_PLAIN / 8
    bits = R === BigFloat ? precision(BigFloat) + 32 : 128
    cancelled = false
    zero_checked = false
    for _ in 1:8
        v, relb, κ = setprecision(BigFloat, bits) do
            _weighted_level_pass(s, w, e2, k, BigFloat)
        end
        isfinite(relb) && relb <= tol && return Complex{R}(_avalue(v))
        # A nonzero sum can stay below the noise floor at several widths. Never
        # replace it by zero without an exact check of the q-weighted sum.
        small = !isfinite(κ) || exponent(κ) >= bits - 64
        if small && cancelled && !zero_checked
            _weighted_level_zero(s, w, k) && return zero(Complex{R})
            zero_checked = true
        end
        cancelled = small
        bits *= 2
    end
    throw(ErrorException("level evaluation of a weighted rule did not converge"))
end

"Evaluate a weighted rule at the requested target; `admissible(k)` gives the level's fusion rule."
@inline function _weighted_value(f, s::FactorialSum, weight::Tuple{Int,Int}, admissible::A, k, q, exact::Bool,
                         ::Type{T}; workspace = nothing, near_edge::N = nothing) where {A,T,N}
    if !isnothing(k)
        exact && throw(_weighted_exact_error(f))
        kk = Int(k)
        qk = _level_q(T, kk)
        admissible(kk) || return zero(promote_type(T, typeof(qk)))
        is_empty_sum(s) && return zero(promote_type(T, typeof(qk)))
        return _weighted_level_value(s, weight..., kk, T;
                                     near_edge = near_edge === nothing ? nothing : () -> near_edge(nothing, kk))
    elseif _is_classical(q)
        return exact ? classical_exact(s) : classical_value(s, T; workspace = workspace)
    end
    ne = (near_edge !== nothing && q isa Real && !(q isa BigFloat) && q > 0) ? () -> near_edge(q, nothing) : nothing
    return analytic_value(s, q, T; workspace = workspace, weight = weight, near_edge = ne)
end

"""
    qcg(j1, m1, j2, m2, j, m = m1 + m2; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

The quantum Clebsch–Gordan coefficient ⟨j₁m₁; j₂m₂|j m⟩_q of U_q(sl₂), with K|m⟩ = q^m|m⟩ and coproduct
Δ(E) = E⊗K + K⁻¹⊗E; the default is its classical value, the ordinary Clebsch–Gordan coefficient. The
argument order follows `clebschgordan` in WignerSymbols.jl.

Use `q` for a real or complex parameter and `k` for the level, q = e^{iπ/(k+2)} (mutually exclusive). At real
q > 0 the coefficients form an orthogonal matrix in (j, m₁) at fixed m; at complex q, including a level,
orthogonality is bilinear, Cᵀ C = 1 (at a level, for complete sectors). Negative real q is evaluated as
`complex(q)`. At a level, labels outside the fusion rule j₁ + j₂ + j ≤ k give zero. `exact = true` and
`Exact()` give the exact classical value; exact level values and `Symbolic()` are not available yet.
[`qcg_matrix`](@ref) and [`qcg_row`](@ref) build whole matrices and rows by recurrence.

Square roots follow the balanced-factor convention: `[n] = Π_{d|n,d>1} Ψ_d(q)`, `Ψ_d(q) = q^{−φ(d)} Φ_d(q²)`,
each `Ψ_d` rooted separately; Δ(F) = F⊗K + K⁻¹⊗F, with F the transpose of E.

The value is the 3j rule with term z weighted by q^{−z(j₁+j₂+j+1)} and the whole by
q^{½(j₁+j₂−j)(j₁+j₂+j+1) + j₁m₂ − j₂m₁}, evaluated by the same scaled, escalating kernel. See also [`q3j`](@ref).
"""
function qcg(j1::Spin, m1::Spin, j2::Spin, m2::Spin, j::Spin, m::Spin = m1 + m2;
             k = nothing, q = nothing, exact::Bool = false, T::Type{TT} = Float64,
             workspace = nothing) where {TT}
    q = _evaluation_q(k, q, exact)
    k isa AbstractVector && return [qcg(j1, m1, j2, m2, j, m; k = kk, exact = exact, T = T) for kk in k]
    J1, M1, J2, M2, J, M = doubled(j1, m1, j2, m2, j, m)
    s, weight = qcg_rule(J1, M1, J2, M2, J, M)
    # a coupling coefficient is one entry of its Casimir column: the near-edge tier (`qcg_columns.jl`)
    ne = is_empty_sum(s) ? nothing : (qq, kk) -> _cg_entry(qq, kk, J1, M1, J2, M2, J; workspace=workspace)
    return _weighted_value(qcg, s, weight, kk -> _qδ(J1, J2, J, kk), k, q, exact, T;
                           workspace = workspace, near_edge = ne)
end

"""
    q3j(j1, j2, j3, m1, m2, m3 = -m1-m2; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

The quantum 3j symbol of U_q(sl₂),

    (j₁ j₂ j₃; m₁ m₂ m₃)_q = (−1)^{j₁−j₂−m₃} [2j₃+1]^{−1/2} ⟨j₁m₁; j₂m₂|j₃, −m₃⟩_q,

with the Clebsch–Gordan coefficient of [`qcg`](@ref). The default is its classical value, the Wigner 3j
symbol, which it reproduces bit for bit. Use `q` for a real or complex parameter and `k` for the level
q = e^{iπ/(k+2)} (mutually exclusive). At q ≠ 1 the symbol is complex in general and is related to its
relabellings by q-powers: (j₂ j₁ j₃; m₂ m₁ m₃)_q = (−1)^{j₁+j₂+j₃} (j₁ j₂ j₃; m₁ m₂ m₃)_{1/q}, the same for
m → −m, and (j₂ j₃ j₁; m₂ m₃ m₁)_q = q^{m₂} (j₁ j₂ j₃; m₁ m₂ m₃)_q.

Until v0.4.0 `q3j` was the classical formula with quantum factorials substituted and no q-power weights.
That function is [`QRecoupling.q3j_factorial`](@ref), not exported; the two agree at q = 1 and differ at every other q.
"""
function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1 - m2;
             k = nothing, q = nothing, exact::Bool = false, T::Type{TT} = Float64,
             workspace = nothing) where {TT}
    q = _evaluation_q(k, q, exact)
    k isa AbstractVector && return [q3j(j1, j2, j3, m1, m2, m3; k = kk, exact = exact, T = T) for kk in k]
    J1, J2, J3, M1, M2, M3 = doubled(j1, j2, j3, m1, m2, m3)
    s, weight = q3j_rule(J1, J2, J3, M1, M2, M3)
    # near-edge tier through the coupling coefficient: (−1)^{j1−j2−m3} [2j3+1]^{−1/2} ⟨j1m1; j2m2|j3, −m3⟩
    ne = is_empty_sum(s) ? nothing : (qq, kk) -> begin
        c = _cg_entry(qq, kk, J1, M1, J2, M2, J3; workspace=workspace)
        c === nothing && return nothing
        d = kk === nothing ? qint(J3 + 1; q = qq) : qint(J3 + 1; k = kk)
        return (iseven((J1 - J2 - M3) ÷ 2) ? c : -c) / sqrt(d)
    end
    return _weighted_value(q3j, s, weight, kk -> _qδ(J1, J2, J3, kk), k, q, exact, T;
                           workspace = workspace, near_edge = ne)
end

for f in (:q3j, :qcg)
    @eval function $f(t::EvalTarget, args...; kw...)
        t isa Symbolic && throw(_weighted_symbolic_error($(string(f))))
        t isa Exact && throw(_weighted_exact_error($(string(f))))
        fixed = target_kwargs(t)
        any(key -> key in (:k, :q) || haskey(fixed, key), keys(kw)) &&
            throw(ArgumentError("evaluation target conflicts with explicit evaluation keywords"))
        return $f(args...; fixed..., kw...)
    end
    # A collection of label tuples (5 or 6 labels each). One workspace per worker keeps the tables of the
    # shared q; the results are complex off the positive real axis, so the element type is taken from them.
    @eval function $f(labels::Union{AbstractVector,Base.Generator,Base.Iterators.Filter,Tuple{Any,Vararg{Any}}};
                      k = nothing, q = nothing, exact::Bool = false, T::Type = Float64, threads = nothing)
        L = _normalize_labels(labels, (5, 6), $(string(f)))
        k isa AbstractVector && return [$f(l...; k = kk, q = q, exact = exact, T = T) for l in L, kk in k]
        isempty(L) && return Any[]
        n = length(L)
        first_value = $f(first(L)...; k = k, q = q, exact = exact, T = T)
        out = Vector{typeof(first_value)}(undef, n)
        out[1] = first_value
        _run(n - 1, threads) do rng
            work = EvaluationWorkspace()
            for i in rng
                out[i + 1] = $f(L[i + 1]...; k = k, q = q, exact = exact, T = T, workspace = work)
            end
        end
        return out
    end
end
