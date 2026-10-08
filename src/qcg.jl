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

The quantum 3j symbol, (-1)^{j₁−j₂−m₃} [2j₃+1]^{−1/2} ⟨j₁m₁; j₂m₂|j₃, −m₃⟩_q, as the 3j rule (which already
carries the phase) with the Clebsch-Gordan weights.
"""
q3j_rule(J1::Int, J2::Int, J3::Int, M1::Int, M2::Int, M3::Int) =
    (threej_sum(J1, J2, J3, M1, M2, M3), _cg_weights(J1, J2, J3, M1, M2))

"q = e^{iπ/(k+2)} in the precision of `T`."
_level_q(::Type{T}, k::Int) where {T} = cispi(one(real(float(T))) / (k + 2))

_weighted_exact_error(f) = ArgumentError(
    "$(f) at a level is complex, with values in ℚ(ζ_{4h}) up to a square root, and has no exact form yet. " *
    "Use `Level(k)` for numerical level values, or `Exact()` for exact classical values.")

_weighted_symbolic_error(f) = ArgumentError(
    "$(f) carries q-power weights that the factorial-rule display cannot hold yet. " *
    "Use `At(q)` or `Level(k)` for numerical values.")

"""
    _level_half_phase_parts(k) -> (hi, lo)

High and low parts of e^{iπr/(2h)} for r = 0, …, 4h − 1, h = k + 2. Integer powers of q = e^{iπ/h}
use even r, and q^{e2/2} can use any r. Compute the first quadrant at 128 bits and obtain the others
by exact quarter turns. Both parts are published together in the level cache.
"""
const LEVEL_HALF_PHASES = LevelCache{Tuple{Vector{ComplexF64},Vector{ComplexF64}}}()
_level_half_phase_parts(k::Int) = get_level!(LEVEL_HALF_PHASES, k) do
    h = k + 2
    hi = Vector{ComplexF64}(undef, 4h)
    lo = similar(hi)
    setprecision(BigFloat, 128) do
        # Evaluate the axes explicitly to retain cispi's signed-zero convention.
        for t in 0:3
            z = cispi(BigFloat(t) / 2)
            hi[t*h+1] = ComplexF64(z)
            lo[t*h+1] = ComplexF64(z - hi[t*h+1])
        end
        for r in 1:h-1
            z = cispi(BigFloat(r) / (2h))
            a = ComplexF64(z); b = ComplexF64(z - a)
            @inbounds begin
                hi[r+1] = a; lo[r+1] = b
                hi[h+r+1] = complex(-imag(a), real(a))
                lo[h+r+1] = complex(-imag(b), real(b))
                hi[2h+r+1] = -a; lo[2h+r+1] = -b
                hi[3h+r+1] = complex(imag(a), -real(a))
                lo[3h+r+1] = complex(imag(b), -real(b))
            end
        end
    end
    return hi, lo
end

"The high parts of the cached half-step phases, rounded once to Float64."
_level_half_phases(k::Int) = first(_level_half_phase_parts(k))

"The low parts of the cached half-step phases, so phase = hi + lo to about u²."
_level_half_phases_lo(k::Int) = last(_level_half_phase_parts(k))

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

"The prefactor Π [n]!^c (rooted if the rule says so) of a level rule, formed in double words."
function _weighted_level_prefactor(s::FactorialSum, tab::QIntTables{Float64})
    ph, pl, pe = 1.0, 0.0, 0
    for (n, c) in s.pre
        ph, pl, pe = _split_mul_dw(ph, pl, pe, tab, Int(n), Int(c))
        ph, pl, pe = _renorm_dw(ph, pl, pe)
    end
    if s.sqrt_pre
        isodd(pe) && ((ph, pl, pe) = (2ph, 2pl, pe - 1))
        ph, pl = _dw_sqrt(ph, pl); pe ÷= 2
    end
    return _ascaled(ph + pl, pe)
end

"""
    _weighted_level_pass_dw(s, w, e2, k) -> (value, relbound, κ)

The compensated evaluation of [`_weighted_level_pass`](@ref) for `Float64` output. The plain pass builds each
term in a double word and rounds it before adding, so its error is about u·κ. Here the term, its phase and
the running sum all stay in double words, and the error is about n·u²·κ for n terms: the tier is accepted
up to κ ≈ 10¹⁵ for typical term counts, where the plain pass stops near 10³. The accumulated sum and
prefactor are rounded for the final scaled product.
"""
function _weighted_level_pass_dw(s::FactorialSum, w::Int, e2::Int, k::Int)
    tab = qint_tables(Float64, k)
    hi, lo = _level_half_phase_parts(k); P = 4(k + 2)
    srh = srl = sih = sil = 0.0             # the sum (real and imaginary double words), in units of 2^E
    mass = 0.0; E = typemin(Int)
    r = 2 * mod(w * s.zlo, 2(k + 2))
    for z in s.zlo:s.zhi
        mh, ml, e = 1.0, 0.0, 0
        for f in s.fac
            mh, ml, e = _split_mul_dw(mh, ml, e, tab, _arg(f, z), Int(f.c))
            mh < 0x1p-500 && ((mh, ml, e) = _renorm_dw(mh, ml, e))
        end
        mh, ml, e = _renorm_dw(mh, ml, e)
        s.alternating && isodd(z) && ((mh, ml) = (-mh, -ml))
        @inbounds c = hi[mod(r, P) + 1]
        @inbounds cl = lo[mod(r, P) + 1]
        trh, trl = _dw_mul(mh, ml, real(c), real(cl))
        tih, til = _dw_mul(mh, ml, imag(c), imag(cl))
        m = abs(mh)
        if e > E
            if E != typemin(Int)            # bring the sum down to the larger exponent (exact scalings)
                d = E - e
                srh = ldexp(srh, d); srl = ldexp(srl, d); sih = ldexp(sih, d); sil = ldexp(sil, d)
                mass = ldexp(mass, d)
            end
            E = e
        elseif e < E
            d = e - E
            trh = ldexp(trh, d); trl = ldexp(trl, d); tih = ldexp(tih, d); til = ldexp(til, d)
            m = ldexp(m, d)
        end
        srh, srl = _dw_add(srh, srl, trh, trl)
        sih, sil = _dw_add(sih, sil, tih, til)
        mass += m
        r += 2w
    end
    accm = complex(srh + srl, sih + sil)
    iszero(accm) && return _ascaled(accm), Inf, Inf
    acc = _ascaled(accm, E)
    pre = _weighted_level_prefactor(s, tab)
    v = _amul(_amul(AnalyticScaled(_half_phase(Float64, k, e2), 0), AnalyticScaled(complex(pre.m), pre.e)), acc)
    s.sign0 < 0 && (v = _aneg(v))
    κ = mass / abs(accm)
    n = s.zhi - s.zlo + 1
    # each term: a product of double-word table entries and one double-word phase; each addition 3u².
    return v, (32 + 2n) * κ * eps(Float64)^2 + 16 * eps(Float64), κ
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
        pre = _weighted_level_prefactor(s, tab)
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
Level value of a weighted rule: try the double-word pass first for Float64 output, then the plain pass,
an optional recurrence and arbitrary precision at doubling widths. Numerical cancellation triggers
an algebraic zero test; agreement near zero at two widths is not a proof.
"""
function _weighted_level_value(s::FactorialSum, w::Int, e2::Int, k::Int, ::Type{T}; near_edge = nothing) where {T}
    R = real(float(T))
    if R !== BigFloat
        # The specialized double-word accumulation is cheaper than the general scaled plain loop.
        # Both passes share their phase table, so this order adds no cold-cache setup.
        if R === Float64
            v, relb, _ = _weighted_level_pass_dw(s, w, e2, k)
            isfinite(relb) && relb <= RTOL_PLAIN && return Complex{R}(_avalue(v))
        end
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
    qcg((j1, m1), (j2, m2), (j, m); kwargs...)

The quantum Clebsch-Gordan coefficient ⟨j₁m₁; j₂m₂|j m⟩_q of U_q(sl₂), with K|m⟩ = q^m|m⟩ and coproduct
Δ(E) = E⊗K + K⁻¹⊗E; the default is its classical value, the ordinary Clebsch-Gordan coefficient. The
arguments are the (jᵢ, mᵢ) pairs in the order of ⟨j₁m₁; j₂m₂|j m⟩, with m = m₁ + m₂ by default, or as
three `(j, m)` tuples.

Use `q` for a real or complex parameter and `k` for the level, q = e^{iπ/(k+2)} (mutually exclusive). At real
q > 0 the coefficients form an orthogonal matrix in (j, m₁) at fixed m; at complex q, including a level,
orthogonality is bilinear, Cᵀ C = 1 (at a level, for complete sectors). Negative real q is evaluated as
`complex(q)`. At a level, labels outside the fusion rule j₁ + j₂ + j ≤ k give zero. `exact = true` and
`Exact()` give the exact classical value; exact level values and `Symbolic()` are not available yet.
[`qcg_matrix`](@ref) and [`qcg_row`](@ref) build whole matrices and rows by recurrence.

Square roots follow the balanced-factor convention: `[n] = Π_{d|n,d>1} Ψ_d(q)`, `Ψ_d(q) = q^{−φ(d)} Φ_d(q²)`,
each `Ψ_d` rooted separately; Δ(F) = F⊗K + K⁻¹⊗F, with F the transpose of E.

The value is the 3j rule with term z weighted by q^{-z(j₁+j₂+j+1)} and the whole by
q^{½(j₁+j₂-j)(j₁+j₂+j+1) + j₁m₂ - j₂m₁}, evaluated by the same scaled, escalating kernel. See also [`q3j`](@ref).
"""
function qcg(j1::Spin, m1::Spin, j2::Spin, m2::Spin, j::Spin, m::Spin = m1 + m2;
             k = nothing, q = nothing, exact::Bool = false, T::Type{TT} = Float64,
             workspace = nothing) where {TT}
    q = _evaluation_q(k, q, exact)
    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps do not accept caller-owned scratch"))
    k isa AbstractVector && return [qcg(j1, m1, j2, m2, j, m; k = kk, exact = exact, T = T) for kk in k]
    J1, M1, J2, M2, J, M = doubled(j1, m1, j2, m2, j, m)
    s, weight = qcg_rule(J1, M1, J2, M2, J, M)
    # a coupling coefficient is one entry of its Casimir column: the near-edge tier (`qcg_columns.jl`)
    ne = is_empty_sum(s) ? nothing : (qq, kk) -> _cg_entry(qq, kk, J1, M1, J2, M2, J; workspace=workspace)
    return _weighted_value(qcg, s, weight, kk -> _qδ(J1, J2, J, kk), k, q, exact, T;
                           workspace = workspace, near_edge = ne)
end

qcg(a::Tuple{Spin,Spin}, b::Tuple{Spin,Spin}, c::Tuple{Spin,Spin}; kw...) = qcg(a..., b..., c...; kw...)

"""
    q3j(j1, j2, j3, m1, m2, m3 = -m1-m2; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

The quantum 3j symbol of U_q(sl₂),

    (j₁ j₂ j₃; m₁ m₂ m₃)_q = (−1)^{j₁−j₂−m₃} [2j₃+1]^{−1/2} ⟨j₁m₁; j₂m₂|j₃, −m₃⟩_q,

with the Clebsch–Gordan coefficient of [`qcg`](@ref). The default is its classical value, the Wigner 3j
symbol, which it reproduces bit for bit. Use `q` for a real or complex parameter and `k` for the level
q = e^{iπ/(k+2)} (mutually exclusive). At q ≠ 1 the symbol is complex in general and is related to its
relabellings by q-powers: (j₂ j₁ j₃; m₂ m₁ m₃)_q = (−1)^{j₁+j₂+j₃} (j₁ j₂ j₃; m₁ m₂ m₃)_{1/q}, the same for
m → −m, and (j₂ j₃ j₁; m₂ m₃ m₁)_q = q^{m₂} (j₁ j₂ j₃; m₁ m₂ m₃)_q.

From v0.5, the deformation-dependent weights are included; classical values are unchanged.
`Exact()` is supported; `Exact(k)` and `Symbolic()` are not available yet.
`eager=true` remains a deprecated alias for the standard evaluator.
"""
function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1 - m2;
             k = nothing, q = nothing, exact::Bool = false, T::Type{TT} = Float64,
             workspace = nothing, eager::Bool = false) where {TT}
    q = _evaluation_q(k, q, exact)
    eager && _deprecated_eager()
    k isa AbstractVector && workspace !== nothing &&
        throw(ArgumentError("workspace is for scalar calls; level sweeps do not accept caller-owned scratch"))
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
        q = _evaluation_q(k,q,exact)
        exact && k !== nothing && throw(_weighted_exact_error($(string(f))))
        L = _normalize_labels(labels, (5, 6), $(string(f)))
        k isa AbstractVector && return [$f(l...; k = kk, q = q, exact = exact, T = T) for l in L, kk in k]
        isempty(L) && return Any[]
        n = length(L)
        first_value = $f(first(L)...; k = k, q = q, exact = exact, T = T)
        out = Vector{typeof(first_value)}(undef, n)
        out[1] = first_value
        _run(n - 1, exact ? 1 : threads) do rng
            work = EvaluationWorkspace()
            for i in rng
                out[i + 1] = $f(L[i + 1]...; k = k, q = q, exact = exact, T = T, workspace = work)
            end
        end
        return out
    end
end
