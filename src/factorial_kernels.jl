# ---------------------------------------------------------------------------------
#  Shared pieces of the two kernels: the prefactor ∏[n]!^c (with any square root folded in) and the
#  rounding cost of a term or ratio. The kernels differ only in how they traverse a segment.
# ---------------------------------------------------------------------------------

"Prefactor in split form `(mantissa, exponent, relative error bound)`, single word."
@inline function _rule_prefactor(s::FactorialSum, tab::QIntTables{T}) where {T}
    mp = one(T); ep = 0; np = 0
    @inbounds for (n, c) in s.pre
        mp, ep = _split_mul(mp, ep, tab, Int(n), c)
        np += abs(Int(c))
    end
    mp, ep = _renorm(mp, ep)
    rel = _gamma(T, 2np)
    if s.sqrt_pre
        isodd(ep) && (mp *= 2; ep -= 1)
        mp = sqrt(mp); ep = ep ÷ 2
        rel = (rel / 2) + _unit(T)    # halving the exponent is exact; the root adds one rounding
    end
    return mp, ep, rel
end

"The same prefactor as a double word."
@inline function _rule_prefactor_dw(s::FactorialSum, tab::QIntTables{T}) where {T}
    ph = one(T); pl = zero(T); ep = 0
    @inbounds for (n, c) in s.pre
        ph, pl, ep = _split_mul_dw(ph, pl, ep, tab, Int(n), c)
    end
    ph, pl, ep = _renorm_dw(ph, pl, ep)
    if s.sqrt_pre
        isodd(ep) && (ph *= 2; pl *= 2; ep -= 1)
        ph, pl = _dw_sqrt(ph, pl); ep = ep ÷ 2
    end
    return ph, pl, ep
end

"Factorials multiplied into a segment's first term: the roundings it carries."
@inline _first_term_cost(s::FactorialSum) = sum(f -> abs(Int(f.c)), s.fac; init = 0)

"A segment's first term ∏[arg(f,z0)]!^c in split form, single word."
@inline function _first_term(s::FactorialSum, z0::Int, tab::QIntTables{T}) where {T}
    m0 = one(T); e0 = 0
    @inbounds for f in s.fac
        m0, e0 = _split_mul(m0, e0, tab, _arg(f, z0), f.c)
    end
    return _renorm(m0, e0)
end

"The same first term as a double word."
@inline function _first_term_dw(s::FactorialSum, z0::Int, tab::QIntTables{T}) where {T}
    mh = one(T); ml = zero(T); e0 = 0
    @inbounds for f in s.fac
        mh, ml, e0 = _split_mul_dw(mh, ml, e0, tab, _arg(f, z0), f.c)
    end
    return _renorm_dw(mh, ml, e0)
end

# ---------------------------------------------------------------------------------
#  Summation kernels for a factorial rule, with certified error bounds: the plain ratio loop,
#  compensated Horner, the acceptance policy, and the escalation ladder (zero test, K-word, BigFloat).
# ---------------------------------------------------------------------------------

"""
    _sum_at_level(s, segs, k, tab) -> (value, bound)

Σ over the contributing segments in type T by the term ratio, with the prefactor. `bound` is a rigorous
bound on |value − exact|: a term reached after j ratio steps carries at most j(2K+1) relative roundings
(K factorials per ratio), and summation adds at most u Σ|partial sums|. Bookkeeping is in split form,
mantissa · 2^exponent, so no `log` or `exp` is taken.
"""
function _sum_at_level(s::FactorialSum, segs, k::Int, tab::QIntTables{T},
                       buf::B = nothing) where {T,B<:Union{Nothing,Vector}}
    q = tab.q

    ib = 0     # buffer: (rh, rl) per step, segment by segment
    u = _unit(T)
    mp, ep, relpre = _rule_prefactor(s, tab)
    intr = _integer_ratios(s, segs, tab)
    qq, qql = tab.qq, tab.qql
    unit, plan = _ratio_plan(s, length(q))
    rsign = s.alternating ? -one(T) : one(T)
    Kratio = _ratio_cost(s)
    cstep = intr ? 2 : 2Kratio + 1     # roundings per step: a/b and t·r, or K entries + K products
    relfirst = _gamma(T, 2 * _first_term_cost(s))
    big_ = ldexp(one(T), SCALE_BITS)
    am = zero(T); ae = 0; started = false          # the sum, and its error bound in the same frame
    Eacc = zero(T)
    @inbounds for seg in segs
        z0 = first(seg)
        nsteps = length(seg) - 1
        m0, e0 = _first_term(s, z0, tab)
        t = one(T); ssum = one(T); sc = 0
        W = zero(T); SS = zero(T); j = 0           # Σ j|t_j| and Σ|partial sums|: the bound needs only these
        for z in z0:last(seg)-1
            j += 1
            if intr
                a, b = _int_ratio(s, z)
                fa = T(a); fb = T(b)
                r = fa / fb
                if buf !== nothing
                    buf[ib+1] = r; buf[ib+2] = fma(-r, fb, fa) / fb; ib += 2
                end
            elseif unit && buf === nothing
                r = rsign
                for (a, o) in plan
                    r *= qq[a*z+o]
                end
            elseif unit                                         # same, low part kept for later
                r = rsign; rl = zero(T)
                for (a, o) in plan
                    i = a * z + o
                    xh = qq[i]; xl = qql[i]
                    p, π = _two_prod(r, xh)
                    rl = fma(rl, xh, fma(r, xl, π)); r = p
                end
                buf[ib+1] = r; buf[ib+2] = rl; ib += 2
            elseif buf === nothing
                r = _general_ratio(s, z, tab)
            else
                r, rl = _general_ratio_dw(s, z, tab)
                buf[ib+1] = r; buf[ib+2] = rl; ib += 2
            end
            t *= r
            ssum += t
            at = abs(t)
            W = fma(T(j), at, W); SS += abs(ssum)
            if at > big_
                t = ldexp(t, -SCALE_BITS); ssum = ldexp(ssum, -SCALE_BITS)
                W = ldexp(W, -SCALE_BITS); SS = ldexp(SS, -SCALE_BITS)
                sc += SCALE_BITS
            end
        end
        # error of this segment in units of 2^(sc+e0): terms, summation, first term, final product
        g = _gamma(T, cstep * nsteps)
        Eraw = abs(m0) * (cstep * u * W / (1 - g)^2 + u / (1 - u) * SS) + abs(ssum * m0) * (relfirst + u)
        sgn = (s.alternating && isodd(z0)) ? -one(T) : one(T)
        vm, ve = frexp(sgn * ssum * m0)
        ve = iszero(vm) ? sc + e0 : ve + sc + e0
        Em = ldexp(Eraw, sc + e0 - ve)
        if !started
            am, ae, Eacc, started = vm, ve, Em, true
        else
            if ve > ae
                am = ldexp(am, ae - ve); Eacc = ldexp(Eacc, ae - ve); ae = ve
                am += vm; Eacc += Em
            else
                am += ldexp(vm, ve - ae); Eacc += ldexp(Em, ve - ae)
            end
            Eacc += u * abs(am)
        end
        if !iszero(am)
            fr2, ex2 = frexp(am)
            am = fr2; ae += ex2
            Eacc = ldexp(Eacc, -ex2)
        end
    end
    started || return zero(T), zero(T)
    bound = ldexp(Eacc * abs(mp) * (1 + relpre) + abs(am * mp) * (relpre + u), ae + ep) * BOUND_SLACK(T)
    iszero(am) && return zero(T), bound
    return s.sign0 * ldexp(am * mp, ae + ep), bound
end

"Covers the rounding of the bound's own evaluation; the bound itself is computed in working precision."
BOUND_SLACK(::Type{T}) where {T} = one(T) + T(1e-6)

"""
    _digits_lost(value, bound) -> Float64

Decimal digits cancellation has consumed, log₁₀(B/(u|v|)) from the certified bound B. It overstates the
loss by 1.3–3.9 digits, the safe direction for sizing an escalated pass. `Inf` when the bound does not
separate the value from zero.
"""
function _digits_lost(value::T, bound) where {T}
    v = abs(Float64(value)); b = Float64(bound)
    (iszero(v) || !isfinite(b)) && return Inf
    return max(0.0, log10(b / v) - log10(Float64(_unit(T))))
end

"""
    _bound_pessimism(s, segs, tab) -> Float64

The factor `c·n` (roundings charged per ratio step, over every segment) by which the rigorous bound exceeds
the sharp O(u Σ|terms|) error. Only estimates divide it out; callers keep the rigorous bound.
"""
function _bound_pessimism(s::FactorialSum, segs, tab::QIntTables)
    nsteps = sum(seg -> length(seg) - 1, segs; init = 0)
    cstep = _integer_ratios(s, segs, tab) ? 2 : 2 * _ratio_cost(s) + 1
    return max(1.0, Float64(cstep * nsteps))
end

"""
    _surviving_digits(value, bound, pessimism) -> Float64

Decimal digits of `value` that its `bound` leaves standing, with the bound's pessimism removed. An estimate,
used to stop escalation: a pass that lost every digit has `bound/|value| ≈ c·n`, which says nothing about
the cancellation, so sizing from the raw ratio overshoots.
"""
function _surviving_digits(value, bound, pessimism)
    # Logarithms in the value's own type: a BigFloat value below the Float64 range (e.g. 2e-344 at a large
    # level) must not read as zero, or the doubling loop that relies on this estimate never terminates.
    (iszero(value) || !isfinite(value) || !isfinite(bound) || bound < 0) && return -Inf
    iszero(bound) && return Inf
    return log10(pessimism) - Float64(log10(abs(bound)) - log10(abs(value)))
end

"""
    _sum_compensated(s, segs, k, tab) -> (value, bound)

Compensated Horner (Graillat–Langlois–Louvet 2005) over ratios that are products of double-word entries:

    Y_i = 1 + r_i Y_{i+1}  from the top of each segment,  S_segment = t_first · Y_0,

with TwoProd/TwoSum errors propagated in a working-precision correction. Relative error
≤ u + O((c n u)²) κ, as the plain loop in doubled precision; `bound` is its a posteriori form.
"""
function _sum_compensated(s::FactorialSum, segs, k::Int, tab::QIntTables{T},
                          buf::B = nothing) where {T,B<:Union{Nothing,Vector}}
    q = tab.q
    base = 0                                       # buffer offset of the current segment
    u = _unit(T)
    ph, pl, ep = _rule_prefactor_dw(s, tab)
    intr = _integer_ratios(s, segs, tab)
    qq, qql = tab.qq, tab.qql
    unit, plan = _ratio_plan(s, length(q))
    rsign = s.alternating ? -one(T) : one(T)
    big_ = ldexp(one(T), SCALE_BITS)
    ah_ = zero(T); al_ = zero(T); ae = 0; started = false
    Aacc = zero(T); Bacc = zero(T)                 # Σ|terms| and the error bound, in the frame 2^ae
    @inbounds for seg in segs
        z0 = first(seg); z1 = last(seg)
        nseg = z1 - z0
        mh, ml, e0 = _first_term_dw(s, z0, tab)
        y = one(T); c = zero(T); one_s = one(T); A = one(T); EA = zero(T); sc = 0
        for z in (z1 - 1):-1:z0
            if buf !== nothing
                ib = base + 2 * (z - z0)
                rh = buf[ib+1]; rl = buf[ib+2]
            elseif intr
                a, b = _int_ratio(s, z)
                fa = T(a); fb = T(b)
                rh = fa / fb
                rl = fma(-rh, fb, fa) / fb                  # a − rh·b is exact; one rounding on the low part
            elseif unit
                rh = rsign; rl = zero(T)
                for (a, o) in plan
                    i = a * z + o
                    xh = qq[i]; xl = qql[i]
                    p, π = _two_prod(rh, xh)
                    rl = fma(rl, xh, fma(rh, xl, π)); rh = p
                end
            else
                rh, rl = _general_ratio_dw(s, z, tab)
            end
            p, π = _two_prod(rh, y)
            sy, σ = _two_sum(one_s, p)
            c = fma(rh, c, fma(rl, y, π + σ))
            EA = fma(abs(rh), EA, abs(π) + abs(σ) + abs(rl * y))   # sizes of the captured error terms
            y = sy
            A = fma(abs(rh), A, one_s)
            if abs(y) > big_ || A > big_
                y = ldexp(y, -SCALE_BITS); c = ldexp(c, -SCALE_BITS); EA = ldexp(EA, -SCALE_BITS)
                A = ldexp(A, -SCALE_BITS); one_s = ldexp(one_s, -SCALE_BITS); sc += SCALE_BITS
            end
        end
        base += 2 * (z1 - z0)
        yh, yl = _two_sum(y, c)
        vh, vl = _dw_mul(mh, ml, yh, yl)
        if s.alternating && isodd(z0)
            vh = -vh; vl = -vl
        end
        if iszero(vh)
            ve = sc + e0
        else
            vh, vl, ve = _renorm_dw(vh, vl, sc + e0)
        end
        # error of the correction's own evaluation (γ on the captured terms), plus the dropped second-order
        # parts of the ratios (≤ 2K u² per step on every term), in units of 2^(sc+e0)
        Kr = intr ? 1 : _ratio_cost(s)                                # low-part error per ratio: one rounding, or K
        Braw = abs(mh) * (_gamma(T, (2Kr + 4) * max(nseg, 1)) * EA + 2 * Kr * max(nseg, 1) * u^2 * A)
        Am = ldexp(abs(mh) * A, sc + e0 - ve)
        Bm = ldexp(Braw, sc + e0 - ve)
        if !started
            ah_, al_, ae, Aacc, Bacc, started = vh, vl, ve, Am, Bm, true
        else
            if ve > ae
                ah_ = ldexp(ah_, ae - ve); al_ = ldexp(al_, ae - ve)
                Aacc = ldexp(Aacc, ae - ve); Bacc = ldexp(Bacc, ae - ve); ae = ve
            else
                vh = ldexp(vh, ve - ae); vl = ldexp(vl, ve - ae); Am = ldexp(Am, ve - ae); Bm = ldexp(Bm, ve - ae)
            end
            s1, e1 = _two_sum(ah_, vh)
            ah_, al_ = _fast_two_sum(s1, e1 + (al_ + vl))
            Aacc += Am; Bacc += Bm
        end
        if !iszero(ah_)
            ah_, al_, ex2 = _renorm_dw(ah_, al_, 0)
            ae += ex2; Aacc = ldexp(Aacc, -ex2); Bacc = ldexp(Bacc, -ex2)
        end
    end
    started || return zero(T), T(Inf)
    rh_, rl_ = _dw_mul(ah_, al_, ph, pl)
    value = s.sign0 * ldexp(rh_ + rl_, ae + ep)
    # final rounding + double-word products (prefactor, first terms, combination: ≲ 128 u² relative)
    # + the computed bound on the correction, carried through the prefactor
    bound = u * abs(value) + 128 * u^2 * abs(value) + ldexp(Bacc * abs(ph) * (1 + 4u), ae + ep)
    return value, bound * (1 + T(1e-6))
end

"""
Accuracy policy for `Float64` values: keep the plain pass when its certified relative bound is at most
`RTOL_PLAIN`, else the compensated pass at `RTOL_CERTIFIED`, else the zero test and higher precision. The
plain bound is typically 10–1000× pessimistic.
"""
const RTOL_PLAIN = 2.0^-40          # ≈ 9.1e-13 certified
const RTOL_CERTIFIED = 2.0^-44

"Decimal digits a level value in type T must keep; below this the sum is redone at higher precision."
_target_digits(::Type{T}) where {T} = T === BigFloat ? _decimal_digits(T) - 4 : min(12.0, _decimal_digits(T) - 4)

"""
    value_at_level(s, k, T; fallback) -> T

Value at q = e^{iπ/(k+2)}. Poles (a `DomainError`) and vanishing terms are decided from valuations. When
cancellation leaves too few digits, identities or exact confirmation decide zeros and nonzero values
escalate precision. `fallback()` handles factorials outside the level tables. With `labels` (a 6j's doubled
labels), escalation first tries the column recurrence.
"""
function value_at_level(s::FactorialSum, k::Int, ::Type{T}; fallback, labels = nothing, family = nothing, workspace = nothing) where {T}
    v, status, segs = level_pass1(s, k, qint_tables(T, k); family=family, workspace=workspace)
    status === :done && return v
    status === :fallback && return fallback()
    return level_escalate(s, segs, k, T, level_zero_table(k); labels=labels, workspace=workspace)
end

"""
    level_pass1(s, k, tab) -> (value, status, segments)

First pass in the table's own type, lock-free and without `BigFloat`, so safe on many threads. `status` is
`:done`, `:escalate` (cancellation: zero test and higher precision follow), or `:fallback` (a factorial
outside the level tables).
"""
const NO_SEGMENTS = UnitRange{Int}[]          # shared empty result; never mutated

function level_pass1(s::FactorialSum, k::Int, tab::QIntTables{T}; family=nothing, workspace=nothing) where {T}
    is_empty_sum(s) && return (zero(T), :done, NO_SEGMENTS)
    if family !== nothing
        # Only callers that checked level admissibility may select a standard-symbol shortcut.
        segs = (s.zlo:(family === Val(:threej) ? s.zhi : min(s.zhi,k)),)
        v, st = _certified_value(s, segs, k, tab, workspace)
        return (v, st, segs)
    end
    if _valuation_free(s, k + 2)                  # one segment, nothing to classify
        seg = s.zlo:s.zhi
        v, st = _certified_value(s, (seg,), k, tab, workspace)
        return (v, st, st === :done ? NO_SEGMENTS : [seg])
    end
    small = _classify_small(s, k)
    small === nothing || return _pass1_segments(s, k, tab, small...; workspace=workspace)
    st, segs = classify_at_level(s, k)
    st === :pole && throw(DomainError(k, "Topological pole at level k=$k."))
    (st === :empty || st === :zero) && return (zero(T), :done, UnitRange{Int}[])
    _within_level_tables(s, segs, k) || return (zero(T), :fallback, segs)
    v, st = _certified_value(s, segs, k, tab, workspace)
    return (v, st, segs)
end

"level_pass1 from a classification held in a `SegList` (no allocation unless the value escalates)."
function _pass1_segments(s::FactorialSum, k::Int, tab::QIntTables{T}, st::Symbol, segs::SegList; workspace=nothing) where {T}
    st === :pole && throw(DomainError(k, "Topological pole at level k=$k."))
    st === :zero && return (zero(T), :done, NO_SEGMENTS)
    _within_level_tables(s, segs, k) || return (zero(T), :fallback, collect(segs))
    v, st2 = _certified_value(s, segs, k, tab, workspace)
    return (v, st2, st2 === :done ? NO_SEGMENTS : collect(segs))
end

"""
    _certified_value(s, segs, k, tab, workspace) -> (value, :done | :escalate)

Plain pass, then the compensated pass if its bound is too loose; a value is kept when its own bound
certifies the type's target. Thread safe.
"""
function _certified_value(s::FactorialSum, segs, k::Int, tab::QIntTables{T}, workspace=nothing) where {T}
    if T !== Float64                     # a wider type keeps its own digit budget, as before
        v, B = _sum_at_level(s, segs, k, tab)
        kept = _surviving_digits(v, B, _bound_pessimism(s, segs, tab))
        return v, (kept >= _target_digits(T) ? :done : :escalate)
    end
    mode = POLICY[]
    if mode === :compensated_only
        vc, Bc = _sum_compensated(s, segs, k, tab)
        return vc, (_certifies(vc, Bc, RTOL_CERTIFIED) ? :done : :escalate)
    end
    rtol = mode === :strict || mode === :strict_lazy ? RTOL_CERTIFIED : RTOL_PLAIN
    nsteps = sum(seg -> length(seg) - 1, segs; init = 0)
    # A one-step sum recomputes its ratio on fallback: cheaper than allocating a buffer on every call.
    lazy = (mode === :lazy || mode === :strict_lazy) && nsteps >= LAZY_MIN_STEPS
    slot = nothing
    if !lazy
        buf = nothing
    elseif workspace === nothing
        slot, buf = _borrow_ratio_buffer(2 * nsteps)    # no allocation on the common path
    else
        buf = _ratio_buffer(workspace, T, 2 * nsteps)
    end
    v, B = _sum_at_level(s, segs, k, tab, buf)
    if _certifies(v, B, rtol)
        _return_ratio_buffer(slot)
        return v, :done
    end
    # A zero never certifies, so a pairwise-cancelling sum would pay for the compensated pass and the
    # escalation before being recognised. The test is a proof and costs a few comparisons, and it runs
    # only here, after the plain pass has already failed — never on a value that certifies.
    if pairwise_zero(s)
        _return_ratio_buffer(slot)
        return zero(v), :done
    end
    vc, Bc = _sum_compensated(s, segs, k, tab, buf)
    _return_ratio_buffer(slot)
    _certifies(vc, Bc, RTOL_CERTIFIED) && return vc, :done
    return vc, :escalate
end

"Does a finite, nonzero value satisfy the relative error estimate? Invalid bounds are rejected."
@inline _certifies(value, bound, rtol) =
    !iszero(value) && isfinite(value) && isfinite(bound) && 0 <= bound <= rtol * abs(value)

"""
Evaluation policy for `Float64` level and classical values:

- `:lazy` (default): plain pass kept at a certified 2⁻⁴⁰; it stores each ratio's low part, so the
  compensated fallback reuses them (15–28% cheaper fallback, no cost on easy symbols).
- `:strict_lazy`: every value certified to 2⁻⁴⁴; +0–48% on cancellation-heavy single symbols.
- `:compensated_only`: always the compensated pass; ≤ 1e-16 everywhere, +10–18% on easy symbols.
- `:default`, `:strict`: the same thresholds without ratio reuse (kept for comparison).

A configuration knob, not per-call state: set it before computing, not concurrently.
"""
const POLICY = Ref(:lazy)

"""
Below this many ratio steps the lazy policy keeps no ratio buffer (the values are identical either way).
Allocating one costs ~9 ns and recomputing a step's ratio ~10–20 ns, so only a one-step sum skips it.
"""
const LAZY_MIN_STEPS = 2

"""
    level_escalate(s, segs, k, T, ztab) -> T

Finishes an `:escalate` case: identities or exact confirmation resolve zeros; nonzero sums go to multiword
tiers, then `BigFloat`, doubling until the estimate meets the target. Keep off worker threads where
`BigFloat` precision is shared.
"""
function level_escalate(s::FactorialSum, segs, k::Int, ::Type{T}, ztab::LevelZeroTable;
                        labels = nothing, workspace = nothing) where {T}
    # The recurrence comes first when available; it accepts entries using its own error estimate.
    if T === Float64 && labels !== nothing       # the symbol as one entry of its column: O(distance), no κ
        v = sixj_entry(labels, k, _column_workspace(workspace))
        v === nothing || return T(v)
    end
    # the pairwise test is a proof and costs a few comparisons; the modular screen is a pass over every term
    (pairwise_zero(s) || reflection_zero(s, segs, k)) && return zero(T)
    if is_cancellation_zero(s, segs, k, ztab) !== false
        iszero(exact_x(s, k)) && return zero(T)
    end
    target = _target_digits(T)
    if T === Float64                             # K-word tiers: κ up to ~1e30 (K = 3) and ~1e46 (K = 4)
        for tier in (Val(3), Val(4))
            v, B = _sum_at_level(s, segs, k, mw_tables(tier, k))
            _certifies(v, B, RTOL_CERTIFIED) && return T(Float64(v))
        end
    end
    # Reaching here, every fixed-precision pass has lost all of its digits, so none of them can say how much
    # cancellation there is: start from the minimum and let the doubling find it. Growing geometrically from
    # below costs at most twice the pass that finally succeeds; starting from a guess that is neither
    # informed nor minimal costs more, which is exactly what a bound-derived "loss" does here.
    pess = _bound_pessimism(s, segs, qint_tables(Float64, k))
    bits = 64 * cld(ceil(Int, (target + 26) * log2(10)), 64)
    while true
        w, B = setprecision(BigFloat, bits) do
            _sum_at_level(s, segs, k, BigFloat)
        end
        _surviving_digits(w, B, pess) >= target + 2 && return T(w)
        bits *= 2
    end
end
