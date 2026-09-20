# ---------------------------------------------------------------------------------
#  Summation kernels for a factorial rule, with certified error bounds
#
#  The plain ratio loop with a running error bound (Theorem 1 of `dev/results/certified_compensated.md`),
#  compensated Horner over double-word ratios (Theorem 2), the policy that decides which of them is kept,
#  and the escalation ladder above them: exact zero test, K-word tiers, BigFloat.
# ---------------------------------------------------------------------------------

"""
    _sum_at_level(s, segs, k, tab) -> (value, loss, bound, kappa)

Σ over the contributing segments in type T by the term ratio, with the prefactor.

- `loss`: decimal digits lost to cancellation, log₁₀(max |term| / |Σ|) (the old estimate).
- `bound`: a rigorous bound on |value − exact| (see `dev/results/certified_compensated.md`, Theorem 1).
  Every table entry is the correct rounding of its value, so a term reached after j ratio steps carries at
  most j(2K+1) relative roundings (K factorials per ratio); summation adds at most u Σ|partial sums|; the
  prefactor and first terms are products of correctly rounded split entries.
- `kappa`: Σ|terms| / |Σ|, the condition number of the sum that every fixed-precision method pays.

All bookkeeping is in split form, mantissa · 2^exponent: rescaling is by exact powers of two and segments
are combined by aligning exponents. No `log` or `exp` is taken.
"""
function _sum_at_level(s::FactorialSum, segs, k::Int, tab::QIntTables{T},
                       buf::B = nothing) where {T,B<:Union{Nothing,Vector}}
    q, qi = tab.q, tab.qi
    ql, qil = tab.ql, tab.qil
    ib = 0                                         # buffer cursor: (rh, rl) per step, segment by segment
    u = _unit(T)
    mp = one(T); ep = 0; np = 0
    @inbounds for (n, c) in s.pre
        mp, ep = _split_mul(mp, ep, tab, Int(n), c)
        np += abs(Int(c))
    end
    mp, ep = _renorm(mp, ep)
    relpre = _gamma(T, 2np)
    if s.sqrt_pre
        isodd(ep) && (mp *= 2; ep -= 1)
        mp = sqrt(mp); ep = ep ÷ 2
        relpre = relpre / 2 + u
    end
    K = 0
    for f in s.fac
        K += abs(Int(f.c))
    end
    intr = _integer_ratios(s, segs, tab)
    qq, qql = tab.qq, tab.qql
    unit, plan = _ratio_plan(s, length(q))
    rsign = s.alternating ? -one(T) : one(T)
    Kratio = _ratio_cost(s)
    cstep = intr ? 2 : 2Kratio + 1                      # roundings per step: a/b and t·r, or K entries + K products
    relfirst = _gamma(T, 2K)
    big_ = ldexp(one(T), SCALE_BITS)
    am = zero(T); ae = 0; started = false          # the sum, and its error bound and Σ|terms| in the same frame
    Eacc = zero(T); Aacc = zero(T)
    tm = zero(T); te = 0; seen = false             # the largest |term|
    @inbounds for seg in segs
        z0 = first(seg)
        nsteps = length(seg) - 1
        m0 = one(T); e0 = 0
        for f in s.fac
            m0, e0 = _split_mul(m0, e0, tab, _arg(f, z0), f.c)
        end
        m0, e0 = _renorm(m0, e0)
        t = one(T); ssum = one(T); tmax = one(T); sc = 0
        W = zero(T); SA = one(T); SS = zero(T); j = 0
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
            at > tmax && (tmax = at)
            W = fma(T(j), at, W); SA += at; SS += abs(ssum)
            if at > big_
                t = ldexp(t, -SCALE_BITS); ssum = ldexp(ssum, -SCALE_BITS); tmax = ldexp(tmax, -SCALE_BITS)
                W = ldexp(W, -SCALE_BITS); SA = ldexp(SA, -SCALE_BITS); SS = ldexp(SS, -SCALE_BITS)
                sc += SCALE_BITS
            end
        end
        fr, ex = frexp(tmax * m0)
        ex += sc + e0
        if !seen || ex > te || (ex == te && fr > tm)
            tm, te, seen = fr, ex, true
        end
        # error of this segment in units of 2^(sc+e0): terms, summation, first term, final product
        g = _gamma(T, cstep * nsteps)
        Eraw = abs(m0) * (cstep * u * W / (1 - g)^2 + u / (1 - u) * SS) + abs(ssum * m0) * (relfirst + u)
        Araw = abs(m0) * SA
        sgn = (s.alternating && isodd(z0)) ? -one(T) : one(T)
        vm, ve = frexp(sgn * ssum * m0)
        ve = iszero(vm) ? sc + e0 : ve + sc + e0
        Em = ldexp(Eraw, sc + e0 - ve); Am = ldexp(Araw, sc + e0 - ve)
        if !started
            am, ae, Eacc, Aacc, started = vm, ve, Em, Am, true
        else
            if ve > ae
                am = ldexp(am, ae - ve); Eacc = ldexp(Eacc, ae - ve); Aacc = ldexp(Aacc, ae - ve); ae = ve
                am += vm; Eacc += Em; Aacc += Am
            else
                am += ldexp(vm, ve - ae); Eacc += ldexp(Em, ve - ae); Aacc += ldexp(Am, ve - ae)
            end
            Eacc += u * abs(am)
        end
        if !iszero(am)
            fr2, ex2 = frexp(am)
            am = fr2; ae += ex2
            Eacc = ldexp(Eacc, -ex2); Aacc = ldexp(Aacc, -ex2)
        end
    end
    started || return zero(T), Inf, zero(T), T(Inf)
    bound = ldexp(Eacc * abs(mp) * (1 + relpre) + abs(am * mp) * (relpre + u), ae + ep) * (1 + T(1e-6))
    iszero(am) && return zero(T), Inf, bound, T(Inf)
    value = s.sign0 * ldexp(am * mp, ae + ep)
    loss = log10(Float64(tm) / abs(Float64(am))) + (te - ae) * log10(2.0)
    kappa = Aacc / abs(am)
    return value, loss, bound, kappa
end

"""
    _sum_compensated(s, segs, k, tab) -> (value, bound)

Compensated Horner evaluation of the same sum (after Graillat–Langlois–Louvet 2005), extended to ratios
that are themselves products of double-word table entries:

    Y_i = 1 + r_i Y_{i+1}  from the top of each segment,  S_segment = t_first · Y_0,

where every product r_i·y and sum 1 + p is split exactly by TwoProd/TwoSum, the low part of each ratio is
tracked to first order, and the rounding errors are propagated in a correction c_i = r_i c_{i+1} + … in
working precision. Prefactor, first terms, segment combination and the final product are double-word.
The result is as accurate as the plain loop run in doubled precision: relative error ≤ u + O((c n u)²) κ
(Theorem 2 of the write-up). `bound` is the a posteriori version of that estimate.
"""
function _sum_compensated(s::FactorialSum, segs, k::Int, tab::QIntTables{T},
                          buf::B = nothing) where {T,B<:Union{Nothing,Vector}}
    q, qi, ql, qil = tab.q, tab.qi, tab.ql, tab.qil
    base = 0                                       # buffer offset of the current segment
    u = _unit(T)
    ph = one(T); pl = zero(T); ep = 0
    @inbounds for (n, c) in s.pre
        ph, pl, ep = _split_mul_dw(ph, pl, ep, tab, Int(n), c)
    end
    ph, pl, ep = _renorm_dw(ph, pl, ep)
    if s.sqrt_pre
        isodd(ep) && (ph *= 2; pl *= 2; ep -= 1)
        ph, pl = _dw_sqrt(ph, pl); ep = ep ÷ 2
    end
    K = 0
    for f in s.fac
        K += abs(Int(f.c))
    end
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
        mh = one(T); ml = zero(T); e0 = 0
        for f in s.fac
            mh, ml, e0 = _split_mul_dw(mh, ml, e0, tab, _arg(f, z0), f.c)
        end
        mh, ml, e0 = _renorm_dw(mh, ml, e0)
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
Accuracy policy for `Float64` values. A plain-pass value is kept when its *certified* relative error bound is
at most `RTOL_PLAIN`; otherwise the compensated pass runs and is kept when its certified bound is at most
`RTOL_CERTIFIED`; only then does evaluation fall back to the exact zero test and higher precision. The plain
bound is pessimistic by ~10–1000× (`dev/results/certified_compensated.md`), so kept plain values are in
practice accurate to a few units in the last place.
"""
const RTOL_PLAIN = 2.0^-40          # ≈ 9.1e-13 certified
const RTOL_CERTIFIED = 2.0^-44

"Decimal digits a level value in type T must keep; below this the sum is redone at higher precision."
_target_digits(::Type{T}) where {T} = T === BigFloat ? _decimal_digits(T) - 4 : min(12.0, _decimal_digits(T) - 4)

"""
    value_at_level(s, k, T; fallback) -> T

Value at q = e^{iπ/(k+2)}. Poles throw a `DomainError` and vanishing terms are skipped, both decided from
valuations before any arithmetic. The sum runs in T and reports how many digits cancellation consumed.
If that leaves fewer than `_target_digits(T)`, the modular test first decides whether the value is exactly
zero; otherwise the sum is redone in BigFloat with enough bits for the loss (doubling until the digits
are there). The digit count is an estimate from max|term|/|Σ|, so a `Float64` result that is accepted
carries about 11 significant digits or better; escalated results carry the full target. `fallback()` is used if a factorial falls outside the level tables. With `labels` (the doubled labels of a
6j), the escalation first tries the column recurrence, which has no condition number at all.
"""
function value_at_level(s::FactorialSum, k::Int, ::Type{T}; fallback, labels = nothing, family = nothing, workspace = nothing) where {T}
    v, status, segs = level_pass1(s, k, qint_tables(T, k); family=family, workspace=workspace)
    status === :done && return v
    status === :fallback && return fallback()
    return level_escalate(s, segs, k, T, level_zero_table(k); labels=labels, workspace=workspace)
end

"""
    level_pass1(s, k, tab) -> (value, status, segments)

First pass in the table's own type, with no locking and no `BigFloat`, so it is safe to run on many
threads. `status` is `:done` when the value keeps enough digits, `:escalate` when cancellation consumed
them (the caller then runs the zero test and higher precision), or `:fallback` when a factorial falls
outside the level tables.
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

Plain pass with its certified bound; if the bound is too loose, the compensated pass with its own. Thread
safe (no global state). For types other than `Float64` the heuristic digit budget is used, as before.
"""
function _certified_value(s::FactorialSum, segs, k::Int, tab::QIntTables{T}, workspace=nothing) where {T}
    if T !== Float64
        v, loss, _, _ = _sum_at_level(s, segs, k, tab)
        return v, (_decimal_digits(T) - loss >= _target_digits(T) ? :done : :escalate)
    end
    mode = POLICY[]
    if mode === :compensated_only
        vc, Bc = _sum_compensated(s, segs, k, tab)
        return vc, ((isfinite(Bc) && !iszero(vc) && Bc <= RTOL_CERTIFIED * abs(vc)) ? :done : :escalate)
    end
    rtol = mode === :strict || mode === :strict_lazy ? RTOL_CERTIFIED : RTOL_PLAIN
    nsteps = sum(seg -> length(seg) - 1, segs; init = 0)
    # A one-step sum recomputes its ratio on fallback: cheaper than allocating a buffer on every call.
    lazy = (mode === :lazy || mode === :strict_lazy) && nsteps >= LAZY_MIN_STEPS
    buf = lazy ? _ratio_buffer(workspace, T, 2 * nsteps) : nothing
    v, _, B, _ = _sum_at_level(s, segs, k, tab, buf)
    (isfinite(B) && B <= rtol * abs(v)) && return v, :done
    vc, Bc = _sum_compensated(s, segs, k, tab, buf)
    (isfinite(Bc) && !iszero(vc) && Bc <= RTOL_CERTIFIED * abs(vc)) && return vc, :done
    return vc, :escalate
end

"""
Evaluation policy for `Float64` level and classical values (`dev/results/kfold_lazy_families.md`):

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

Finishes an `:escalate` case: the modular test first decides whether the value is exactly zero, otherwise
the sum is redone in `BigFloat` with enough bits for the loss, doubling until the digits are there. Changes
the global `BigFloat` precision while it runs, so callers keep it off worker threads.
"""
function level_escalate(s::FactorialSum, segs, k::Int, ::Type{T}, ztab::LevelZeroTable;
                        labels = nothing, workspace = nothing) where {T}
    # The recurrence comes first: when it returns a value its own estimate has certified that the entry is
    # far above the noise, so the entry cannot be an exact zero and the modular test (which costs a pass over
    # every term) is not needed.
    if T === Float64 && labels !== nothing       # the symbol as one entry of its column: O(distance), no κ
        v = sixj_entry(labels, k, _column_workspace(workspace))
        v === nothing || return T(v)
    end
    is_cancellation_zero(s, segs, k, ztab) === true && return zero(T)
    if T === Float64                             # K-word tiers: κ up to ~1e30 (K = 3) and ~1e46 (K = 4)
        for tier in (Val(3), Val(4))
            v, _, B, _ = _sum_at_level(s, segs, k, mw_tables(tier, k))
            (isfinite(B) && !iszero(v) && B <= RTOL_CERTIFIED * abs(v)) && return T(Float64(v))
        end
    end
    target = _target_digits(T)
    _, loss = _sum_at_level(s, segs, k, qint_tables(Float64, k))
    loss = max(loss, 16.0)                       # reaching here means Float64 compensation was not enough
    bits = 64 * cld(ceil(Int, (target + min(loss, 1e6) + 10) * log2(10)), 64)
    while true
        w, lossw = setprecision(BigFloat, bits) do
            _sum_at_level(s, segs, k, BigFloat)
        end
        bits * log10(2.0) - lossw >= target + 2 && return T(w)
        bits *= 2
    end
end
