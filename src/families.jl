# ---------------------------------------------------------------------------------
#  Families of level-k 6j symbols by a three-term recurrence
#
#  A column f(x) = {x j2 j3; l1 l2 l3}_q, x over its admissible range at level k, satisfies
#
#    up(x) f(x+1) + di(x) f(x) + up(x−1) f(x−1) = 0,        up(x) = E(x+1)/[2x+2],
#    di(x) = [2x+1] (B(x) + D(x) + C),
#    E(x)² = Dn(x) Bn(x−1),
#    Bn(x) = [x+1−j2+j3][x+2+j2+j3][x+1−l2+l3][l2+l3−x],        B(x) = Bn(x) / ([2x+1][2x+2]),
#    Dn(x) = [x+j2−j3][j2+j3+1−x][x+l2−l3][x+l2+l3+1],          D(x) = Dn(x) / ([2x][2x+1]),
#    C     = [j3−l1−l2][j3+l1−l2+1],
#
#  the q-analogue of the Schulten–Gordon recurrence (found numerically against the package's certified
#  values and verified to 1e-75; `dev/results/q_families.md`). The level truncation needs no special
#  case: E vanishes at the top of the level range because [h] = 0.
#
#  The Racah sum of a single symbol cancels (condition number κ, exponential in spin); the recurrence does
#  not, if each part of the column is computed in the direction in which the solution grows. So the column
#  is run forward from the left end to the first local maximum and backward from the right end to that
#  point, matched there, normalised by orthogonality Σ_x [2x+1][2l1+1] f(x)² = 1, and signed by the
#  certified value of one edge entry.
#
#  Arithmetic is double-word throughout. Rescaling f = h/P with P(x+1) = P(x) up(x) makes both directions
#  free of divisions and square roots, because up(x−1)² = E(x)²/[2x]² = Bn(x−1) Dn(x)/[2x]² is a product of
#  table entries:
#      forward   h(x+1) = −di(x) h(x) − up(x−1)² h(x−1),
#      backward  g(x−1) = −di(x) g(x) − up(x)² g(x+1)          (f = g/R, R(x−1) = R(x) up(x−1)),
#  and a square root is taken once per entry to leave the gauge (only the running values stay in it). Entrywise error ≈ n·u² relative to the
#  column's largest entry, so every entry that is not within ~10⁻¹⁵ of a node comes out to full precision.
# ---------------------------------------------------------------------------------

const DWord = Tuple{Float64,Float64}

@inline _dwm(a::DWord, b::DWord) = _dw_mul(a[1], a[2], b[1], b[2])
@inline function _dwa(a::DWord, b::DWord)
    s, e = _two_sum(a[1], b[1])
    return _fast_two_sum(s, e + (a[2] + b[2]))
end
@inline _dwn(a::DWord) = (-a[1], -a[2])
@inline _dws(a::DWord, t::Float64) = (a[1] * t, a[2] * t)           # exact scaling by a power of two
@inline function _dwdiv(a::DWord, b::DWord)
    q = a[1] / b[1]
    p = _dw_mul(q, 0.0, b[1], b[2])
    r = _dwa(a, _dwn(p))
    return _fast_two_sum(q, r[1] / b[1])
end
@inline _dwsqrt(a::DWord) = a[1] <= 0 ? (0.0, 0.0) : _dw_sqrt(a[1], a[2])

"""
Where the q-integers come from: `LevelQ` at q = e^{iπ/h} (level tables, any integer argument, since
[n+h] = −[n] and [−n] = −[n]), `ClassicalQ` at q = 1 ([n] = n, exact in a double word).
"""
struct LevelQ
    tab::QIntTables{Float64}
    k::Int
end
struct ClassicalQ end

@inline function _qint_dw(Q::LevelQ, n::Int)
    h = Q.k + 2
    r = mod(n, 2h)
    (r == 0 || r == h) && return (0.0, 0.0)
    @inbounds return r < h ? (Q.tab.q[r], Q.tab.ql[r]) : (-Q.tab.q[r-h], -Q.tab.ql[r-h])
end
@inline function _qinv_dw(Q::LevelQ, n::Int)                       # n not a multiple of h
    h = Q.k + 2
    r = mod(n, 2h)
    @inbounds return r < h ? (Q.tab.qi[r], Q.tab.qil[r]) : (-Q.tab.qi[r-h], -Q.tab.qil[r-h])
end
@inline _qint_dw(::ClassicalQ, n::Int) = (Float64(n), 0.0)          # exact for |n| < 2⁵³
@inline function _qinv_dw(::ClassicalQ, n::Int)
    x = Float64(n); r = 1 / x
    return (r, -fma(r, x, -1.0) / x)
end

"Admissible doubled x for {x j2 j3; l1 l2 l3} (doubled labels) at level k, or classically; empty if none."
function sixj_column_range(J2::Int, J3::Int, L1::Int, L2::Int, L3::Int, k::Int)
    lo = max(abs(J2 - J3), abs(L2 - L3))
    hi = min(J2 + J3, L2 + L3, 2k - J2 - J3, 2k - L2 - L3)
    ok = _qδ(L1, J2, L3, k) && _qδ(L1, L2, J3, k) && iseven(lo + J2 + J3) && iseven(lo + L2 + L3)
    return ok ? (lo:2:hi) : (lo:2:lo-2)
end
function sixj_column_range(J2::Int, J3::Int, L1::Int, L2::Int, L3::Int, ::Nothing)
    lo = max(abs(J2 - J3), abs(L2 - L3))
    hi = min(J2 + J3, L2 + L3)
    ok = _δ(L1, J2, L3) && _δ(L1, L2, J3) && iseven(lo + J2 + J3) && iseven(lo + L2 + L3)
    return ok ? (lo:2:hi) : (lo:2:lo-2)
end
_column_range(Q::LevelQ, J...) = sixj_column_range(J..., Q.k)
_column_range(::ClassicalQ, J...) = sixj_column_range(J..., nothing)

"Work arrays for one column (reused across columns by a worker)."
struct ColumnWork
    di::Vector{DWord}      # di(x)
    e::Vector{DWord}       # up(x−1)² at entry i (e[1] unused)
    w::Vector{DWord}       # the column, unnormalised
end
ColumnWork() = ColumnWork(DWord[], DWord[], DWord[])

function _resize!(w::ColumnWork, n::Int)
    for v in (w.di, w.e, w.w)
        length(v) < n && resize!(v, n)
    end
    return w
end


"""
    sixj_column!(out, J2, J3, L1, L2, L3, k, tab, work) -> range of doubled x
    sixj_column!(out, J2, J3, L1, L2, L3, Q, work)

Fills `out[i]` with {x j2 j3; l1 l2 l3} for the i-th admissible x (doubled labels, `Float64`), at level k
(`tab` the level's q-integer tables) or, with `Q = ClassicalQ()`, at q = 1. `work.w[i]` holds the same
values as double words.

With `xdim = true` each entry carries the extra factor √[2x+1] and with `cfac` a constant factor, which
together give the tetrahedrally symmetric normalisation (the G-symbol, and the F-symbol when its running
label is one of the two that carry a dimension). Both are folded into the square root that leaves the gauge
and into the normalising scale, so they cost nothing per entry — where multiplying each finished entry by
√(Π[2j+1]) costs six double-word products and a square root. The orthogonality sum is unchanged by the
folding, because Σ_x [2x+1] f² = Σ_x (√[2x+1] f)². Exact zeros inside the column come out as tiny values of the order of the
column's rounding error; callers that promise exact zeros test those entries.
"""
sixj_column!(out, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int, k::Int, tab::QIntTables{Float64},
             work::ColumnWork = ColumnWork(); kwargs...) =
    sixj_column!(out, J2, J3, L1, L2, L3, LevelQ(tab, k), work; kwargs...)

function sixj_column!(out, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int, Q::Union{LevelQ,ClassicalQ},
                      work::ColumnWork = ColumnWork(); xdim::Bool = false, cfac::DWord = (1.0, 0.0))
    X = _column_range(Q, J2, J3, L1, L2, L3)
    n = length(X)
    n == 0 && return X
    _resize!(work, n)
    di, e, w = work.di, work.e, work.w
    q(m) = _qint_dw(Q, m)
    C = _dwm(q((J3 - L1 - L2) ÷ 2), q((J3 + L1 - L2) ÷ 2 + 1))
    # ---- coefficients: Bn(x), Dn(x) products of four q-integers; E(x)² = Dn(x) Bn(x−1) ----
    Bprev = (0.0, 0.0)
    @inbounds for i in 1:n
        x2 = X[i]
        Bn = _dwm(_dwm(q((x2 - J2 + J3) ÷ 2 + 1), q((x2 + J2 + J3) ÷ 2 + 2)),
                  _dwm(q((x2 - L2 + L3) ÷ 2 + 1), q((L2 + L3 - x2) ÷ 2)))
        # [2x+1] B(x) = Bn / [2x+2]. At x = k/2 (possible only when j2+j3 = l2+l3 = k/2) [2x+2] = [h] = 0,
        # but Bn holds the factor [l2+l3−x] = [0] identically, so B = 0 there.
        t = iszero(Bn[1]) ? (0.0, 0.0) : _dwm(Bn, _qinv_dw(Q, x2 + 2))
        if x2 > 0
            Dn = _dwm(_dwm(q((x2 + J2 - J3) ÷ 2), q((J2 + J3 - x2) ÷ 2 + 1)),
                      _dwm(q((x2 + L2 - L3) ÷ 2), q((x2 + L2 + L3) ÷ 2 + 1)))
            iv = _qinv_dw(Q, x2)
            t = _dwa(t, _dwm(Dn, iv))                                  # [2x+1] D(x) = Dn / [2x]
            i > 1 && (e[i] = _dwm(_dwm(Dn, Bprev), _dwm(iv, iv)))      # up(x−1)² = E(x)² / [2x]²
        end
        di[i] = _dwa(t, _dwm(q(x2 + 1), C))
        Bprev = Bn
    end
    # ---- forward from the left end in the h-gauge, to the first local maximum of |f| ----
    # Only the running h(x−1), h(x), P(x)² are kept in the gauge; each entry leaves it as it is made, so a
    # rescaling (P² spans far more than the exponent range over a long column) touches three numbers.
    w[1] = xdim ? _dwsqrt(q(X[1] + 1)) : (1.0, 0.0)   # f(x₁) = 1 in the gauge, times the folded factor
    hm = (0.0, 0.0); hc = (1.0, 0.0); P2 = (1.0, 0.0)
    m = n
    @inbounds for i in 1:n-1
        hn = _dwn(_dwm(di[i], hc))
        i > 1 && (hn = _dwa(hn, _dwn(_dwm(e[i], hm))))
        P2 = _dwm(P2, e[i+1])                                          # P(x+1)² = P(x)² up(x)²
        w[i+1] = _leave_gauge(hn, P2, xdim ? q(X[i+1] + 1) : nothing)  # f = h / P, with √[2x+1] folded in
        hm, hc = hc, hn
        if !(2.0^-600 <= P2[1] <= 2.0^600)
            sc = P2[1] > 1 ? 2.0^-300 : 2.0^300
            hm = _dws(hm, sc); hc = _dws(hc, sc); P2 = _dws(P2, sc * sc)
        end
        if abs(w[i+1][1]) < abs(w[i][1])
            m = i
            break
        end
        if abs(w[i+1][1]) > 2.0^500                                    # the column spans a huge range
            _shrink!(w, 1:i+1)
            hm = _dws(hm, 2.0^-500); hc = _dws(hc, 2.0^-500)
        end
    end
    _divide!(w, 1:m, w[m])                                             # f(m) = 1 on the forward part
    # ---- backward from the right end in the g-gauge, to m, matched to the forward value there ----
    if m < n
        @inbounds begin
            w[n] = xdim ? _dwsqrt(q(X[n] + 1)) : (1.0, 0.0)
            gp = (0.0, 0.0); gc = (1.0, 0.0); R2 = (1.0, 0.0)           # g(x+1), g(x), R(x)²
            for i in n:-1:m+1
                gn = _dwn(_dwm(di[i], gc))
                i < n && (gn = _dwa(gn, _dwn(_dwm(e[i+1], gp))))         # up(x)² = e at x + 1
                R2 = _dwm(R2, e[i])                                      # R(x−1)² = R(x)² up(x−1)²
                fv = _leave_gauge(gn, R2, xdim ? q(X[i-1] + 1) : nothing)
                if i - 1 == m
                    _divide!(w, m+1:n, fv)                               # backward f(m) = 1 as well
                else
                    w[i-1] = fv
                    gp, gc = gc, gn
                    if !(2.0^-600 <= R2[1] <= 2.0^600)
                        s2 = R2[1] > 1 ? 2.0^-300 : 2.0^300
                        gp = _dws(gp, s2); gc = _dws(gc, s2); R2 = _dws(R2, s2 * s2)
                    end
                    if abs(fv[1]) > 2.0^500
                        _shrink!(w, i-1:n)
                        gp = _dws(gp, 2.0^-500); gc = _dws(gc, 2.0^-500)
                    end
                end
            end
        end
    end
    # ---- normalise by orthogonality, sign from a certified entry ----
    tot = (0.0, 0.0)
    @inbounds for i in 1:n
        sq = _dwm(w[i], w[i])
        tot = _dwa(tot, xdim ? sq : _dwm(q(X[i] + 1), sq))   # Σ [2x+1] f² either way
    end
    sc = _dwm(cfac, _dwdiv((1.0, 0.0), _dwsqrt(_dwm(q(L1 + 1), tot))))
    ie, ref = _sign_reference(w, X, J2, J3, L1, L2, L3, Q)
    signbit(ref) == signbit(w[ie][1]) || (sc = _dwn(sc))
    @inbounds for i in 1:n
        v = _dwm(w[i], sc)
        w[i] = v
        out === nothing || (out[i] = v[1] + v[2])
    end
    return X
end

"""
Leave the gauge: f = h/P, optionally with the factor √[2x+1] folded into the same square root. Passing the
`q`-integer rather than multiplying afterwards saves a double-word square root per entry.
"""
@inline _leave_gauge(h::DWord, P2::DWord, ::Nothing) = _dwdiv(h, _dwsqrt(P2))
@inline _leave_gauge(h::DWord, P2::DWord, xq::DWord) = _dwm(h, _dwsqrt(_dwdiv(xq, P2)))

"Scale entries by 2⁻⁵⁰⁰ (exact; entries far below the column's peak may underflow)."
@inline function _shrink!(w, rng)
    @inbounds for t in rng
        w[t] = _dws(w[t], 2.0^-500)
    end
    return nothing
end

"Divide entries by d (one double-word reciprocal, then products)."
@inline function _divide!(w, rng, d::DWord)
    r = _dwdiv((1.0, 0.0), d)
    @inbounds for t in rng
        w[t] = _dwm(w[t], r)
    end
    return nothing
end

"""
An entry of the column and the sign of its value, for the overall sign: the end with the shorter Racah
sum, else the other end, else the largest entry. A one-term sum whose factorials are all positive ([n] > 0
for 0 < n < h at a level; always classically) has the sign (−1)^z of its term, with no arithmetic.
"""
function _sign_reference(w, X, J2, J3, L1, L2, L3, Q)
    nt(x2) = ((α, β) = racah_sums(x2, J2, J3, L1, L2, L3); β[1] - α[4])
    ends = nt(first(X)) <= nt(last(X)) ? (1, length(X)) : (length(X), 1)
    for ie in ends
        sg = _entry_sign(Q, X[ie], J2, J3, L1, L2, L3)
        iszero(sg) || return ie, sg
    end
    ie = argmax(i -> abs(w[i][1]), eachindex(X))
    return ie, sign(_entry_value(Q, X[ie], J2, J3, L1, L2, L3))
end

function _entry_sign(Q, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    s = sixj_sum(X2, J2, J3, L1, L2, L3)
    is_empty_sum(s) && return 0.0
    if s.zlo == s.zhi && (Q isa ClassicalQ || _valuation_free(s, Q.k + 2))
        return (s.alternating && isodd(s.zlo)) ? -Float64(s.sign0) : Float64(s.sign0)
    end
    return sign(_entry_value(Q, X2, J2, J3, L1, L2, L3))
end

"The certified single-symbol value of one entry (doubled labels)."
function _entry_value(Q::LevelQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)::Float64
    k = Q.k
    s = sixj_sum(X2, J2, J3, L1, L2, L3)
    v, st, segs = level_pass1(s, k, Q.tab)
    st === :done && return v
    if st === :fallback
        js = canonical_spins(X2 // 2, J2 // 2, J3 // 2, L1 // 2, L2 // 2, L3 // 2)
        return real(project_discrete(q6j_dcr(js...), k, Float64))
    end
    return level_escalate(s, segs, k, Float64, level_zero_table(k))
end
_entry_value(::ClassicalQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int) =
    classical_value(sixj_sum(X2, J2, J3, L1, L2, L3), Float64)

# ---------------------------------------------------------------------------------
#  One hard symbol from the nearer end of its column
#
#  A symbol whose Racah sum has lost its digits (κ beyond the compensated range) is one entry of a column,
#  and the column recurrence has no κ at all. Computing the *whole* column would be O(n); we only need the
#  stretch between a seed and the target, so the cost is O(distance to the seed) with no dependence on κ.
#
#  Two seeds are needed for a three-term recurrence, and at a column end the boundary condition supplies the
#  second for free: E vanishes there, so the recurrence degenerates to f(edge+1) = −di(edge) f(edge)/up(edge).
#  One certified value at the end is therefore enough — and at the end the Racah sum usually has a single
#  term, so that value is exact and cheap (and is taken in split form, mantissa and binary exponent, because
#  at large spin an end value can be far below the smallest Float64 while the target is of order one).
#
#  The number of Racah terms, and with it κ, grows as x moves inward from the end. So the seed does not have
#  to be the end: `_seed_index` walks inward while the certified single-symbol kernel still succeeds, and the
#  recurrence then runs only the remaining distance (two certified seeds, no boundary condition needed).
# ---------------------------------------------------------------------------------

"""
Tetrahedral symmetries that bring position p of {j1 j2 j3; j4 j5 j6} to the first, so that any of the six
spins can play the role of the running label x of a column.
"""
const _TO_FIRST = ((1, 2, 3, 4, 5, 6), (2, 1, 3, 5, 4, 6), (3, 2, 1, 6, 5, 4),
                   (4, 5, 3, 1, 2, 6), (5, 4, 3, 2, 1, 6), (6, 5, 1, 3, 2, 4))

"How often the single-entry recurrence leaves the gauge to sample |f| for its error estimate."
const SAMPLE = 8

"""
Racah terms an entry may have and still be certified by the single-symbol kernel. Along a column the term
count grows with the distance from either end, and the certified region ends where the compensated pass runs
out (κ ≈ 10¹⁵); measured, that boundary sits at 60–82 terms over k = 1000–2000 and spins 100–250
(`dev/prototypes/check_turning_points.jl`), so a conservative threshold places the seed analytically instead
of probing for it.
"""
const SEED_TERMS = 48

"""
Largest rise above the target that the single-entry recurrence will accept, as a power of two. The error
estimate alone would tolerate a ratio of ~10¹⁶, which is too permissive for an entry that is an exact zero or
sits on a node: seeded from nearby, the run never sees the column's scale. Declining above 2²⁰ keeps those
entries with the modular test and the precision tiers, where they belong.
"""
const MAX_RISE_LOG2 = 20

"Steps saved must be worth the two certified seed evaluations."
const SEED_MIN_SAVING = 16

"Number of Racah terms of one entry of a column (doubled labels), in closed form."
@inline function _nterms(X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    α, β = racah_sums(X2, J2, J3, L1, L2, L3)
    return β[1] - α[4] + 1
end

"""
    _seed_index(X, i, from_left, J...) -> index or 0

The entry closest to the target from which the certified kernel can still start: the farthest index from the
chosen end whose Racah sum stays under `SEED_TERMS`. The term count is unimodal along a column, so a binary
search over the prefix (or suffix) finds it in a handful of O(1) evaluations. Returns 0 when the whole run
from the end is short anyway, in which case the end seed is used.
"""
function _seed_index(X, i::Int, from_left::Bool, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    n = length(X)
    nt(j) = _nterms(X[j], J2, J3, L1, L2, L3)
    if from_left
        lo, hi = 2, i - 1                       # need seeds at s−1 and s, and s < i
        hi < lo + 1 && return 0
        nt(hi) <= SEED_TERMS && return hi        # the certified region reaches the target's neighbour
        while lo < hi                            # largest s in [lo, hi] with nt(s) ≤ SEED_TERMS
            mid = (lo + hi + 1) ÷ 2
            nt(mid) <= SEED_TERMS ? (lo = mid) : (hi = mid - 1)
        end
        return nt(lo) <= SEED_TERMS ? lo : 0
    else
        lo, hi = i + 1, n - 1
        hi < lo && return 0
        nt(lo) <= SEED_TERMS && return lo
        while lo < hi
            mid = (lo + hi) ÷ 2
            nt(mid) <= SEED_TERMS ? (hi = mid) : (lo = mid + 1)
        end
        return nt(hi) <= SEED_TERMS ? hi : 0
    end
end

"A certified value of one entry, only if the kernel certifies it without escalating (`nothing` otherwise)."
function _certified_entry(Q::LevelQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    s = sixj_sum(X2, J2, J3, L1, L2, L3)
    is_empty_sum(s) && return nothing
    v, st, _ = level_pass1(s, Q.k, Q.tab)
    return st === :done ? v : nothing
end
function _certified_entry(Q::ClassicalQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    s = sixj_sum(X2, J2, J3, L1, L2, L3)
    is_empty_sum(s) && return nothing
    N = max_argument(s)
    v, st = _certified_value(s, (s.zlo:s.zhi,), 0, classical_tables(Float64, N))
    return st === :done ? v : nothing
end

"A value in split form: mantissa (double word) times 2^exp, so tiny end values do not underflow."
struct SplitValue
    m::DWord
    e::Int
end

_value(v::SplitValue) = ldexp(v.m[1] + v.m[2], v.e)

"""
    _edge_split(Q, X2, J2, J3, L1, L2, L3) -> SplitValue or nothing

The symbol at one end of a column, in split form. A one-term Racah sum is a product of split-table entries,
so it is exact to a few ulp at any size; otherwise the certified kernel is used, and `nothing` is returned
when its value has underflowed (then the caller cannot use this seed).
"""
function _edge_split(Q::LevelQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    s = sixj_sum(X2, J2, J3, L1, L2, L3)
    is_empty_sum(s) && return nothing
    tab = Q.tab
    if s.zlo == s.zhi && _valuation_free(s, Q.k + 2)          # single term: exact product of table entries
        mh = 1.0; ml = 0.0; ep = 0
        for (n, c) in s.pre
            mh, ml, ep = _split_mul_dw(mh, ml, ep, tab, Int(n), c)
        end
        mh, ml, ep = _renorm_dw(mh, ml, ep)
        if s.sqrt_pre
            isodd(ep) && (mh *= 2; ml *= 2; ep -= 1)
            mh, ml = _dw_sqrt(mh, ml); ep = ep ÷ 2
        end
        mt = 1.0; mtl = 0.0; et = 0
        for f in s.fac
            mt, mtl, et = _split_mul_dw(mt, mtl, et, tab, _arg(f, s.zlo), f.c)
        end
        m = _dwm((mh, ml), (mt, mtl))
        sg = s.sign0 * ((s.alternating && isodd(s.zlo)) ? -1 : 1)
        return SplitValue((sg * m[1], sg * m[2]), ep + et)
    end
    v = _entry_value(Q, X2, J2, J3, L1, L2, L3)
    return (iszero(v) || !isfinite(v)) ? nothing : SplitValue((v, 0.0), 0)
end

function _edge_split(Q::ClassicalQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    v = _entry_value(Q, X2, J2, J3, L1, L2, L3)
    return (iszero(v) || !isfinite(v)) ? nothing : SplitValue((v, 0.0), 0)
end

"Coefficients di and e of the recurrence for column indices `lo:hi` (1-based), into `work`."
function _coefficients!(work::ColumnWork, Q, X, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int,
                        lo::Int, hi::Int)
    di, e = work.di, work.e
    q(m) = _qint_dw(Q, m)
    C = _dwm(q((J3 - L1 - L2) ÷ 2), q((J3 + L1 - L2) ÷ 2 + 1))
    Bn(x2) = _dwm(_dwm(q((x2 - J2 + J3) ÷ 2 + 1), q((x2 + J2 + J3) ÷ 2 + 2)),
                  _dwm(q((x2 - L2 + L3) ÷ 2 + 1), q((L2 + L3 - x2) ÷ 2)))
    Dn(x2) = _dwm(_dwm(q((x2 + J2 - J3) ÷ 2), q((J2 + J3 - x2) ÷ 2 + 1)),
                  _dwm(q((x2 + L2 - L3) ÷ 2), q((x2 + L2 + L3) ÷ 2 + 1)))
    Bprev = lo > 1 ? Bn(X[lo-1]) : (0.0, 0.0)
    @inbounds for i in lo:min(hi + 1, length(X))
        x2 = X[i]
        B = Bn(x2)
        t = iszero(B[1]) ? (0.0, 0.0) : _dwm(B, _qinv_dw(Q, x2 + 2))
        if x2 > 0
            D = Dn(x2)
            iv = _qinv_dw(Q, x2)
            t = _dwa(t, _dwm(D, iv))
            i > 1 && (e[i] = _dwm(_dwm(D, Bprev), _dwm(iv, iv)))
        end
        di[i] = _dwa(t, _dwm(q(x2 + 1), C))
        Bprev = B
    end
    return nothing
end

"""
    sixj_entry(J, k, work; rtol) -> value or nothing

One level-k 6j symbol (doubled labels, x in the first position) by recursion along its column from the nearer
seed: cost O(distance), independent of the Racah condition number, and accurate to a few ulp unless the target
sits very close to a node of its column. Returns `nothing` when no usable seed exists or when the estimated
relative error exceeds `rtol`, so the caller can fall back to the precision tiers.
"""
sixj_entry(J::NTuple{6,Int}, k::Int, work::ColumnWork = ColumnWork(); rtol::Float64 = 1e-14) =
    sixj_entry(J, LevelQ(qint_tables(Float64, k), k), work; rtol = rtol)

"""
    _best_column(J, Q) -> (labels with x first, steps)

Any of the six spins can be the running label of a column, and the six columns differ in length and in how
far the symbol sits from an end. This picks the one with the fewest recurrence steps.
"""
function _best_column(J::NTuple{6,Int}, Q)
    best = J; bsteps = typemax(Int)
    for σ in _TO_FIRST
        c = (J[σ[1]], J[σ[2]], J[σ[3]], J[σ[4]], J[σ[5]], J[σ[6]])
        X = _column_range(Q, c[2], c[3], c[4], c[5], c[6])
        n = length(X)
        n <= 1 && continue
        d = c[1] - first(X)
        (iseven(d) && 0 <= d <= last(X) - first(X)) || continue
        i = d ÷ 2 + 1
        st = min(i - 1, n - i)
        if st < bsteps
            bsteps = st; best = c
        end
    end
    return best, bsteps
end

function sixj_entry(J0::NTuple{6,Int}, Q::Union{LevelQ,ClassicalQ}, work::ColumnWork = ColumnWork();
                    rtol::Float64 = 1e-14)
    J, _ = _best_column(J0, Q)                                # the shortest run among the six columns
    X2, J2, J3, L1, L2, L3 = J
    X = _column_range(Q, J2, J3, L1, L2, L3)
    n = length(X)
    n == 0 && return nothing
    d = X2 - first(X)
    (iseven(d) && 0 <= d <= last(X) - first(X)) || return nothing
    i = d ÷ 2 + 1
    n == 1 && return nothing                                  # the symbol *is* its own edge
    _resize!(work, n)
    from_left = i - 1 <= n - i
    # Where to start. Two certified entries near the end of the certified region start the recursion directly
    # and shorten the run; failing that, one value at the column end does, because the boundary condition
    # supplies the second. `s == 0` means "use the end".
    s = _seed_index(X, i, from_left, J2, J3, L1, L2, L3)
    steps_end = from_left ? i - 1 : n - i
    s == 0 || (steps_end - (from_left ? i - s : s - i) >= SEED_MIN_SAVING) || (s = 0)
    seedm = seedc = (0.0, 0.0)                                # f at the two seed entries
    if s != 0
        s2 = from_left ? s - 1 : s + 1
        v1 = _certified_entry(Q, X[s], J2, J3, L1, L2, L3)
        v2 = v1 === nothing ? nothing : _certified_entry(Q, X[s2], J2, J3, L1, L2, L3)
        if v1 === nothing || v2 === nothing || iszero(v1) || iszero(v2)
            s = 0                                             # not usable: fall back to the end
        else
            seedc = (v1, 0.0); seedm = (v2, 0.0)
        end
    end
    seed = s == 0 ? _edge_split(Q, from_left ? first(X) : last(X), J2, J3, L1, L2, L3) : nothing
    (s == 0 && seed === nothing) && return nothing
    lo = from_left ? (s == 0 ? 1 : s - 1) : i
    hi = from_left ? i : (s == 0 ? n : s + 1)
    _coefficients!(work, Q, X, J2, J3, L1, L2, L3, lo, hi)
    di, e = work.di, work.e
    # Gauge recurrence: forward h(x+1) = −di h(x) − up(x−1)² h(x−1), backward g(x−1) = −di g(x) − up(x)² g(x+1).
    # f = h·2^eH / √(P²·2^eP) (eP kept even), so rescaling is exact and nothing over- or underflows.
    hm = (0.0, 0.0); hc = (1.0, 0.0); P2 = (1.0, 0.0); eP = 0; eH = 0
    fv = (1.0, 0.0); fexp = 0; fmax_log = 0
    if s != 0
        # start at the seed pair: with P = 1 at the outer seed, h = f there and h = f·up at the inner one,
        # where up² is the coefficient the recursion already carries
        up2 = from_left ? e[s] : e[s+1]
        up2[1] > 0 || return nothing
        up = _dwsqrt(up2)
        hm = seedm; hc = _dwm(seedc, up); P2 = up2
        fv = seedc
        fmax_log = max(exponent(abs(seedm[1])), exponent(abs(seedc[1])))
    end
    steps = from_left ? ((s == 0 ? 1 : s):i-1) : ((s == 0 ? n : s):-1:i+1)
    isempty(steps) && return s == 0 ? nothing : ldexp(seedc[1] + seedc[2], 0)
    @inbounds for j in steps
        hn = _dwn(_dwm(di[j], hc))
        if from_left
            j > 1 && (hn = _dwa(hn, _dwn(_dwm(e[j], hm))))
            P2 = _dwm(P2, e[j+1])
        else
            j < n && (hn = _dwa(hn, _dwn(_dwm(e[j+1], hm))))
            P2 = _dwm(P2, e[j])
        end
        hm = hc; hc = hn
        if !(2.0^-600 <= P2[1] <= 2.0^600)
            grew = P2[1] > 1
            P2 = _dws(P2, grew ? 2.0^-600 : 2.0^600)
            eP += grew ? 600 : -600
        end
        if abs(hc[1]) > 2.0^500
            hm = _dws(hm, 2.0^-500); hc = _dws(hc, 2.0^-500); eH += 500
        end
        # Only the last entry is wanted, so the gauge is left every SAMPLE steps (and at the end), just to
        # track how far the path rises above the target — the factor in the error estimate.
        if j == last(steps) || rem(j, SAMPLE) == 0
            fv = _dwdiv(hc, _dwsqrt(P2))
            fexp = eH - eP ÷ 2
            iszero(fv[1]) || (fmax_log = max(fmax_log, exponent(abs(fv[1])) + fexp))
        end
    end
    iszero(fv[1]) && return nothing
    # How far the path rose above the target. A large rise means the target is near a node (or is an exact
    # zero), and then this tier must not answer at all, whatever its rounding estimate says.
    rise_log = fmax_log - (exponent(abs(fv[1])) + fexp)
    rise_log > MAX_RISE_LOG2 && return nothing
    # the double-word error of the run, relative to the target entry
    est = length(steps) * 2.0^-104 * exp2(rise_log)
    # a certified seed carries its own rounding, which the run then propagates
    s == 0 || (est += 4 * 2.0^-53)
    est > rtol && return nothing
    if s != 0
        return ldexp(fv[1] + fv[2], fexp)                     # the seeds are already normalised values
    end
    val = _dwm(fv, seed.m)
    iszero(val[1]) && return nothing
    return ldexp(val[1] + val[2], seed.e + fexp)
end
