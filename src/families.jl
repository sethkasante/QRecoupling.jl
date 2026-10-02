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
@inline _dwn(a::Tuple{T,T}) where {T<:Number} = (-a[1], -a[2])
@inline _dws(a::Tuple{T,T}, t::Float64) where {T<:Number} = (a[1] * t, a[2] * t)  # exact by a power of two
@inline function _dwdiv(a::DWord, b::DWord)
    q = a[1] / b[1]
    p = _dw_mul(q, 0.0, b[1], b[2])
    r = _dwa(a, _dwn(p))
    return _fast_two_sum(q, r[1] / b[1])
end
@inline _dwsqrt(a::DWord) = a[1] <= 0 ? (0.0, 0.0) : _dw_sqrt(a[1], a[2])

# ---------------------------------------------------------------------------------
#  The complex double word
#
#  A column at complex q needs the same arithmetic in ℂ. The representation is the obvious one — a high
#  and a low `ComplexF64`, so `a[1]` still means "the value" and `a[2]` "its correction", and every
#  magnitude test, rescaling and array in the kernel keeps working unchanged.
#
#  What it is *not* is a complex-valued Dekker word: `fma` has no complex counterpart, so there is no
#  error-free transformation to build on directly. Instead each of the real and imaginary parts is an
#  ordinary real double word, carried in the same pair; `_mk` reassembles them and `_re_dw`/`_im_dw` take
#  them apart. So a complex product is four real double-word products and two sums, and every component
#  inherits the real word's `u²` accuracy — the recurrence's error analysis carries over term by term
#  rather than having to be redone.
#
#  Division goes through the conjugate, `a/b = a·b̄/|b|²`, and the square root is the standard stable form
#  `√z = t + i y/2t` (or `|y|/2t + i·sign(y)·t` when Re z < 0) with `t = √((|z| + |Re z|)/2)`, which avoids
#  the cancellation the naive formula suffers on one side of the imaginary axis. Both form `|b|²` or
#  `|z|`, so they need their argument within the square root of the exponent range; the kernel's rescaling
#  keeps every running value inside 2^±600, which is exactly what that requires.
# ---------------------------------------------------------------------------------

const CDWord = Tuple{ComplexF64,ComplexF64}

@inline _re_dw(a::CDWord) = (real(a[1]), real(a[2]))
@inline _im_dw(a::CDWord) = (imag(a[1]), imag(a[2]))
@inline _mk(re::DWord, im::DWord) = (complex(re[1], im[1]), complex(re[2], im[2]))

@inline function _dwm(a::CDWord, b::CDWord)
    ar = _re_dw(a); ai = _im_dw(a); br = _re_dw(b); bi = _im_dw(b)
    return _mk(_dwa(_dwm(ar, br), _dwn(_dwm(ai, bi))),
               _dwa(_dwm(ar, bi), _dwm(ai, br)))
end
@inline _dwa(a::CDWord, b::CDWord) =
    _mk(_dwa(_re_dw(a), _re_dw(b)), _dwa(_im_dw(a), _im_dw(b)))

@inline function _dwdiv(a::CDWord, b::CDWord)
    ar = _re_dw(a); ai = _im_dw(a); br = _re_dw(b); bi = _im_dw(b)
    d = _dwa(_dwm(br, br), _dwm(bi, bi))
    return _mk(_dwdiv(_dwa(_dwm(ar, br), _dwm(ai, bi)), d),
               _dwdiv(_dwa(_dwm(ai, br), _dwn(_dwm(ar, bi))), d))
end

@inline _dwabs(a::DWord) = a[1] < 0 ? _dwn(a) : a

"Drop an imaginary part that is rounding noise, so that `_dwsqrt` lands on a definite branch."
@inline _re_only(a::CDWord) = (ComplexF64(real(a[1]), 0.0), ComplexF64(real(a[2]), 0.0))

"Principal square root of a complex double word."
@inline function _dwsqrt(a::CDWord)
    x = _re_dw(a); y = _im_dw(a)
    (iszero(x[1]) && iszero(y[1])) && return (ComplexF64(0), ComplexF64(0))
    if iszero(y[1]) && iszero(y[2])                       # exactly real: keep the real branch exactly
        return x[1] > 0 ? _mk(_dwsqrt(x), (0.0, 0.0)) : _mk((0.0, 0.0), _dwsqrt(_dwn(x)))
    end
    m = _dwsqrt(_dwa(_dwm(x, x), _dwm(y, y)))             # |z|
    t = _dwsqrt(_dws(_dwa(m, _dwabs(x)), 0.5))            # the larger of the two parts
    if x[1] >= 0
        return _mk(t, _dwdiv(y, _dws(t, 2.0)))
    end
    return _mk(_dwdiv(_dwabs(y), _dws(t, 2.0)), y[1] < 0 ? _dwn(t) : t)
end

"The word type a column uses, and its neutral elements."
_wordtype(::Any) = DWord
@inline _dwzero(::Type{DWord}) = (0.0, 0.0)
@inline _dwone(::Type{DWord}) = (1.0, 0.0)
@inline _dwzero(::Type{CDWord}) = (ComplexF64(0), ComplexF64(0))
@inline _dwone(::Type{CDWord}) = (ComplexF64(1), ComplexF64(0))

"""
Magnitude of a high part, for the kernel's rescaling tests. `Float64` is returned unchanged on purpose:
the real path's comparisons, including their behaviour on a negative `P²`, are left exactly as measured.
"""
@inline _wmag(x::Float64) = x
@inline _wmag(z::ComplexF64) = abs(z)

"`ldexp` on whichever scalar a word carries."
@inline _wldexp(x::Float64, e::Int) = ldexp(x, e)
@inline _wldexp(z::ComplexF64, e::Int) = complex(ldexp(real(z), e), ldexp(imag(z), e))

"""
Whether `up(x−1)² = e` admits a usable `up`. On the real axis the recurrence needs it *positive*: a
negative one is a forbidden interval the gauge cannot enter. Off the axis there is no sign to test — the
square root is fixed by the package's branch either way — so only an exact zero disqualifies it.
"""
@inline _w_has_root(a::DWord) = a[1] > 0
@inline _w_has_root(a::CDWord) = !iszero(a[1])

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

"""
    RealQ(q, N)

A column at generic **real** `q`, built from double-word tables of `[n]` and `1/[n]` for `n = 0…N`. The
recurrence itself is q-generic — nothing in the gauge form of §`sixj_column!` mentions a level — so the only
thing that tied families to levels and to `q = 1` was the absence of a table. With this, a column at real q
costs `O(n)` instead of n independent sums, the same as at a level.

Entries are computed once in 128-bit arithmetic and split into double words, so each is the correct rounding
of its value with a low part good to about `u²`, matching what `LevelQ` reads from `QIntTables`.

Off the real axis the same job is done by [`ComplexQ`](@ref), which carries the complex double word.
"""
struct RealQ
    q::Float64
    ints::Vector{DWord}      # [n], n = 0…N
    invs::Vector{DWord}      # 1/[n]
    cancel::Float64          # worst cancellation seen building the table; see `column_accuracy`
end

function RealQ(q::Real, N::Int)
    q > 0 || throw(DomainError(q, "RealQ needs a positive real q"))
    isapprox(q, 1) && throw(ArgumentError("q ≈ 1: use ClassicalQ for the classical column"))
    N >= 0 || throw(DomainError(N, "table length must be nonnegative"))
    ints = Vector{DWord}(undef, N + 1)
    invs = Vector{DWord}(undef, N + 1)
    worst = 1.0
    setprecision(BigFloat, 128) do
        qb = BigFloat(q)
        x = qb + inv(qb)                        # [n+1] = x[n] − [n−1]: one multiply and one subtract per
        prev = zero(BigFloat)                   # entry, against two `^` — and no repeated rounding of a
        cur = one(BigFloat)                     # power. |[n]| grows either side of q = 1, so the forward
        for n in 0:N                            # recurrence follows the dominant solution and is stable.
            v = iszero(n) ? zero(BigFloat) : cur
            hi = Float64(v)
            ints[n+1] = (hi, Float64(v - hi))
            if iszero(v)
                invs[n+1] = (Inf, 0.0)          # [0]: a valid rule never divides by it
            else
                w = inv(v); h = Float64(w)
                invs[n+1] = (h, Float64(w - h))
            end
            if n > 0                        # (a tuple assignment cannot hide inside `||`)
                t = x * cur
                nxt = t - prev
                worst = iszero(nxt) ? Inf : max(worst, Float64((abs(t) + abs(prev)) / abs(nxt)))
                prev, cur = cur, nxt
            end
        end
    end
    return RealQ(Float64(q), ints, invs, worst)
end

"""
Tables for a column at real `q`, cached on `(q, N)`. A column, an F-matrix or a batch at one `q` builds the
table once; without the cache the BigFloat build dominates and costs more than the loop it replaces
(measured 26–297 µs against a 0.5–2.9 µs fill).
"""
const REALQ_CACHE = LRU{Tuple{Float64,Int},Any}(maxsize = 64)
const REALQ_LOCK = ReentrantLock()

function real_q_tables(q::Real, N::Int)
    key = (Float64(q), N)
    lock(REALQ_LOCK) do
        get!(REALQ_CACHE, key) do
            RealQ(q, N)
        end
    end::RealQ
end

real_q_tables(q::Real, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int) =
    real_q_tables(q, _realq_table_size(sixj_column_range(J2, J3, L1, L2, L3, nothing),
                                       J2, J3, L1, L2, L3))

"Table length a column over `X` with these labels can reach."
_realq_table_size(X, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int) =
    max(last(X) + 2, (last(X) + J2 + J3) ÷ 2 + 2, (last(X) + L2 + L3) ÷ 2 + 2,
        (J2 + J3 + L1 + L2 + L3) ÷ 2 + 2, L1 + 2)

RealQ(q::Real, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int) =
    RealQ(q, _realq_table_size(sixj_column_range(J2, J3, L1, L2, L3, nothing), J2, J3, L1, L2, L3))

"""
    ComplexQ(q, N)

A column at generic **complex** `q`, from complex double-word tables of `[n]` and `1/[n]`. The recurrence
is q-generic — nothing in the gauge form mentions a level, a sign, or the reality of `q` — so the only
thing that kept families real was the word type. With [`CDWord`](@ref) in place a column off the real axis
costs `O(n)` like any other.

Entries are built in 256-bit `Complex{BigFloat}` from `[n+1] = x[n] − [n−1]`, `x = q + q⁻¹`, and split into
complex double words. The extra precision over [`RealQ`](@ref)'s 128 bits is not decoration: on `|q| = 1`
the q-integers neither grow nor decay, so that step is a genuine subtraction of comparable numbers and
loses digits near a root of unity. The build tracks exactly that loss — `(|x·[n]| + |[n−1]|)/|[n+1]|` — and
**refuses** the table when it exceeds `CANCEL_LIMIT`, rather than returning a column whose accuracy it
cannot stand behind. A root of unity is a level; use a level target there.
"""
struct ComplexQ
    q::ComplexF64
    ints::Vector{CDWord}      # [n], n = 0…N
    invs::Vector{CDWord}      # 1/[n]
    psi::Vector{CDWord}       # Ψ_d, d = 1…N, with [n] = ∏_{d | n} Ψ_d
    spsi::Vector{CDWord}      # √Ψ_d, one principal root each — the package's branch convention
    divs::Vector{Vector{Int}} # divisors of n, precomputed
    cancel::Float64           # worst cancellation seen building the table; see `CANCEL_LIMIT`
end

"""
    column_accuracy(Q::ComplexQ) -> Float64

Relative accuracy a column built from this table can promise, `u·cancel`.

`cancel` is the worst ratio `(|x[n]| + |[n−1]|)/|[n+1]|` met building the q-integers, so it is how many
digits the recurrence coefficients lose before the column starts. Away from a root of unity it is a few
units and the promise is the usual `1e−15`; at `q = e^{0.1111iπ}`, which is `1.1e−5` from a ninth root of
unity, it is large enough that a nearly-imaginary entry keeps only about eleven digits of its real part.
That is a property of the point, not a defect of the recurrence, and it is reported rather than hidden —
the alternative is a promise the table cannot keep.
"""
column_accuracy(Q::Union{RealQ,ComplexQ}) = max(eps(Float64), eps(Float64) * Q.cancel)

"""
Accuracy floor a near-edge run inherits from its table, on top of the double-word error of the run itself.

The table entries are built at 256 bits and rounded to double words, so they are accurate wherever `q` is.
What a run cannot do better than is the conditioning of everything assembled from them — the recurrence
coefficients, the Ψ-split roots, the seeds — and near a root of unity that is not small: at
`q = e^{0.1111iπ}` a j ≈ 38 entry came out `1.6e−12` from a 512-bit reference while the run's own
double-word estimate said `1e−28`, because the estimate charges the arithmetic and not the inputs.

`u·cancel` is the measured scale of it: over a sweep of 14 `q` and 60 label sets the near-edge error tracks
this within a factor of 2, and it is the same number `column_accuracy` reports; the 2 is that factor.
Tuned by measurement in both directions: at 1 the tier still answered two `e^{0.61iπ}` columns at 1.5–3.0e−14
against a promised 1e−14, and at 4 it began declining well-conditioned real-`q` runs that were accurate to
3.5e−17. Like the rest of the
estimate this is a measured heuristic, not a derived bound — `sixj_entry` has always said so — but with it
in place every answer the tier gives is inside the `rtol` it promises, and the cases it cannot are declined
to the arbitrary-precision ladder instead. A level and the classical limit have exact tables and add
nothing.
"""
@inline _table_floor(Q::Union{RealQ,ComplexQ}) = 2 * eps(Float64) * Q.cancel
@inline _table_floor(::Any) = 0.0

"""
Cancellation a `ComplexQ` table will tolerate in `[n+1] = x[n] − [n−1]`.

The 256-bit build is not what this protects: it could absorb far more. What cannot is the *column*. A
cancellation of `C` means some `[n]` is `C` times smaller than its neighbours, so `1/[n]` is `C` times
larger, and the recurrence coefficients built from it lose `log₂C` bits to the subtraction in `di`. A
complex double word carries 106, so 2⁴⁰ still leaves 66 — comfortably more than `Float64`.

The threshold is where it is because of what happens past it, measured on `fmatrix(1,1,1,1)`: at
`q = e^{iπ/4}` orthogonality came out at 1.0 and at `q = e^{iπ/3}` at 4.5e15, both silently. Those q are
roots of unity, where the table has a `[n]` that is only rounding noise; the cancellation ratio there is
about 2⁵², against about 2⁶ at a generic point on the unit circle. A root of unity is a level and has
exact tables of its own.
"""
const CANCEL_LIMIT = 2.0^40

function ComplexQ(q::Number, N::Int)
    qq = ComplexF64(q)
    isfinite(qq) && !iszero(qq) || throw(DomainError(q, "ComplexQ needs a finite nonzero q"))
    N >= 0 || throw(DomainError(N, "table length must be nonnegative"))
    ints = Vector{CDWord}(undef, N + 1)
    invs = Vector{CDWord}(undef, N + 1)
    worst = 1.0
    setprecision(BigFloat, 256) do
        qb = Complex{BigFloat}(qq)
        x = qb + inv(qb)
        prev = zero(qb); cur = one(qb)
        for n in 0:N
            v = iszero(n) ? zero(qb) : cur
            hi = ComplexF64(v)
            ints[n+1] = (hi, ComplexF64(v - hi))
            if iszero(v)
                invs[n+1] = (ComplexF64(Inf), ComplexF64(0))
            else
                w = inv(v); h = ComplexF64(w)
                invs[n+1] = (h, ComplexF64(w - h))
            end
            if n > 0
                t = x * cur
                nxt = t - prev
                worst = iszero(nxt) ? Inf : max(worst, Float64((abs(t) + abs(prev)) / abs(nxt)))
                prev, cur = cur, nxt
            end
        end
    end
    # The Ψ basis, and its roots. `[n] = ∏_{d | n} Ψ_d` — the same balanced factorisation
    # `analytic_rules.jl` uses — so dividing out the proper divisors in increasing order leaves Ψ.
    psi = Vector{CDWord}(undef, max(N, 1))
    @inbounds for n in 1:N
        psi[n] = ints[n+1]
    end
    @inbounds for d in 2:N÷2, n in 2d:d:N
        psi[n] = _dwdiv(psi[n], psi[d])
    end
    # On the unit circle every Ψ_d is exactly real, and a negative one sits on the branch cut. Which of
    # ±i√|Ψ_d| the root returns would then be decided by the rounding error in |q| − 1 rather than by q,
    # so the imaginary part is discarded and the branch fixed at +i√|Ψ_d| — the same rule, and the same
    # predicate, that `analytic_rules.jl` applies to the scalar path. A column and `q6j` must agree on it:
    # they disagreed by an exact sign on 8 of the suite's unit-circle entries before this was added.
    circle = _on_unit_circle(qq) || _negative_real_axis(qq)
    spsi = Vector{CDWord}(undef, max(N, 1))
    @inbounds for d in 1:N
        spsi[d] = _dwsqrt(circle ? _re_only(psi[d]) : psi[d])
    end
    divs = [Int[] for _ in 0:N]
    @inbounds for d in 1:N, n in d:d:N
        push!(divs[n+1], d)
    end
    if worst > CANCEL_LIMIT
        near = isapprox(abs(qq), 1; atol = 1e-8)
        hint = near ? "a level target — `k = $(round(Int, pi / max(abs(angle(qq)), eps())) - 2)` if that " *
                      "is the level you mean — has exact tables for it" :
                      "move q away from it, or evaluate symbol by symbol with `q6j`/`fsymbol`"
        throw(ArgumentError(
            "q = $qq is a root of unity to within the precision of a complex double word: the q-integers " *
            "cancel by a factor $(worst) building the table, so a column at this q cannot be certified. " *
            hint))
    end
    return ComplexQ(qq, ints, invs, psi, spsi, divs, worst)
end

"`√[n]` at complex q over the Ψ basis. Every Ψ_d divides `[n]` exactly once, so all of them are rooted."
function _sqrt_qint(Q::ComplexQ, n::Int)
    n == 0 && return _dwzero(CDWord)
    n < 0 && throw(DomainError(n, "the dimension factor √[n] needs n ≥ 0"))
    acc = _dwone(CDWord)
    @inbounds for d in Q.divs[n+1]
        acc = _dwm(acc, Q.spsi[d])
    end
    return acc
end

"""
    _sqrt_qint_product(Q, args, cnt) -> CDWord

`√(∏ᵢ [nᵢ])` at complex `q`, in the branch convention the rest of the package uses: rewrite the product
over the Ψ basis, `[n] = ∏_{d | n} Ψ_d`, take the even part of each exponent out of the root exactly, and
root what is left **one Ψ_d at a time**.

Why not simply `sqrt` of the product, or a root per q-integer: neither is the same function. `√` is not
multiplicative across its cut, so each factorisation of the same number gives a different branch, and only
one of them agrees with `analytic_rules.jl`'s prefactor — which splits over exactly this Ψ basis. Measured,
rooting per q-integer instead disagreed with `q6j` by a sign on 9 of 12 sampled complex `q`, including
points far from any cut; over the Ψ basis the disagreement is gone.
"""
function _sqrt_qint_product(Q::ComplexQ, args, cnt::Vector{Int})
    acc = _dwone(CDWord)
    @inbounds for n in args
        n == 0 && return _dwzero(CDWord)
        for d in Q.divs[n+1]
            cnt[d] += 1
        end
    end
    @inbounds for n in args, d in Q.divs[n+1]
        c = cnt[d]
        c == 0 && continue                      # already consumed
        cnt[d] = 0
        for _ in 1:(c ÷ 2)
            acc = _dwm(acc, Q.psi[d])
        end
        isodd(c) && (acc = _dwm(acc, Q.spsi[d]))
    end
    return acc
end

"""
Tables for a column at complex `q`, cached on `(q, N)` exactly as [`real_q_tables`](@ref) is: the 256-bit
build costs far more than the fill it replaces, so a matrix or a batch at one `q` must pay it once.
"""
const COMPLEXQ_CACHE = LRU{Tuple{ComplexF64,Int},Any}(maxsize = 64)
const COMPLEXQ_LOCK = ReentrantLock()

function complex_q_tables(q::Number, N::Int)
    key = (ComplexF64(q), N)
    lock(COMPLEXQ_LOCK) do
        get!(COMPLEXQ_CACHE, key) do
            ComplexQ(q, N)
        end
    end::ComplexQ
end

complex_q_tables(q::Number, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int) =
    complex_q_tables(q, _realq_table_size(sixj_column_range(J2, J3, L1, L2, L3, nothing),
                                          J2, J3, L1, L2, L3))

_wordtype(::ComplexQ) = CDWord

@inline function _qint_dw(Q::ComplexQ, n::Int)
    n < 0 && return _dwn(Q.ints[-n+1])
    return Q.ints[n+1]
end
@inline function _qinv_dw(Q::ComplexQ, n::Int)
    n < 0 && return _dwn(Q.invs[-n+1])
    return Q.invs[n+1]
end
_column_range(::ComplexQ, J...) = sixj_column_range(J..., nothing)

@inline function _qint_dw(Q::RealQ, n::Int)
    n < 0 && return _dwn(Q.ints[-n+1])
    return Q.ints[n+1]
end
@inline function _qinv_dw(Q::RealQ, n::Int)
    n < 0 && return _dwn(Q.invs[-n+1])
    return Q.invs[n+1]
end
_column_range(::RealQ, J...) = sixj_column_range(J..., nothing)

_column_range(Q::LevelQ, J...) = sixj_column_range(J..., Q.k)
_column_range(::ClassicalQ, J...) = sixj_column_range(J..., nothing)

"Work arrays for one column (reused across columns by a worker)."
struct ColumnWork{W}
    di::Vector{W}          # di(x)
    e::Vector{W}           # up(x−1)² at entry i (e[1] unused)
    w::Vector{W}           # the column, unnormalised
    up::Vector{W}          # up(x−1) itself, rooted factor by factor — complex q only
end
ColumnWork(::Type{W} = DWord) where {W} = ColumnWork{W}(W[], W[], W[], W[])
"Scratch of the word type this `Q` needs, reusing `work` when it already has it."
_column_work(Q, work::ColumnWork{W}) where {W} = W === _wordtype(Q) ? work : ColumnWork(_wordtype(Q))

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

function sixj_column!(out, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int,
                      Q::Union{LevelQ,ClassicalQ,RealQ,ComplexQ},
                      work::ColumnWork = ColumnWork(); xdim::Bool = false,
                      cfac = _dwone(_wordtype(Q)))
    X = _column_range(Q, J2, J3, L1, L2, L3)
    n = length(X)
    n == 0 && return X
    work = _column_work(Q, work)
    _resize!(work, n)
    _column_coefficients!(work, Q, X, J2, J3, L1, L2, L3)
    return _sixj_column_core!(out, J2, J3, L1, L2, L3, Q, work, X; xdim = xdim, cfac = cfac)
end

"""
    _column_shared!(shared, work, Q, X, J2, J3, L2, L3)

The part of the recurrence coefficients that does not involve `L1`: `work.e[i] = up(x−1)²` and
`shared[i] = [2x+1](B(x) + D(x))`. Everything here is built from four-factor products of q-integers in
(J2, J3, L2, L3), so for a *family* that varies only `L1` — an F-matrix column index, say — this is computed
once for the whole family and only the scalar `C` changes afterwards (see [`_column_di!`](@ref)).
"""
function _column_shared!(shared, work::ColumnWork, Q, X, J2::Int, J3::Int, L2::Int, L3::Int)
    e = work.e
    q(m) = _qint_dw(Q, m)
    W = _wordtype(Q)
    Z = _dwzero(W)
    cplx = W === CDWord
    cplx && length(work.up) < length(X) && resize!(work.up, length(X))
    cnt = cplx ? zeros(Int, length(Q.psi) + 1) : Int[]                 # Ψ exponents, reused per entry
    Bprev = Z
    bargs = (0, 0, 0, 0)                                               # arguments of Bn(x−1)'s factors
    @inbounds for i in eachindex(X)
        x2 = X[i]
        ba = ((x2 - J2 + J3) ÷ 2 + 1, (x2 + J2 + J3) ÷ 2 + 2,
              (x2 - L2 + L3) ÷ 2 + 1, (L2 + L3 - x2) ÷ 2)
        Bn = _dwm(_dwm(q(ba[1]), q(ba[2])), _dwm(q(ba[3]), q(ba[4])))
        # [2x+1] B(x) = Bn / [2x+2]. At x = k/2 (possible only when j2+j3 = l2+l3 = k/2) [2x+2] = [h] = 0,
        # but Bn holds the factor [l2+l3−x] = [0] identically, so B = 0 there.
        t = iszero(Bn[1]) ? Z : _dwm(Bn, _qinv_dw(Q, x2 + 2))
        if x2 > 0
            dargs = ((x2 + J2 - J3) ÷ 2, (J2 + J3 - x2) ÷ 2 + 1,
                     (x2 + L2 - L3) ÷ 2, (x2 + L2 + L3) ÷ 2 + 1)
            Dn = _dwm(_dwm(q(dargs[1]), q(dargs[2])), _dwm(q(dargs[3]), q(dargs[4])))
            iv = _qinv_dw(Q, x2)
            t = _dwa(t, _dwm(Dn, iv))                                  # [2x+1] D(x) = Dn / [2x]
            if i > 1
                e[i] = _dwm(_dwm(Dn, Bprev), _dwm(iv, iv))             # up(x−1)² = E(x)² / [2x]²
                # …and, at complex q, up(x−1) = E(x)/[2x] itself, rooted over the Ψ basis so that its
                # branch is the one the rest of the package takes. `_gauge_step` says why it is needed
                # and `_sqrt_qint_product` why the Ψ basis is the right place to split.
                cplx && (work.up[i] = _dwm(_sqrt_qint_product(Q,
                    (dargs[1], dargs[2], dargs[3], dargs[4], bargs[1], bargs[2], bargs[3], bargs[4]),
                    cnt), iv))
            end
        end
        shared[i] = t
        Bprev = Bn
        bargs = ba
    end
    return shared
end

"`di` from the shared part: only `C = [(J3−L1−L2)/2][(J3+L1−L2)/2+1]` depends on `L1`."
function _column_di!(work::ColumnWork, Q, X, J3::Int, L1::Int, L2::Int, shared)
    q(m) = _qint_dw(Q, m)
    C = _dwm(q((J3 - L1 - L2) ÷ 2), q((J3 + L1 - L2) ÷ 2 + 1))
    @inbounds for i in eachindex(X)
        work.di[i] = _dwa(shared[i], _dwm(q(X[i] + 1), C))
    end
    return work
end

"Both halves at once, for a single column."
function _column_coefficients!(work::ColumnWork, Q, X, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    shared = Vector{_wordtype(Q)}(undef, length(X))
    _column_shared!(shared, work, Q, X, J2, J3, L2, L3)
    return _column_di!(work, Q, X, J3, L1, L2, shared)
end

"""
The column recurrence itself, given coefficients already in `work`: forward to the first local maximum,
backward from the right end, matched there, normalised by orthogonality and signed by one certified entry.
Split out so that a family sharing its coefficients can call it directly.
"""
function _sixj_column_core!(out, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int, Q,
                            work::ColumnWork, X; xdim::Bool = false,
                            cfac = _dwone(_wordtype(Q)))
    n = length(X)
    n == 0 && return X
    di, e, w = work.di, work.e, work.w
    q(m) = _qint_dw(Q, m)
    W = _wordtype(Q); Z = _dwzero(W); U = _dwone(W)
    # ---- forward from the left end in the h-gauge, to the first local maximum of |f| ----
    # Only the running h(x−1), h(x), P(x)² are kept in the gauge; each entry leaves it as it is made, so a
    # rescaling (P² spans far more than the exponent range over a long column) touches three numbers.
    w[1] = xdim ? _xdim_root(Q, X[1] + 1) : U  # f(x₁) = 1 in the gauge, times the folded factor
    hm = Z; hc = U; P2 = U
    m = n
    @inbounds for i in 1:n-1
        hn = _dwn(_dwm(di[i], hc))
        i > 1 && (hn = _dwa(hn, _dwn(_dwm(e[i], hm))))
        P2 = _gauge_step(P2, work, i + 1)                              # P(x+1)² = P(x)² up(x)²
        w[i+1] = _leave_gauge(hn, P2, xdim ? _xdim_pass(Q, X[i+1] + 1) : nothing)   # √[2x+1] folded in
        hm, hc = hc, hn
        if !(2.0^-600 <= _wmag(P2[1]) <= 2.0^600)
            sc = _wmag(P2[1]) > 1 ? 2.0^-300 : 2.0^300
            hm = _dws(hm, sc); hc = _dws(hc, sc); P2 = _gauge_rescale(P2, sc)
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
            w[n] = xdim ? _xdim_root(Q, X[n] + 1) : U
            gp = Z; gc = U; R2 = U                                       # g(x+1), g(x), R(x)²
            for i in n:-1:m+1
                gn = _dwn(_dwm(di[i], gc))
                i < n && (gn = _dwa(gn, _dwn(_dwm(e[i+1], gp))))         # up(x)² = e at x + 1
                R2 = _gauge_step(R2, work, i)                            # R(x−1)² = R(x)² up(x−1)²
                fv = _leave_gauge(gn, R2, xdim ? _xdim_pass(Q, X[i-1] + 1) : nothing)
                if i - 1 == m
                    _divide!(w, m+1:n, fv)                               # backward f(m) = 1 as well
                else
                    w[i-1] = fv
                    gp, gc = gc, gn
                    if !(2.0^-600 <= _wmag(R2[1]) <= 2.0^600)
                        s2 = _wmag(R2[1]) > 1 ? 2.0^-300 : 2.0^300
                        gp = _dws(gp, s2); gc = _dws(gc, s2); R2 = _gauge_rescale(R2, s2)
                    end
                    if abs(fv[1]) > 2.0^500
                        _shrink!(w, i-1:n)
                        gp = _dws(gp, 2.0^-500); gc = _dws(gc, 2.0^-500)
                    end
                end
            end
        end
    end
    # ---- normalise: by orthogonality on the real axis, from a certified entry at complex q ----
    sc = _dwm(cfac, _column_scale(Q, W, w, X, n, J2, J3, L1, L2, L3, q, U, Z, xdim))
    @inbounds for i in 1:n
        v = _dwm(w[i], sc)
        w[i] = v
        out === nothing || (out[i] = v[1] + v[2])
    end
    return X
end

"""
    _column_scale(Q, W, w, X, n, J2, J3, L1, L2, L3, q, U, Z, xdim)

The one overall factor the recurrence leaves free, and the two ways of pinning it.

**On the real axis and at a level**, orthogonality does it: `Σ_x [2x+1] f(x)² = 1/[2L1+1]`, so the scale is
one square root of a sum of *positive* terms, and only the sign is left for a certified entry to settle.
Nothing cancels and the result is as accurate as the entries. That path is unchanged.

**At complex q the terms are not positive.** `Σ [2x+1] f²` then cancels — measured `4.6e6` on a j = 4
column at `q = e^{0.108iπ}`, and far more nearer a root of unity — and the cancellation multiplies the
entries' own relative error straight into the scale, hence into every entry. That column came out `7.7e−10`
from a 512-bit reference, uniformly, while `q6j` of the same symbols was accurate to `6e−15`: the
recurrence was right and the normalisation threw it away.

So at complex q the scale comes from the certified symbol itself, at the entry where the column peaks:
`f(x*) / w(x*)` is well conditioned by construction, costs the one `analytic_value` call that
`_sign_reference` was already making to settle the sign, and pins magnitude and phase together rather than
pinning the magnitude badly and then correcting the sign. Orthogonality remains the fallback if the peak
entry's certified value is zero or not finite.
"""
function _column_scale(Q, ::Type{W}, w, X, n::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int,
                       q, U, Z, xdim::Bool) where {W}
    if W === CDWord
        ie = 1
        @inbounds for i in 2:n
            abs(w[i][1]) > abs(w[ie][1]) && (ie = i)
        end
        if !iszero(w[ie][1])
            ref = _entry_value(Q, X[ie], J2, J3, L1, L2, L3)
            tgt = _phase_ref(W, ref, xdim ? _xdim_pass(Q, X[ie] + 1) : nothing)
            if isfinite(tgt) && !iszero(tgt)
                return _dwdiv((tgt, zero(tgt)), w[ie])
            end
        end
    end
    tot = Z
    @inbounds for i in 1:n
        sq = _dwm(w[i], w[i])
        tot = _dwa(tot, xdim ? sq : _dwm(q(X[i] + 1), sq))   # Σ [2x+1] f² either way
    end
    s0 = _dwdiv(U, _dwsqrt(_dwm(q(L1 + 1), tot)))
    ie, ref = _sign_reference(w, X, J2, J3, L1, L2, L3, Q)
    return _fix_phase(s0, w[ie], _phase_ref(W, ref, xdim ? _xdim_pass(Q, X[ie] + 1) : nothing))
end

"""
Leave the gauge: f = h/P, optionally with the factor √[2x+1] folded into the same square root. Passing the
`q`-integer rather than multiplying afterwards saves a double-word square root per entry.
"""
@inline _leave_gauge(h::DWord, P2::DWord, ::Nothing) = _dwdiv(h, _dwsqrt(P2))
@inline _leave_gauge(h::DWord, P2::DWord, xq::DWord) = _dwm(h, _dwsqrt(_dwdiv(xq, P2)))

# At complex q the running gauge is P itself, not P²: see `_gauge_step`.
@inline _leave_gauge(h::CDWord, P::CDWord, ::Nothing) = _dwdiv(h, P)
@inline _leave_gauge(h::CDWord, P::CDWord, sq::CDWord) = _dwm(_dwdiv(h, P), sq)

"""
The `√[2x+1]` that `xdim` folds into each entry. It varies along the column, so unlike the overall scale
its branch is not something the final sign fix can absorb: it has to be the package's, i.e. split over the
Ψ basis. On the real axis it is the positive root of a positive number and nothing changes — the real path
still hands `_leave_gauge` the *unrooted* `[2x+1]` so that one square root covers both it and the gauge.
"""
@inline _xdim_pass(Q, n::Int) = _qint_dw(Q, n)
@inline _xdim_pass(Q::ComplexQ, n::Int) = _sqrt_qint(Q, n)
@inline _xdim_root(Q, n::Int) = _dwsqrt(_qint_dw(Q, n))
@inline _xdim_root(Q::ComplexQ, n::Int) = _sqrt_qint(Q, n)

"""
    _gauge_step(P, e) -> P

Advance the gauge by one step. On the real axis the running quantity is `P²`, multiplied by `up(x)² = e`,
and the square root is taken once when an entry leaves the gauge — `√` of a positive number, so nothing
can go wrong. That is left exactly as it was.

At complex `q` it can, and did. `√` is not multiplicative across its cut, so the principal root of an
*accumulated* `P²` is not the accumulation of the roots: as the product's argument wanders past ±π the
root jumps, and the entries after the jump come out with the opposite sign. Measured on
`q = e^{0.2718iπ}`, the first entry of a column matched `q6j` and the rest were negated; rooting `e` as a
whole fixed that case and broke another, where two consecutive steps flipped and only the entry *between*
them came out wrong. So the complex path accumulates `P ← P·up(x)` with `up(x)` assembled from the roots
of its **eight individual q-integers**, `√Dn·√Bn/[2x]` factor by factor — exactly what
`analytic_rules.jl` does for the prefactor. The branch then depends only on which factors occur, and both
the column and `q6j` see the same ones.
"""
@inline _gauge_step(P::DWord, work::ColumnWork, i::Int) = _dwm(P, work.e[i])
@inline _gauge_step(P::CDWord, work::ColumnWork, i::Int) = _dwm(P, work.up[i])

"Rescale the running gauge to match a rescaling of `h` by `sc` — `P²` needs `sc²`, `P` needs `sc`."
@inline _gauge_rescale(P::DWord, sc::Float64) = _dws(P, sc * sc)
@inline _gauge_rescale(P::CDWord, sc::Float64) = _dws(P, sc)

"""
Resolve the one sign the orthogonality normalisation leaves free. Orthogonality fixes `|f|`; the overall
sign comes from one certified entry. On the real axis the scale is positive, so comparing the *unscaled*
entry with the reference settles it, and that is left exactly as it was. At complex `q` the scale carries
a phase of its own, so the comparison has to be made after it is applied — and "same sign" becomes
"same half-plane", `Re(conj(ref)·v) ≥ 0`, which is the only ambiguity a square root can introduce.
"""
@inline _fix_phase(s0::DWord, w_ie::DWord, ref) =
    signbit(ref) == signbit(w_ie[1]) ? s0 : _dwn(s0)

"""
The reference an entry of the column is actually compared against. `_sign_reference` returns the bare 6j
symbol, but with `xdim = true` the column carries `√[2x+1]` as well. On the real axis that factor is
positive and invisible to a sign test; at complex `q` it is a phase, and leaving it out rotates the
comparison — measured, it flipped whole columns of an F-matrix. The root is the one the column folded in,
[`_xdim_pass`](@ref), not a fresh principal root of the product.
"""
@inline _phase_ref(::Type{DWord}, ref, _) = ref
@inline _phase_ref(::Type{CDWord}, ref, ::Nothing) = ComplexF64(ref)
@inline _phase_ref(::Type{CDWord}, ref, sq::Tuple) = ComplexF64(ref) * sq[1]
@inline function _fix_phase(s0::CDWord, w_ie::CDWord, ref)
    v = _dwm(w_ie, s0)[1]
    return real(conj(ComplexF64(ref)) * v) >= 0 ? s0 : _dwn(s0)
end

"Scale entries by 2⁻⁵⁰⁰ (exact; entries far below the column's peak may underflow)."
@inline function _shrink!(w, rng)
    @inbounds for t in rng
        w[t] = _dws(w[t], 2.0^-500)
    end
    return nothing
end

"Divide entries by d (one double-word reciprocal, then products)."
@inline function _divide!(w, rng, d::Tuple{T,T}) where {T<:Number}
    r = _dwdiv(_dwone(Tuple{T,T}), d)
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

"""
At complex `q` no shortcut applies: the q-integers are not positive, so a one-term sum carries the phase of
its factorials rather than a bare `(−1)^z`. The certified value is the reference.
"""
_entry_sign(Q::ComplexQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int) =
    _entry_value(Q, X2, J2, J3, L1, L2, L3)

function _entry_sign(Q, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    s = sixj_sum(X2, J2, J3, L1, L2, L3)
    is_empty_sum(s) && return 0.0
    if s.zlo == s.zhi && (Q isa ClassicalQ || Q isa RealQ || _valuation_free(s, Q.k + 2))
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
"One certified symbol at generic real q, for the column's overall sign."
_entry_value(Q::RealQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)::Float64 =
    real(analytic_value(sixj_sum(X2, J2, J3, L1, L2, L3), Q.q))
"And at complex q, where the reference is a phase rather than a sign."
_entry_value(Q::ComplexQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)::ComplexF64 =
    ComplexF64(analytic_value(sixj_sum(X2, J2, J3, L1, L2, L3), Q.q))

# ---------------------------------------------------------------------------------
#  One hard symbol from the nearer end of its column
#
#  A symbol whose Racah sum has lost its digits (κ beyond the compensated range) is one entry of a column,
#  and a column recurrence can avoid that cancellation. It has its own directional stability constraints.
#  We compute only the stretch between a seed and the target: O(distance) arithmetic, O(n) workspace.
#
#  Two seeds are needed for a three-term recurrence, and at a column end the boundary condition supplies the
#  second for free: E vanishes there, so the recurrence degenerates to f(edge+1) = −di(edge) f(edge)/up(edge).
#  One accurate value at the end is therefore enough — and at the end the Racah sum usually has a single
#  term, so that value avoids sum cancellation (and is taken in split form, mantissa and binary exponent, because
#  at large spin an end value can be far below the smallest Float64 while the target is of order one).
#
#  The number of Racah terms suggests candidate interior seeds. Their actual compensated error bounds
#  must pass a separate check and are transported through the remaining recurrence. Term count alone
#  does not determine conditioning. Two interior seeds need no boundary condition.
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
Empirical term-count threshold for proposing an interior seed. Measurements found
the compensated boundary at 60–82 terms over k = 1000–2000 and spins 100–250
(`dev/prototypes/check_turning_points.jl`). This is a placement heuristic;
the actual seed error bound must still pass and be propagated to the target.
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

"""
Steps an interior-seeded run may take. A boundary seed is a product of split table entries, accurate to a
double word, so its uncertainty is negligible and the run's length is limited only by the rise. An interior
seed is a *compensated sum*, and its certified bound is ~2e-16 relative — measured, and irreducible here:
half of it is the final rounding to `Float64`, the other half the compensated pass's own bound, so carrying
the double word through would gain a factor of two, not the factors of 10¹⁴ that transporting the error over
a long run costs. With `rtol = 1e-14` the seed's error may be amplified by about 50, i.e. the run may rise
about 5 bits, which at the 0.2–0.5 bits per step observed in a forbidden stretch is a few dozen steps.

Past that the transport fails *after* the run has been paid for, and the fall back to the boundary seed pays
for a second one. Measured over eight escalated symbols, the interior seed won (1.7–1.8×) at 3 and 13 steps
and lost (1.1–2.0×) at 33, 53, 153, 253 and 853; this threshold keeps exactly the cases that win.
"""
const SEED_MAX_STEPS = 24

"Number of Racah terms of one entry of a column (doubled labels), in closed form."
@inline function _nterms(X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    α, β = racah_sums(X2, J2, J3, L1, L2, L3)
    return β[1] - α[4] + 1
end

"""
    _seed_index(X, i, from_left, J...) -> index or 0

An interior-seed candidate closest to the target: the farthest index from the
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

"Compensated seed and its absolute error bound; no modular zero decision or precision escalation."
function _entry_seed(Q::LevelQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    s = sixj_sum(X2, J2, J3, L1, L2, L3)
    is_empty_sum(s) && return nothing
    # For an admissible 6j, terms beyond k have a vanishing numerator;
    # prefactor and denominator arguments in the remaining interval are < h.
    v,B = _sum_compensated(s,(s.zlo:min(s.zhi,Q.k),),Q.k,Q.tab)
    return _entry_seed_result(v,B)
end
function _entry_seed(Q::ClassicalQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    s = sixj_sum(X2, J2, J3, L1, L2, L3)
    is_empty_sum(s) && return nothing
    N = max_argument(s)
    v,B = _sum_compensated(s, (s.zlo:s.zhi,), 0, classical_tables(Float64, N))
    return _entry_seed_result(v,B)
end
_entry_seed_result(v,B) = isfinite(v) && abs(v)>=floatmin(Float64) && isfinite(B) ?
    (value=v,bound=max(B,eps(abs(v)))) : nothing

"A value in split form: mantissa (a double word, real or complex) times 2^exp, so tiny end values do not underflow."
struct SplitValue{W}
    m::W
    e::Int
end

_value(v::SplitValue) = _wldexp(v.m[1] + v.m[2], v.e)

"""
    _edge_split(Q, X2, J2, J3, L1, L2, L3) -> SplitValue or nothing

The symbol at one end of a column, in split form. A one-term Racah sum is a product
of split-table entries, retaining tiny seeds without Float64 underflow. Otherwise
the single-symbol kernel is used; nonfinite and subnormal results are rejected.
"""
function _edge_split(Q::LevelQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    s = sixj_sum(X2, J2, J3, L1, L2, L3)
    is_empty_sum(s) && return nothing
    tab = Q.tab
    if s.zlo == s.zhi && _valuation_free(s, Q.k + 2)          # single term: exact product of table entries
        return _one_term_split(s,tab)
    end
    v = _entry_value(Q, X2, J2, J3, L1, L2, L3)
    return (!isfinite(v) || abs(v)<floatmin(Float64)) ? nothing : SplitValue((v, 0.0), 0)
end

"One-term seed retained as double-word mantissa and binary exponent, even below Float64 range."
function _one_term_split(s,tab)
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

function _edge_split(Q::ClassicalQ, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    s=sixj_sum(X2,J2,J3,L1,L2,L3)
    is_empty_sum(s) && return nothing
    s.zlo==s.zhi && return _one_term_split(s,classical_tables(Float64,max_argument(s)))
    v = _entry_value(Q, X2, J2, J3, L1, L2, L3)
    return (!isfinite(v) || abs(v)<floatmin(Float64)) ? nothing : SplitValue((v, 0.0), 0)
end

# ---------------------------------------------------------------------------------
#  Seeds at generic q
#
#  The near-edge tier needs two things the generic-q route did not expose: a value in split form (an end
#  entry can sit far below the smallest Float64 while the target is of order one — narrowing it first
#  loses the symbol), and an honest bound on it. Both come straight out of the double-word analytic pass,
#  whose mantissa *is* a double word and whose exponent is already carried separately. So a seed off the
#  level axis costs one pass and nothing else.
# ---------------------------------------------------------------------------------

# The double-word mantissa of an analytic pass, as this file's `DWord`/`CDWord`. Left untyped because
# `DWNum` is defined in `analytic_rules.jl`, which `workspace.jl` forces to be included after this file.
@inline _dword(m) = (m.hi, m.lo)
@inline _dword(m::Complex) =
    (ComplexF64(real(m).hi, imag(m).hi), ComplexF64(real(m).lo, imag(m).lo))

"""
    _analytic_split(s, q) -> (SplitValue, relative bound) or nothing

The double-word pass in split form. `nothing` when the sum has lost its own digits: a seed must be far
better than the run it starts, so a cancelling one is refused rather than transported.
"""
function _analytic_split(s::FactorialSum, q)
    is_empty_sum(s) && return nothing
    N = max_argument(s)
    v, kappa, = _analytic_pass(s, _analytic_table(_dwnum(q), N, nothing))
    (isfinite(kappa) && !iszero(v.m)) || return nothing
    rel = Float64(_SEED_OPS * eps(DWNum) * kappa)
    (isfinite(rel) && rel <= SEED_MAX_REL) || return nothing
    return SplitValue(_dword(v.m), v.e), rel
end

"Operation count charged to a double-word seed pass; generous, since a seed is cheap and a bad one is not."
const _SEED_OPS = 1024.0

"How inaccurate a generic-q seed may be before the recurrence is not worth starting from it."
const SEED_MAX_REL = 1e-22

function _entry_seed(Q::Union{RealQ,ComplexQ}, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    r = _analytic_split(sixj_sum(X2, J2, J3, L1, L2, L3), Q.q)
    r === nothing && return nothing
    sv, rel = r
    v = _value(sv)
    (isfinite(v) && !iszero(v) && abs(v) >= floatmin(Float64)) || return nothing
    return (value = v, bound = max(rel * abs(v), eps(abs(v))))
end

function _edge_split(Q::Union{RealQ,ComplexQ}, X2::Int, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    r = _analytic_split(sixj_sum(X2, J2, J3, L1, L2, L3), Q.q)
    r === nothing && return nothing
    sv, _ = r
    return isfinite(_wmag(sv.m[1])) ? sv : nothing
end

"Respect the relative-error budget when a split result reaches the Float64 output range."
function _entry_output(m,e::Int,estimate::Float64,rtol::Float64)
    v=_wldexp(m[1]+m[2],e)
    (!isfinite(v) || iszero(v)) && return nothing
    # Include final rounding explicitly; dividing first avoids underflow in rtol*|v|.
    estimate + eps(abs(v))/abs(v) <= rtol ? v : nothing
end

"""
Coefficients `di` and `e` of the recurrence for column indices `lo:hi` (1-based), into `work` — and, at
complex `q`, `up(x−1)` itself alongside `up(x−1)² = e`.

`up` is not redundant with `e`. Off the real axis the gauge cannot be left by rooting an accumulated `P²`:
`√` is not multiplicative across its cut, so the principal root of the product is not the product of the
roots and the entries come out with the wrong sign — measured as an exact factor of −1 on two of four
sampled complex `q`. `sixj_column!` solved this by accumulating `P ← P·up(x)` with `up` assembled from the
roots of its eight individual q-integers over the Ψ basis, and the near-edge tier needs the same thing;
see [`_gauge_step`](@ref) and [`_sqrt_qint_product`](@ref).
"""
function _coefficients!(work::ColumnWork, Q, X, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int,
                        lo::Int, hi::Int)
    di, e = work.di, work.e
    W = _wordtype(Q)
    Z = _dwzero(W)
    cplx = W === CDWord
    cplx && length(work.up) < length(X) && resize!(work.up, length(X))
    cnt = cplx ? zeros(Int, length(Q.psi) + 1) : Int[]
    q(m) = _qint_dw(Q, m)
    C = _dwm(q((J3 - L1 - L2) ÷ 2), q((J3 + L1 - L2) ÷ 2 + 1))
    bargs_of(x2) = ((x2 - J2 + J3) ÷ 2 + 1, (x2 + J2 + J3) ÷ 2 + 2,
                    (x2 - L2 + L3) ÷ 2 + 1, (L2 + L3 - x2) ÷ 2)
    dargs_of(x2) = ((x2 + J2 - J3) ÷ 2, (J2 + J3 - x2) ÷ 2 + 1,
                    (x2 + L2 - L3) ÷ 2, (x2 + L2 + L3) ÷ 2 + 1)
    prod4(a) = _dwm(_dwm(q(a[1]), q(a[2])), _dwm(q(a[3]), q(a[4])))
    bprev = lo > 1 ? bargs_of(X[lo-1]) : (0, 0, 0, 0)
    Bprev = lo > 1 ? prod4(bprev) : Z
    @inbounds for i in lo:min(hi + 1, length(X))
        x2 = X[i]
        ba = bargs_of(x2)
        B = prod4(ba)
        t = iszero(B[1]) ? Z : _dwm(B, _qinv_dw(Q, x2 + 2))
        if x2 > 0
            da = dargs_of(x2)
            D = prod4(da)
            iv = _qinv_dw(Q, x2)
            t = _dwa(t, _dwm(D, iv))
            if i > 1
                e[i] = _dwm(_dwm(D, Bprev), _dwm(iv, iv))
                cplx && (work.up[i] = _dwm(_sqrt_qint_product(Q,
                    (da[1], da[2], da[3], da[4], bprev[1], bprev[2], bprev[3], bprev[4]), cnt), iv))
            end
        end
        di[i] = _dwa(t, _dwm(q(x2 + 1), C))
        Bprev = B
        bprev = ba
    end
    return nothing
end

"""
The gauge on the `P²` scale, so one rescaling rule serves both word types: the real path accumulates `P²`
and the complex path accumulates `P`, and `eP` counts powers of two in `P²` either way.
"""
@inline _gauge_span(P::DWord) = P[1]
@inline _gauge_span(P::CDWord) = abs(P[1])^2

"`up(x−1)` at index `i` — the Ψ-split one at complex q, the ordinary root of `up²` on the real axis."
@inline _gauge_root(work::ColumnWork{DWord}, i::Int, up2::DWord) = _dwsqrt(up2)
@inline _gauge_root(work::ColumnWork{CDWord}, i::Int, ::CDWord) = @inbounds work.up[i]

"The gauge's own starting value from `up(x−1)`: `P² = up²` on the real axis, `P = up` off it."
@inline _gauge_seed(::ColumnWork{DWord}, up::DWord, up2::DWord) = up2
@inline _gauge_seed(::ColumnWork{CDWord}, up::CDWord, ::CDWord) = up

"""
    sixj_entry(J, k, work; rtol) -> value or nothing

One level-k 6j symbol (doubled labels) by recursion along a short column path.
The arithmetic cost is O(distance), with O(column length) workspace. This avoids
Racah-sum cancellation but has its own directional stability and node sensitivity.
Seed uncertainty and output rounding are checked; coefficient and recurrence
roundoff use a heuristic estimate, not a rigorous accuracy certificate. Returns
`nothing` on unsafe paths or insufficient estimated accuracy, allowing sum fallback.
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

function sixj_entry(J0::NTuple{6,Int}, Q::Union{LevelQ,ClassicalQ,RealQ,ComplexQ},
                    work::ColumnWork = ColumnWork(_wordtype(Q)); rtol::Float64 = 1e-14)
    isfinite(rtol) && rtol>=0 || throw(ArgumentError("rtol must be finite and nonnegative"))
    J, _ = _best_column(J0, Q)                                # the shortest run among the six columns
    return _sixj_entry_path(J,Q,_column_work(Q,work);rtol=rtol)
end

"""
    sixj_entry(J, q::Number, work; rtol) -> value or nothing

The same tier at generic `q`, which is what makes it reachable from `analytic_value`. The recurrence never
mentioned a level; what was missing was a table (`RealQ`/`ComplexQ` supply it) and seeds in split form
(`_analytic_split` does). The table is sized for the worst of the six candidate columns, since
`_best_column` has not chosen yet.
"""
function sixj_entry(J::NTuple{6,Int}, q::Number, work::ColumnWork = ColumnWork(); rtol::Float64 = 1e-14)
    N = sum(J) ÷ 2 + 4                       # every q-integer index `_coefficients!` can reach
    qq = _analytic_q(q)
    Q = qq isa Real ? real_q_tables(qq, N) : complex_q_tables(ComplexF64(qq), N)
    return sixj_entry(J, Q, work; rtol = rtol)
end

"Internal path kernel. Explicit direction/seed switches are for stability experiments."
function _sixj_entry_path(J::NTuple{6,Int},Q,work::ColumnWork=ColumnWork();
                         rtol::Float64=1e-14,from_left::Union{Nothing,Bool}=nothing,
                         interior_seeds::Bool=true)
    W = _wordtype(Q)
    X2, J2, J3, L1, L2, L3 = J
    X = _column_range(Q, J2, J3, L1, L2, L3)
    n = length(X)
    n == 0 && return nothing
    d = X2 - first(X)
    (iseven(d) && 0 <= d <= last(X) - first(X)) || return nothing
    i = d ÷ 2 + 1
    n == 1 && return nothing                                  # the symbol *is* its own edge
    _resize!(work, n)
    from_left = isnothing(from_left) ? i - 1 <= n - i : from_left
    # Where to start. Two certified entries near the end of the certified region start the recursion directly
    # and shorten the run; failing that, one value at the column end does, because the boundary condition
    # supplies the second. `s == 0` means "use the end".
    s = interior_seeds ? _seed_index(X, i, from_left, J2, J3, L1, L2, L3) : 0
    steps_end = from_left ? i - 1 : n - i
    if s != 0
        steps_seed = from_left ? i - s : s - i
        # Worth the two seed evaluations, and short enough that the seeds' own 2e-16 uncertainty survives
        # being transported to the target. A long interior run fails that test only after paying for itself.
        (steps_end - steps_seed >= SEED_MIN_SAVING && steps_seed <= SEED_MAX_STEPS) || (s = 0)
    end
    seedm = seedc = _dwzero(W)                                # f at the two seed entries
    errm = errc = 0.0
    if s != 0
        s2 = from_left ? s - 1 : s + 1
        v1 = _entry_seed(Q, X[s], J2, J3, L1, L2, L3)
        v2 = v1 === nothing ? nothing : _entry_seed(Q, X[s2], J2, J3, L1, L2, L3)
        if v1 === nothing || v2 === nothing || v1.bound>rtol/8*abs(v1.value) || v2.bound>rtol/8*abs(v2.value)
            s = 0                                             # not usable: fall back to the end
        else
            seedc = (v1.value, zero(v1.value)); seedm = (v2.value, zero(v2.value))
            errc=v1.bound; errm=v2.bound
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
    hm = _dwzero(W); hc = _dwone(W); P2 = _dwone(W); eP = 0; eH = 0
    fv = _dwone(W); fexp = 0; fmax_log = 0
    if s != 0
        # start at the seed pair: with P = 1 at the outer seed, h = f there and h = f·up at the inner one,
        # where up² is the coefficient the recursion already carries
        ixs = from_left ? s : s+1
        up2 = e[ixs]
        _w_has_root(up2) || return nothing
        up = _gauge_root(work, ixs, up2)
        hm = seedm; hc = _dwm(seedc, up); P2 = _gauge_seed(work, up, up2)
        errc=nextfloat(errc*nextfloat(abs(up[1])+abs(up[2])))
        fv = seedc
        fmax_log = max(exponent(abs(seedm[1])), exponent(abs(seedc[1])))
    end
    steps = from_left ? ((s == 0 ? 1 : s):i-1) : ((s == 0 ? n : s):-1:i+1)
    isempty(steps) && return s == 0 ? nothing : _entry_output(seedc,0,4eps(Float64),rtol)
    oscillatory=false
    @inbounds for j in steps
        if 1<j<n && W === DWord
            # Local characteristic discriminant of the symmetric recurrence:
            # d² - 4*up(x-1)*up(x). Past an oscillatory interval, entering a
            # forbidden interval can amplify the unwanted dominant solution.
            # Reject this direction instead of trusting its final magnitude.
            #
            # This is a real-axis statement — "classically allowed" and "forbidden" are what the two
            # regions of a 6j column are — and it does not survive off it: with complex coefficients the
            # magnitudes wander in and out of the inequality and the guard aborts almost every run
            # (measured: every column past j = 10 at `q = 0.8 + 0.3im`). At complex `q` the growth of an
            # unwanted solution is caught instead by `rise_log`, which measures the thing the guard is a
            # proxy for — how far the path rose above the target — and is checked below for every word
            # type. Validated on a sweep: see `dev/results/generic_q_path.md`.
            discriminant_scale=4sqrt(max(0.0,e[j][1]*e[j+1][1]))
            if di[j][1]^2<=discriminant_scale
                oscillatory=true
            elseif oscillatory
                return nothing
            end
        end
        if s != 0
            # Propagate seed uncertainty through the actual recurrence. Outward
            # rounding protects this nonnegative envelope from Float64 rounding.
            ix=from_left ? j : j+1
            a=nextfloat(abs(di[j][1])+abs(di[j][2]))
            b=nextfloat(abs(e[ix][1])+abs(e[ix][2]))
            errn=nextfloat(nextfloat(a*errc)+nextfloat(b*errm))
            errm=errc; errc=errn
        end
        hn = _dwn(_dwm(di[j], hc))
        if from_left
            j > 1 && (hn = _dwa(hn, _dwn(_dwm(e[j], hm))))
            P2 = _gauge_step(P2, work, j+1)
        else
            j < n && (hn = _dwa(hn, _dwn(_dwm(e[j+1], hm))))
            P2 = _gauge_step(P2, work, j)
        end
        hm = hc; hc = hn
        span = _gauge_span(P2)
        if !(2.0^-600 <= span <= 2.0^600)
            grew = span > 1
            P2 = _gauge_rescale(P2, grew ? 2.0^-300 : 2.0^300)   # 2^∓600 in P², whichever P2 holds
            eP += grew ? 600 : -600
        end
        if abs(hc[1]) > 2.0^500
            hm = _dws(hm, 2.0^-500); hc = _dws(hc, 2.0^-500); eH += 500
            if s != 0
                errm=nextfloat(ldexp(errm,-500)); errc=nextfloat(ldexp(errc,-500))
            end
        end
        # Only the last entry is wanted, so the gauge is left every SAMPLE steps (and at the end), just to
        # track how far the path rises above the target — the factor in the error estimate.
        if j == last(steps) || rem(j, SAMPLE) == 0
            fv = _leave_gauge(hc, P2, nothing)
            fexp = eH - eP ÷ 2
            iszero(fv[1]) || (fmax_log = max(fmax_log, exponent(abs(fv[1])) + fexp))
        end
    end
    iszero(fv[1]) && return nothing
    # How far the path rose above the target. A large rise means the target is near a node (or is an exact
    # zero), and then this tier must not answer at all, whatever its rounding estimate says.
    rise_log = fmax_log - (exponent(abs(fv[1])) + fexp)
    rise_log > MAX_RISE_LOG2 && return nothing
    # the double-word error of the run, plus the floor its table imposes (zero at a level and
    # classically), both relative to the target entry and both magnified by the rise
    est = (length(steps) * 2.0^-104 + _table_floor(Q)) * exp2(rise_log)
    # Interior seed uncertainty propagated through the rounded coefficients.
    s == 0 || (est += errc/abs(hc[1]+hc[2]))
    if !isfinite(est) || est > rtol
        # Interior seeds may be too uncertain even when the split boundary seed
        # supports an accurate run. Retry it once before escalating the sum.
        return s==0 ? nothing : _sixj_entry_path(J,Q,work;rtol=rtol,
            from_left=from_left,interior_seeds=false)
    end
    if s != 0
        return _entry_output(fv,fexp,est,rtol)               # the seeds are already normalised values
    end
    val = _dwm(fv, seed.m)
    iszero(val[1]) && return nothing
    return _entry_output(val,seed.e+fexp,est,rtol)
end
