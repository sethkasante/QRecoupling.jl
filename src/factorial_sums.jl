
# ---------------------------------------------------------------------------------
#  Factorial sums: the rule of a terminating q-hypergeometric sum
#
#  value = σ₀ · (√) Π_n [n]!^{c_n} · Σ_{z = zlo}^{zhi} (±1)^z Π_i [a_i z + b_i]!^{c_i}
#
#  A DCR expands such a sum into one cyclotomic monomial per term ratio. The rule itself is a few
#  integers (16 prefactor and 8 term factorials for a 6j symbol, at any spin), and it answers the
#  questions asked before evaluating anything:
#
#  * Valuations. The Φ_h exponent of [n]! is ⌊n/h⌋ (Legendre), so a term's valuation is
#    V(z) = Σ_i c_i ⌊(a_i z + b_i)/h⌋, piecewise constant between points where an argument crosses a
#    multiple of h. Poles, vanishing terms and contributing z-ranges follow without expanding.
#    Convention (as for DCR projections): with the prefactor under a square root, the doubled valuation
#    of a term is v_pre + 2V(z); negative is a pole, positive means the term vanishes.
#  * Cancellation zeros. The contributing sum vanishes in ℚ(ζ) iff it vanishes at every Galois
#    conjugate; testing two conjugates modulo one prime errs with probability about p⁻².
#  * Numbers. A ratio loop over q-integer tables on the contributing ranges only.
# ---------------------------------------------------------------------------------


"The affine factorial [a z + b]!^c."
struct AffineFactorial
    a::Int32
    b::Int32
    c::Int32
end

@inline _arg(f::AffineFactorial, z::Int) = Int(f.a) * z + Int(f.b)

"""
    FactorialSum

`sign0 · (√)Π pre · Σ_{z=zlo}^{zhi} (−1)^{z·alternating} Π fac(z)`, where `pre` holds pairs n => c meaning
[n]!^c and sits under a square root when `sqrt_pre` is set.
"""
struct FactorialSum
    pre::Vector{Pair{Int32,Int32}}
    sqrt_pre::Bool
    fac::Vector{AffineFactorial}
    alternating::Bool
    sign0::Int8
    zlo::Int
    zhi::Int
end

const EMPTY_FACTORIAL_SUM = FactorialSum(Pair{Int32,Int32}[], false, AffineFactorial[], false, Int8(0), 0, -1)

is_empty_sum(s::FactorialSum) = s.sign0 == 0 || s.zlo > s.zhi


# ---- recoupling symbols as factorial rules (doubled labels) ----

@inline function _sort4(a, b, c, d)
    a, b = minmax(a, b); c, d = minmax(c, d)
    a, c = minmax(a, c); b, d = minmax(b, d)
    b, c = minmax(b, c)
    return a, b, c, d
end

@inline function _sort3(a, b, c)
    a, b = minmax(a, b); b, c = minmax(b, c); a, b = minmax(a, b)
    return a, b, c
end

"Sorted triangle sums α and quadrilateral sums β of a 6j symbol (doubled labels); they fix its value."
@inline function racah_sums(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    α = _sort4((J1 + J2 + J3) ÷ 2, (J1 + J5 + J6) ÷ 2, (J2 + J4 + J6) ÷ 2, (J3 + J4 + J5) ÷ 2)
    β = _sort3((J1 + J2 + J4 + J5) ÷ 2, (J1 + J3 + J4 + J6) ÷ 2, (J2 + J3 + J5 + J6) ÷ 2)
    return α, β
end

"The 6j symbol as a factorial rule."
function sixj_sum(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    _δtet(J1, J2, J3, J4, J5, J6) || return EMPTY_FACTORIAL_SUM
    α, β = racah_sums(J1, J2, J3, J4, J5, J6)
    pre = Pair{Int32,Int32}[]
    for b in β, a in α
        push!(pre, Int32(b - a) => Int32(1))   # [β_j − α_i]!
    end
    for a in α
        push!(pre, Int32(a + 1) => Int32(-1))                   # 1/[α_i + 1]!
    end
    fac = AffineFactorial[AffineFactorial(1, 1, 1)]             # [z + 1]!
    for a in α
        push!(fac, AffineFactorial(1, -a, -1))                  # 1/[z − α_i]!
    end
    for b in β
        push!(fac, AffineFactorial(-1, b, -1))                  # 1/[β_j − z]!
    end
    return FactorialSum(pre, true, fac, true, Int8(1), α[4], β[1])
end

"The 3j symbol as a factorial rule, including the Wigner phase (−1)^{j1−j2−m3}."
function threej_sum(J1::Int, J2::Int, J3::Int, M1::Int, M2::Int, M3::Int = -M1 - M2)
    (_δ(J1, J2, J3) && _mproj_ok(J1, J2, J3, M1, M2, M3)) || return EMPTY_FACTORIAL_SUM
    pre = Pair{Int32,Int32}[]
    for n in ((J1 + J2 - J3) ÷ 2, (J1 - J2 + J3) ÷ 2, (-J1 + J2 + J3) ÷ 2,
              (J1 + M1) ÷ 2, (J1 - M1) ÷ 2, (J2 + M2) ÷ 2, (J2 - M2) ÷ 2, (J3 + M3) ÷ 2, (J3 - M3) ÷ 2)
        push!(pre, Int32(n) => Int32(1))
    end
    push!(pre, Int32((J1 + J2 + J3) ÷ 2 + 1) => Int32(-1))
    α1 = (J3 - J2 + M1) ÷ 2; α2 = (J3 - J1 - M2) ÷ 2
    β1 = (J1 + J2 - J3) ÷ 2; β2 = (J1 - M1) ÷ 2; β3 = (J2 + M2) ÷ 2
    fac = AffineFactorial[AffineFactorial(1, 0, -1), AffineFactorial(1, α1, -1), AffineFactorial(1, α2, -1),
                          AffineFactorial(-1, β1, -1), AffineFactorial(-1, β2, -1), AffineFactorial(-1, β3, -1)]
    sign0 = iseven((J1 - J2 - M3) ÷ 2) ? Int8(1) : Int8(-1)
    return FactorialSum(pre, true, fac, true, sign0, max(0, -α1, -α2), min(β1, β2, β3))
end

"Multiply the square-rooted prefactor by the quantum dimensions [J+1] = [J+1]!/[J]! (doubled labels)."
function with_dimensions(s::FactorialSum, Js::Int...)
    is_empty_sum(s) && return s
    pre = copy(s.pre)
    for J in Js
        push!(pre, Int32(J + 1) => Int32(1))
        push!(pre, Int32(J) => Int32(-1))
    end
    return FactorialSum(pre, s.sqrt_pre, s.fac, s.alternating, s.sign0, s.zlo, s.zhi)
end

negated(s::FactorialSum) = FactorialSum(s.pre, s.sqrt_pre, s.fac, s.alternating, -s.sign0, s.zlo, s.zhi)

"F-symbol (−1)^{j1+j2+j4+j5} √([2j3+1][2j6+1]) {6j} as a factorial rule."
function fsymbol_sum(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    s = with_dimensions(sixj_sum(J1, J2, J3, J4, J5, J6), J3, J6)
    return iseven((J1 + J2 + J4 + J5) ÷ 2) ? s : negated(s)
end

"G-symbol √(Π [2j_i+1]) {6j} as a factorial rule."
gsymbol_sum(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int) =
    with_dimensions(sixj_sum(J1, J2, J3, J4, J5, J6), J1, J2, J3, J4, J5, J6)


# ---- valuations at q = e^{iπ/h} ----

"⌊n/h⌋, the Φ_h exponent of [n]!; negative arguments never reach a contributing term."
@inline fld(n::Int, h::Int) = Base.fld(n, h)

"Φ_h exponent of Π [n]!^c."
prefactor_valuation(s::FactorialSum, h::Int) = sum((Int(c) * fld(Int(n), h) for (n, c) in s.pre); init = 0)

"Φ_h exponent of the term at z."
@inline function term_valuation(s::FactorialSum, z::Int, h::Int)
    v = 0
    @inbounds for f in s.fac
        v += Int(f.c) * fld(_arg(f, z), h)
    end
    return v
end

"Points z in (zlo, zhi] where some ⌊(a z + b)/h⌋ changes value."
function valuation_breakpoints(s::FactorialSum, h::Int)
    pts = Int[]
    lo, hi = s.zlo, s.zhi
    for f in s.fac
        a, b = Int(f.a), Int(f.b)
        if a == 1 || a == -1
            r = a == 1 ? mod(-b, h) : mod(b + 1, h)
            z = lo + 1 + mod(r - (lo + 1), h)
            while z <= hi
                push!(pts, z)
                z += h
            end
        elseif a != 0
            for z in lo+1:hi
                fld(_arg(f, z), h) != fld(_arg(f, z - 1), h) && push!(pts, z)
            end
        end
    end
    sort!(pts)
    unique!(pts)
    return pts
end

function _merge_ranges(rs::Vector{UnitRange{Int}})
    out = UnitRange{Int}[]
    for r in rs
        if !isempty(out) && first(r) == last(out[end]) + 1
            out[end] = first(out[end]):last(r)
        else
            push!(out, r)
        end
    end
    return out
end

"""
    classify_at_level(s, k) -> (status, segments)

Status at q = e^{iπ/(k+2)}: `:empty`, `:pole` (some term has negative valuation), `:zero` (every term
vanishes) or `:finite`; `segments` are the z-ranges of contributing terms. Cancellation zeros are not
detected here; see `is_cancellation_zero`.
"""
function classify_at_level(s::FactorialSum, k::Int)
    is_empty_sum(s) && return (:empty, UnitRange{Int}[])
    h = k + 2
    vp = (s.sqrt_pre ? 1 : 2) * prefactor_valuation(s, h)
    starts = vcat(s.zlo, valuation_breakpoints(s, h))
    segs = UnitRange{Int}[]
    for i in eachindex(starts)
        z0 = starts[i]
        z1 = i < length(starts) ? starts[i+1] - 1 : s.zhi
        v2 = vp + 2 * term_valuation(s, z0, h)
        v2 < 0 && return (:pole, UnitRange{Int}[])
        v2 == 0 && push!(segs, z0:z1)
    end
    isempty(segs) && return (:zero, segs)
    return (:finite, _merge_ranges(segs))
end

"Do all factorials needed on the contributing segments lie in the level tables (arguments 0..k+1)?"
function _within_level_tables(s::FactorialSum, segs, k::Int)
    all(p -> 0 <= p.first <= k + 1, s.pre) || return false
    for seg in segs, f in s.fac
        (0 <= _arg(f, first(seg)) <= k + 1 && 0 <= _arg(f, last(seg)) <= k + 1) || return false
    end
    return true
end


# ---- exact zero test: two Galois conjugates modulo one prime ----

"[n]! and 1/[n]! at two conjugates q ↦ r^a of q = e^{iπ/(k+2)}, modulo a prime p ≡ 1 (mod 2(k+2))."
struct LevelZeroTable
    m::Montgomery
    fact::Matrix{UInt64}       # row n+1, one column per conjugate (Montgomery form)
    invfact::Matrix{UInt64}
end

function LevelZeroTable(k::Int)
    h = k + 2
    n2h = 2h
    N = k + 1
    p = primes_1mod(n2h, 1)[1]
    m = Montgomery(p)
    g = to_mont(m, root_of_unity(p, n2h))
    exps = [1]
    a = 3
    while length(exps) < 2
        gcd(a, n2h) == 1 && push!(exps, a)
        a += 2
    end
    fact = Matrix{UInt64}(undef, N + 1, 2)
    invf = similar(fact)
    qv = Vector{UInt64}(undef, N)
    for (c, e) in enumerate(exps)
        r = mont_pow(m, g, e)
        ri = mont_inv(m, r)
        dinv = mont_inv(m, mont_sub(m, r, ri))
        rn = m.one; rin = m.one
        fact[1, c] = m.one
        for n in 1:N
            rn = mont_mul(m, rn, r); rin = mont_mul(m, rin, ri)
            qv[n] = mont_mul(m, mont_sub(m, rn, rin), dinv)        # [n] = (r^n − r^{−n}) / (r − r^{−1})
            fact[n+1, c] = mont_mul(m, fact[n, c], qv[n])
        end
        invf[N+1, c] = mont_inv(m, fact[N+1, c])
        for n in N:-1:1
            invf[n, c] = mont_mul(m, invf[n+1, c], qv[n])
        end
    end
    return LevelZeroTable(m, fact, invf)
end

# ---- per-level table cache with lock-free reads ----
#
# Levels are small integers, so a table lives in slot k+1 of a vector. Readers take one atomic load of
# the whole vector and index it; writers build the table under a lock, publish a grown copy atomically,
# and never mutate a vector a reader can hold. This matters once many threads evaluate across levels,
# where a lock per lookup was the bottleneck.

mutable struct LevelCache{V}
    @atomic slots::Vector{Union{Nothing,V}}
    lock::ReentrantLock
end

LevelCache{V}() where {V} = LevelCache{V}(Vector{Union{Nothing,V}}(nothing, 0), ReentrantLock())

function get_level!(build::F, c::LevelCache{V}, k::Int) where {F,V}
    sl = @atomic c.slots
    if k + 1 <= length(sl)
        @inbounds t = sl[k+1]
        t !== nothing && return t::V
    end
    @lock c.lock begin
        sl = @atomic c.slots
        if k + 1 <= length(sl)
            @inbounds t = sl[k+1]
            t !== nothing && return t::V
        end
        v = build()::V
        fresh = if k + 1 > length(sl)                 # grow only when the level is out of range
            g = Vector{Union{Nothing,V}}(nothing, max(k + 1, 2 * length(sl), 16))
            copyto!(g, 1, sl, 1, length(sl))
            g
        else
            copy(sl)                                  # never mutate a vector a reader may hold
        end
        fresh[k+1] = v
        @atomic c.slots = fresh
        return v
    end
end

function Base.empty!(c::LevelCache{V}) where {V}
    @lock c.lock begin
        @atomic c.slots = Vector{Union{Nothing,V}}(nothing, 0)
    end
    return c
end

const LEVEL_ZERO_TABLES = LevelCache{LevelZeroTable}()

level_zero_table(k::Int) = get_level!(() -> LevelZeroTable(k), LEVEL_ZERO_TABLES, k)

"""
    is_cancellation_zero(s, segs, k) -> Bool or nothing

Whether the contributing sum vanishes exactly, from two Galois conjugates modulo one prime (false
positives with probability about p⁻² ≈ 2⁻¹²⁴). `nothing` if a factorial falls outside the level tables.
"""
is_cancellation_zero(s::FactorialSum, segs, k::Int) = is_cancellation_zero(s, segs, k, level_zero_table(k))

function is_cancellation_zero(s::FactorialSum, segs, k::Int, tab::LevelZeroTable)
    m = tab.m
    @inbounds for c in axes(tab.fact, 2)
        acc = UInt64(0)
        for seg in segs, z in seg
            t = m.one
            for f in s.fac
                n = _arg(f, z)
                (0 <= n <= k + 1) || return nothing
                b = f.c > 0 ? tab.fact[n+1, c] : tab.invfact[n+1, c]
                for _ in 1:abs(f.c)
                    t = mont_mul(m, t, b)
                end
            end
            acc = (s.alternating && isodd(z)) ? mont_sub(m, acc, t) : mont_add(m, acc, t)
        end
        iszero(acc) || return false
    end
    return true
end

"Is the value exactly zero at level k (empty, all terms vanishing, or cancelling)?"
is_zero_at_level(s::FactorialSum, k::Int) = is_zero_at_level(s, k, level_zero_table(k))

function is_zero_at_level(s::FactorialSum, k::Int, tab::LevelZeroTable)
    st, segs = classify_at_level(s, k)
    (st === :empty || st === :zero) && return true
    st === :pole && return false
    return is_cancellation_zero(s, segs, k, tab) === true
end


# ---- numbers at q = e^{iπ/h} ----

struct QIntTables{T}
    q::Vector{T}      # [n], n = 1..N, correctly rounded from a wide value
    qi::Vector{T}     # 1/[n], correctly rounded
    ql::Vector{T}     # low parts: [n] = q + ql and 1/[n] = qi + qil to about u²
    qil::Vector{T}
    fm::Vector{T}     # [n]! = (fm + fml) · 2^fe, fm ∈ [1/2, 1) correctly rounded
    fe::Vector{Int}
    fml::Vector{T}
    gm::Vector{T}     # 1/[n]! = (gm + gml) · 2^ge, rounded the same way, so inverses are multiplies
    ge::Vector{Int}
    gml::Vector{T}
    classical::Bool   # [n] = n: term ratios are ratios of integers
end

"""
    split_factorials(T, qhi) -> (fm, fe)

Split-exponent factorials from q-integers given in a wider `BigFloat` precision: the running product is
kept exactly as `BigFloat` and each [n]! is rounded to `T` once. Products of ~20 such entries keep ~20 ε of
relative accuracy at any size, where log tables lose |log [n]!| · ε.
"""
function split_factorials(::Type{T}, qhi::Vector{BigFloat}) where {T}
    N = length(qhi)
    ms = Vector{BigFloat}(undef, N + 1); fe = Vector{Int}(undef, N + 1)
    gs = Vector{BigFloat}(undef, N + 1); ge = Vector{Int}(undef, N + 1)
    f = one(BigFloat)
    for n in 0:N
        n > 0 && (f *= qhi[n])
        ms[n+1], fe[n+1] = frexp(f)
        gs[n+1], ge[n+1] = frexp(inv(f))
    end
    return ms, fe, gs, ge
end

"Working precision for building tables of type T: 64 bits beyond what T holds."
_table_precision(::Type{BigFloat}) = precision(BigFloat) + 64
_table_precision(::Type{T}) where {T<:AbstractFloat} = precision(T) + 64

"""
Build `QIntTables{T}` from a generator of the q-integers in wide precision. `qgen()` runs inside the wide
precision (64 guard bits) and returns [1], …, [N] as `BigFloat`. Every stored entry is the correct rounding
of its wide value, and each has a low part so that entry + low part is accurate to about u².
"""
function _qint_tables_from(::Type{T}, qgen::F; classical::Bool = false) where {T,F}
    qw, qiw, ms, fe, gs, ge = setprecision(BigFloat, _table_precision(T)) do
        qhi = qgen()
        ms, fe, gs, ge = split_factorials(T, qhi)
        qhi, BigFloat[inv(x) for x in qhi], ms, fe, gs, ge
    end
    hi(x) = T(x)
    lo(x) = T(x - T(x))                          # the wide value minus its rounding, then rounded
    return QIntTables{T}(hi.(qw), hi.(qiw), lo.(qw), lo.(qiw),
                         hi.(ms), fe, lo.(ms), hi.(gs), ge, lo.(gs), classical)
end

"[1], …, [N] at q = e^{iπ/h} in the current BigFloat precision, by the Chebyshev recurrence
[n+1] = 2cos(π/h)[n] − [n−1] (errors grow linearly, absorbed by the 64 guard bits)."
function _qints_wide(h::Int, N::Int)
    x = Vector{BigFloat}(undef, N)
    N >= 1 && (x[1] = one(BigFloat))
    N >= 2 && (x[2] = 2 * cos(big(pi) / h))
    c2 = N >= 2 ? x[2] : zero(BigFloat)
    for n in 3:N
        x[n] = c2 * x[n-1] - x[n-2]
    end
    return x
end

function QIntTables(::Type{T}, k::Int) where {T}
    h = k + 2
    N = k + 1
    return _qint_tables_from(T, () -> _qints_wide(h, N))
end

const QINT_F64_TABLES = LevelCache{QIntTables{Float64}}()
const QINT_TABLES = Dict{Tuple{DataType,Int,Int},Any}()
const QINT_TABLES_LOCK = ReentrantLock()

qint_tables(::Type{Float64}, k::Int) = get_level!(() -> QIntTables(Float64, k), QINT_F64_TABLES, k)

function qint_tables(::Type{T}, k::Int) where {T}
    key = (T, k, T === BigFloat ? precision(BigFloat) : 0)
    tab = @lock QINT_TABLES_LOCK get!(() -> QIntTables(T, k), QINT_TABLES, key)
    return tab::QIntTables{T}
end

_decimal_digits(::Type{BigFloat}) = precision(BigFloat) * log10(2.0)
_decimal_digits(::Type{T}) where {T<:AbstractFloat} = -log10(eps(T))

"""
Σ over the contributing segments in type T by the term ratio, with the prefactor. Returns the value and
the decimal digits lost to cancellation, log₁₀(max |term| / |Σ|).
"""
_sum_at_level(s::FactorialSum, segs, k::Int, ::Type{T}) where {T} =
    _sum_at_level(s, segs, k, qint_tables(T, k))

"""
Multiply (m, e) by ([n]!)^c using the split tables — multiplies only, no renormalisation. Mantissas lie
in [1/2, 1), so a product of a few dozen of them stays far inside the exponent range; callers renormalise
once with `frexp` after a whole product.
"""
@inline function _split_mul(m::T, e::Int, tab::QIntTables{T}, n::Int, c::Integer) where {T}
    @inbounds if c > 0
        for _ in 1:c
            m *= tab.fm[n+1]; e += tab.fe[n+1]
        end
    else
        for _ in 1:-c
            m *= tab.gm[n+1]; e += tab.ge[n+1]
        end
    end
    return m, e
end

"Renormalise a split value so its mantissa lies in [1/2, 1)."
@inline function _renorm(m, e::Int)
    fr, ex = frexp(m)
    return fr, e + ex
end

"Scale step for the ratio loop: exact powers of two, applied when a term grows past 2^SCALE_BITS."
const SCALE_BITS = 256

"""
    _integer_ratios(s, segs, tab) -> Bool

Classically every term ratio is a ratio of two integers, a/b, built from the factorial arguments. When
both stay below 2^53 they convert to floating point exactly, so a ratio costs one division (one rounding)
instead of K table multiplications (2K roundings), and its low part is exact from one fma.
"""
function _integer_ratios(s::FactorialSum, segs, tab::QIntTables)
    tab.classical || return false
    isempty(segs) && return false
    kn = 0; kd = 0; m = 1
    for f in s.fac
        abs(f.c) == 1 || return false
        (f.a == 1) == (f.c > 0) ? (kn += 1) : (kd += 1)
        for seg in segs
            m = max(m, _arg(f, first(seg)) + 1, _arg(f, last(seg)) + 1)
        end
    end
    return float(m)^max(kn, kd) < 2.0^53
end

"Numerator and denominator of the term ratio t_{z+1}/t_z as integers (sign in the numerator)."
@inline function _int_ratio(s::FactorialSum, z::Int)
    a = s.alternating ? -1 : 1
    b = 1
    @inbounds for f in s.fac
        n = _arg(f, z)
        x = f.a == 1 ? n + 1 : n
        (f.a == 1) == (f.c > 0) ? (a *= x) : (b *= x)
    end
    return a, b
end

"Unit roundoff of T."
_unit(::Type{T}) where {T} = eps(T) / 2

"γ_m = m·u/(1 − m·u), the standard bound on m relative roundings (Inf when m·u ≥ 1/2)."
@inline function _gamma(::Type{T}, m::Integer) where {T}
    mu = m * _unit(T)
    return mu < 0.5 ? T(mu / (1 - mu)) : T(Inf)
end

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
    cstep = intr ? 2 : 2K + 1                      # roundings per step: a/b and t·r, or K entries + K products
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
            elseif buf === nothing
                r = s.alternating ? -one(T) : one(T)
                for f in s.fac
                    n = _arg(f, z)
                    if f.a == 1                                 # [n+1]!/[n]! = [n+1]
                        r *= f.c == 1 ? q[n+1] : f.c == -1 ? qi[n+1] : q[n+1]^f.c
                    elseif f.a == -1                            # [n−1]!/[n]! = 1/[n]
                        r *= f.c == 1 ? qi[n] : f.c == -1 ? q[n] : qi[n]^f.c
                    end
                end
            else                                                # same product, low part kept for later
                r = s.alternating ? -one(T) : one(T); rl = zero(T)
                for f in s.fac
                    n = _arg(f, z)
                    if f.a == 1
                        xh, xl = f.c > 0 ? (q[n+1], ql[n+1]) : (qi[n+1], qil[n+1])
                    else
                        xh, xl = f.c > 0 ? (qi[n], qil[n]) : (q[n], ql[n])
                    end
                    for _ in 1:abs(Int(f.c))
                        p, π = _two_prod(r, xh)
                        rl = fma(rl, xh, fma(r, xl, π)); r = p
                    end
                end
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

# ---- error-free transformations and double-word arithmetic ----

@inline _two_prod(a, b) = (p = a * b; (p, fma(a, b, -p)))
@inline function _two_sum(a, b)
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)
end
@inline _fast_two_sum(a, b) = (s = a + b; (s, b - (s - a)))

"(ah + al)(bh + bl) as a double word."
@inline function _dw_mul(ah, al, bh, bl)
    p, e = _two_prod(ah, bh)
    e = fma(ah, bl, fma(al, bh, e))
    return _fast_two_sum(p, e)
end

"√(h + l) as a double word, h > 0."
@inline function _dw_sqrt(h, l)
    r = sqrt(h)
    p, e = _two_prod(r, r)
    corr = ((h - p) - e + l) / (2r)
    return _fast_two_sum(r, corr)
end

"Multiply (mh + ml)·2^e by ([n]!)^c with double-word split entries."
@inline function _split_mul_dw(mh::T, ml::T, e::Int, tab::QIntTables{T}, n::Int, c::Integer) where {T}
    @inbounds if c > 0
        for _ in 1:c
            mh, ml = _dw_mul(mh, ml, tab.fm[n+1], tab.fml[n+1]); e += tab.fe[n+1]
        end
    else
        for _ in 1:-c
            mh, ml = _dw_mul(mh, ml, tab.gm[n+1], tab.gml[n+1]); e += tab.ge[n+1]
        end
    end
    return mh, ml, e
end

@inline function _renorm_dw(mh, ml, e::Int)
    fr, ex = frexp(mh)
    return fr, ldexp(ml, -ex), e + ex
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
            else
                rh = s.alternating ? -one(T) : one(T); rl = zero(T)
                for f in s.fac
                    n = _arg(f, z)
                    if f.a == 1
                        xh, xl = f.c > 0 ? (q[n+1], ql[n+1]) : (qi[n+1], qil[n+1])
                    else
                        xh, xl = f.c > 0 ? (qi[n], qil[n]) : (q[n], ql[n])
                    end
                    for _ in 1:abs(Int(f.c))
                        p, π = _two_prod(rh, xh)
                        rl = fma(rl, xh, fma(rh, xl, π)); rh = p
                    end
                end
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
        Kr = intr ? 1 : K                                # low-part error per ratio: one rounding, or K
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
carries about 11 significant digits or better; escalated results carry the full target. `fallback()` is used if a factorial falls outside the level tables.
"""
function value_at_level(s::FactorialSum, k::Int, ::Type{T}; fallback) where {T}
    v, status, segs = level_pass1(s, k, qint_tables(T, k))
    status === :done && return v
    status === :fallback && return fallback()
    return level_escalate(s, segs, k, T, level_zero_table(k))
end

"""
    level_pass1(s, k, tab) -> (value, status, segments)

First pass in the table's own type, with no locking and no `BigFloat`, so it is safe to run on many
threads. `status` is `:done` when the value keeps enough digits, `:escalate` when cancellation consumed
them (the caller then runs the zero test and higher precision), or `:fallback` when a factorial falls
outside the level tables.
"""
function level_pass1(s::FactorialSum, k::Int, tab::QIntTables{T}) where {T}
    st, segs = classify_at_level(s, k)
    st === :pole && throw(DomainError(k, "Topological pole at level k=$k."))
    (st === :empty || st === :zero) && return (zero(T), :done, UnitRange{Int}[])
    _within_level_tables(s, segs, k) || return (zero(T), :fallback, segs)
    v, st = _certified_value(s, segs, k, tab)
    return (v, st, segs)
end

"""
    _certified_value(s, segs, k, tab) -> (value, :done | :escalate)

Plain pass with its certified bound; if the bound is too loose, the compensated pass with its own. Thread
safe (no global state). For types other than `Float64` the heuristic digit budget is used, as before.
"""
function _certified_value(s::FactorialSum, segs, k::Int, tab::QIntTables{T}) where {T}
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
    lazy = mode === :lazy || mode === :strict_lazy
    buf = lazy ? Vector{T}(undef, 2 * sum(seg -> length(seg) - 1, segs; init = 0)) : nothing
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
    level_escalate(s, segs, k, T, ztab) -> T

Finishes an `:escalate` case: the modular test first decides whether the value is exactly zero, otherwise
the sum is redone in `BigFloat` with enough bits for the loss, doubling until the digits are there. Changes
the global `BigFloat` precision while it runs, so callers keep it off worker threads.
"""
function level_escalate(s::FactorialSum, segs, k::Int, ::Type{T}, ztab::LevelZeroTable) where {T}
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
