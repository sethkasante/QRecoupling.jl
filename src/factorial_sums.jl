
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

const LEVEL_ZERO_TABLES = Dict{Int,LevelZeroTable}()
const LEVEL_ZERO_TABLES_LOCK = ReentrantLock()

level_zero_table(k::Int) = @lock LEVEL_ZERO_TABLES_LOCK get!(() -> LevelZeroTable(k), LEVEL_ZERO_TABLES, k)

"""
    is_cancellation_zero(s, segs, k) -> Bool or nothing

Whether the contributing sum vanishes exactly, from two Galois conjugates modulo one prime (false
positives with probability about p⁻² ≈ 2⁻¹²⁴). `nothing` if a factorial falls outside the level tables.
"""
function is_cancellation_zero(s::FactorialSum, segs, k::Int)
    tab = level_zero_table(k)
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
function is_zero_at_level(s::FactorialSum, k::Int)
    st, segs = classify_at_level(s, k)
    (st === :empty || st === :zero) && return true
    st === :pole && return false
    return is_cancellation_zero(s, segs, k) === true
end


# ---- numbers at q = e^{iπ/h} ----

struct QIntTables{T}
    q::Vector{T}      # [n], n = 1..k+1
    qi::Vector{T}     # 1/[n]
    lf::Vector{T}     # log [n]! at index n + 1
end

function QIntTables(::Type{T}, k::Int) where {T}
    h = k + 2
    N = k + 1
    s = sin(T(π) / h)
    q = [sin(min(n, h - n) * T(π) / h) / s for n in 1:N]         # [n] = [h − n] avoids rounding near π
    lf = Vector{T}(undef, N + 1)
    lf[1] = zero(T)
    for n in 1:N
        lf[n+1] = lf[n] + log(q[n])
    end
    return QIntTables{T}(q, inv.(q), lf)
end

const QINT_TABLES = Dict{Tuple{DataType,Int,Int},Any}()
const QINT_TABLES_LOCK = ReentrantLock()

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
function _sum_at_level(s::FactorialSum, segs, k::Int, ::Type{T}) where {T}
    tab = qint_tables(T, k)
    lf, q, qi = tab.lf, tab.q, tab.qi
    big_ = T(1e100); small_ = T(1e-100); ln10 = log(T(10))
    lpre = zero(T)
    @inbounds for (n, c) in s.pre
        lpre += c * lf[n+1]
    end
    s.sqrt_pre && (lpre /= 2)
    maxL = zero(T); acc = zero(T); started = false
    logtmax = zero(T); seen = false
    @inbounds for seg in segs
        z0 = first(seg)
        l0 = zero(T)
        for f in s.fac
            l0 += f.c * lf[_arg(f, z0)+1]
        end
        t = one(T); ssum = one(T); tmax = one(T); sc = 0
        for z in z0:last(seg)-1
            r = s.alternating ? -one(T) : one(T)
            for f in s.fac
                n = _arg(f, z)
                if f.a == 1                                     # [n+1]!/[n]! = [n+1]
                    r *= f.c == 1 ? q[n+1] : f.c == -1 ? qi[n+1] : q[n+1]^f.c
                elseif f.a == -1                                # [n−1]!/[n]! = 1/[n]
                    r *= f.c == 1 ? qi[n] : f.c == -1 ? q[n] : qi[n]^f.c
                end
            end
            t *= r
            ssum += t
            at = abs(t)
            at > tmax && (tmax = at)
            if at > big_
                t *= small_; ssum *= small_; tmax *= small_; sc += 100
            end
        end
        Lt = l0 + sc * ln10 + log(tmax)
        logtmax = seen ? max(logtmax, Lt) : Lt
        seen = true
        iszero(ssum) && continue
        L = l0 + sc * ln10 + log(abs(ssum))
        sgn = ((s.alternating && isodd(z0)) ? -one(T) : one(T)) * sign(ssum)
        if !started
            acc = sgn; maxL = L; started = true
        elseif L > maxL
            acc = acc * exp(maxL - L) + sgn
            maxL = L
        else
            acc += sgn * exp(L - maxL)
        end
    end
    (!started || iszero(acc)) && return zero(T), Inf
    value = s.sign0 * acc * exp(maxL + lpre)
    loss = Float64((logtmax - (maxL + log(abs(acc)))) / ln10)
    return value, loss
end

"Decimal digits a level value in type T must keep; below this the sum is redone at higher precision."
_target_digits(::Type{T}) where {T} = T === BigFloat ? _decimal_digits(T) - 4 : min(12.0, _decimal_digits(T) - 4)

"""
    value_at_level(s, k, T; fallback) -> T

Value at q = e^{iπ/(k+2)}. Poles throw a `DomainError` and vanishing terms are skipped, both decided from
valuations before any arithmetic. The sum runs in T and reports how many digits cancellation consumed.
If that leaves fewer than `_target_digits(T)`, the modular test first decides whether the value is exactly
zero; otherwise the sum is redone in BigFloat with enough bits for the loss (doubling until the digits
are there). `fallback()` is used if a factorial falls outside the level tables.
"""
function value_at_level(s::FactorialSum, k::Int, ::Type{T}; fallback) where {T}
    st, segs = classify_at_level(s, k)
    st === :pole && throw(DomainError(k, "Topological pole at level k=$k."))
    (st === :empty || st === :zero) && return zero(T)
    _within_level_tables(s, segs, k) || return fallback()
    v, loss = _sum_at_level(s, segs, k, T)
    target = _target_digits(T)
    _decimal_digits(T) - loss >= target && return v
    is_cancellation_zero(s, segs, k) === true && return zero(T)
    bits = 64 * cld(ceil(Int, (target + min(loss, 1e6) + 10) * log2(10)), 64)
    while true
        w, lossw = setprecision(BigFloat, bits) do
            _sum_at_level(s, segs, k, BigFloat)
        end
        bits * log10(2.0) - lossw >= target + 2 && return T(w)
        bits *= 2
    end
end
