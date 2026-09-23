# ---------------------------------------------------------------------------------
#  Per-level tables: exact zero tests, q-integer and split-factorial tables, ratio helpers 
#
#  Everything a level needs before any symbol is evaluated. Tables are built once per level and read
#  without a lock (`LevelCache`), which is what makes batches at one level cheap: the tables are the
#  expensive part and every symbol at the level shares them.
# ---------------------------------------------------------------------------------

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
        dinv = mont_inv(m, mont_sub(m, r, ri))   # 1/(r − r^{−1}
        rn = m.one; rin = m.one
        fact[1, c] = m.one
        for n in 1:N
            rn = mont_mul(m, rn, r); rin = mont_mul(m, rin, ri)
            qv[n] = mont_mul(m, mont_sub(m, rn, rin), dinv)    # [n] = (r^n − r^{−n}) / (r − r^{−1})
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
    k >= 0 || throw(DomainError(k, "level/cache index must be nonnegative"))
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
    qq::Vector{T}     # [q; qi] and [ql; qil] in one array, so a ratio factor is one load at an affine index
    qql::Vector{T}
end

QIntTables{T}(q, qi, ql, qil, fm, fe, fml, gm, ge, gml, classical::Bool) where {T} =
    QIntTables{T}(q, qi, ql, qil, fm, fe, fml, gm, ge, gml, classical, vcat(q, qi), vcat(ql, qil))

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
    @inbounds if c == 1                          # the usual exponents, without a loop
        return m * tab.fm[n+1], e + tab.fe[n+1]
    elseif c == -1
        return m * tab.gm[n+1], e + tab.ge[n+1]
    elseif c > 0
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
        (abs(f.a) == 1 && abs(f.c) == 1) || return false
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

"""
    _ratio_plan(s, N) -> (ok, plan)

The ratio t_{z+1}/t_z as one load per factor. For [a z + b]!^{±1} with a = ±1 the factor is [n+1]^{±1}
(a = 1) or [n]^{∓1} (a = −1), n = a z + b, i.e. entry `qq[a z + o]` of the combined table `[q; qi]` (length
2N) for a fixed offset o. `plan` holds the pairs (a, o); `ok` is false when some factor has another slope
or exponent, and the kernels then use the general branchy loop. Same entries in the same order, so the
result is bitwise the same — the gain is that the loop has no data-dependent branches (1.3 ns instead of
~10 ns per step for a 6j).
"""
@inline function _ratio_plan(s::FactorialSum, N::Int)
    ok = all(f -> (f.a == 1 || f.a == -1) && (f.c == 1 || f.c == -1), s.fac)
    plan = map(f -> _plan_entry(f, N), s.fac)
    return ok, plan
end

@inline function _plan_entry(f::AffineFactorial, N::Int)
    b = Int(f.b)
    return f.a == 1 ? (1, f.c > 0 ? b + 1 : b + 1 + N) : (-1, f.c > 0 ? b + N : b)
end

"Adjacent-term ratio for arbitrary affine slopes and integer powers."
function _general_ratio(s::FactorialSum, z::Int, tab::QIntTables{T}) where {T}
    r = s.alternating ? -one(T) : one(T)
    for f in s.fac
        lo, hi, c = _factor_step(f,z)
        c == 0 && continue
        q = c > 0 ? tab.q : tab.qi
        @inbounds for n in lo:hi, _ in 1:abs(c)
            r *= q[n]
        end
    end
    return r
end

"Double-word adjacent-term ratio; the standard unit-slope plan remains the fast path."
function _general_ratio_dw(s::FactorialSum, z::Int, tab::QIntTables{T}) where {T}
    rh = s.alternating ? -one(T) : one(T)
    rl = zero(T)
    for f in s.fac
        lo, hi, c = _factor_step(f,z)
        c == 0 && continue
        qh, ql = c > 0 ? (tab.q,tab.ql) : (tab.qi,tab.qil)
        @inbounds for n in lo:hi, _ in 1:abs(c)
            xh, xl = qh[n], ql[n]
            p, e = _two_prod(rh,xh)
            rl = fma(rl,xh,fma(rh,xl,e))
            rh = p
        end
    end
    return rh, rl
end

"Unit roundoff of T."
_unit(::Type{T}) where {T} = eps(T) / 2

"γ_m = m·u/(1 − m·u), the standard bound on m relative roundings (Inf when m·u ≥ 1/2)."
@inline function _gamma(::Type{T}, m::Integer) where {T}
    mu = m * _unit(T)
    return mu < 0.5 ? T(mu / (1 - mu)) : T(Inf)
end
