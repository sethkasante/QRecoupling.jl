
# ---------------------------------------------------------------------------------
#  Per-level tables: modular zero screening, q-integer and split-factorial tables, ratio helpers
#
#  Built once per level and read without a lock (`LevelCache`); every symbol at the level shares them.
# ---------------------------------------------------------------------------------

# ---- modular zero screen: two Galois conjugates modulo one prime ----

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
#
# A table grows with its level, so a sweep over many levels (one symbol at k = 1, 2, …, or spins and level
# scaled together) would otherwise keep every table it has built: tens of gigabytes by k ~ 10⁴. Each cache
# therefore holds at most `LEVEL_CACHE_LIMIT[]` bytes. When a new table would pass that, the cache starts
# again from this table alone; callers that still hold an older table keep it until they are done.

"Bytes one level cache may hold before it is emptied to make room (set before computing; default 512 MiB)."
const LEVEL_CACHE_LIMIT = Ref(512 * 2^20)

mutable struct LevelCache{V}
    @atomic slots::Vector{Union{Nothing,V}}
    lock::ReentrantLock
    bytes::Int                      # size of the cached tables, kept under the lock
end

LevelCache{V}() where {V} = LevelCache{V}(Vector{Union{Nothing,V}}(nothing, 0), ReentrantLock(), 0)

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
        size = Base.summarysize(v)
        fresh = if c.bytes + size > LEVEL_CACHE_LIMIT[] && c.bytes > 0
            c.bytes = 0                               # over the limit: keep only the new table
            Vector{Union{Nothing,V}}(nothing, max(k + 1, 16))
        elseif k + 1 > length(sl)                     # grow only when the level is out of range
            g = Vector{Union{Nothing,V}}(nothing, max(k + 1, 2 * length(sl), 16))
            copyto!(g, 1, sl, 1, length(sl))
            g
        else
            copy(sl)                                  # never mutate a vector a reader may hold
        end
        fresh[k+1] = v
        c.bytes += size
        @atomic c.slots = fresh
        return v
    end
end

function Base.empty!(c::LevelCache{V}) where {V}
    @lock c.lock begin
        @atomic c.slots = Vector{Union{Nothing,V}}(nothing, 0)
        c.bytes = 0
    end
    return c
end

const LEVEL_ZERO_TABLES = LevelCache{LevelZeroTable}()

level_zero_table(k::Int) = get_level!(() -> LevelZeroTable(k), LEVEL_ZERO_TABLES, k)


"""
    is_cancellation_zero(s, segs, k) -> Bool or nothing

Screen the contributing sum at two Galois conjugates modulo one prime. `false` excludes a zero;
`true` is a candidate requiring exact confirmation, without a universal false-positive probability.
Return `nothing` if a factorial falls outside the level tables.
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

# ---- an exact zero in the integer group ring ----
#
# With t_{z+1}/t_z = ±N_z/D_z, products of q-integers, a sum over one segment is its first term times
# A / Π_y D_y, where A = Σ_z (±ζ^w)^(z-z0) (Π_{y<z} N_y)(Π_{y≥z} D_y). When numerator and denominator have
# equally many factors, each [x] can be replaced by ζ^x − ζ^(−x): the common power of ζ − ζ^(−1) drops out
# of the question whether A vanishes. A is then an integer combination of powers of ζ = e^{iπ/h}, held as a
# vector of 2h integers, and multiplying by ζ^x − ζ^(−x) is two shifts and a subtraction. The nesting
# G_z = L_z ± ζ^w N_z G_{z+1}, L_z = D_z L_{z+1}, gives A = G_{z0}, and A = 0 exactly when the polynomial
# with these coefficients is divisible by the cyclotomic polynomial Φ_{2h}.

const LEVEL_PHI = LevelCache{Vector{Int}}()
"Coefficients of Φ_{2h}, h = k + 2, from the constant term up."
_level_phi(k::Int) = get_level!(LEVEL_PHI, k) do
    p = phi_q(2 * (k + 2))
    Int[Int(coeff(p, i)) for i in 0:degree(p)]
end

"v ← (ζ^x − ζ^(−x)) v on the cyclic group of order length(v), through the scratch vector `tmp`."
@inline function _binomial_mul!(v::Vector{Int}, tmp::Vector{Int}, x::Int)
    n = length(v)
    x = mod(x, n)
    @inbounds for j in 1:n
        a = j - x; a < 1 && (a += n)
        b = j + x; b > n && (b -= n)
        tmp[j] = v[a] - v[b]
    end
    copyto!(v, tmp)
    return v
end

"""
    _level_integer_zero(s, z0, z1, k, w = 0) -> Bool or nothing

Whether Σ_{z=z0}^{z1} (±1)^z ζ^{wz} t_z vanishes at ζ = e^{iπ/(k+2)}, decided in integers. The rule needs
unit slopes, exponents ±1, as many numerator as denominator factors in its term ratio, nonzero denominators,
and few enough steps for the coefficients to stay inside 64 bits; otherwise `nothing`. Both answers are
proofs.
"""
function _level_integer_zero(s::FactorialSum, z0::Int, z1::Int, k::Int, w::Int = 0)
    h = k + 2; n = 2h
    nsteps = z1 - z0
    nsteps >= 0 || return nothing
    kn = 0; kd = 0
    for f in s.fac
        (abs(f.a) == 1 && abs(f.c) == 1) || return nothing
        (f.a == 1) == (f.c > 0) ? (kn += 1) : (kd += 1)
    end
    kn == kd || return nothing
    # each factor at most doubles the coefficients: |G_z| ≤ (z1 − z + 1)·2^(kn (z1 − z))
    kn * nsteps + (8 * sizeof(Int) - leading_zeros(nsteps + 1)) <= 61 || return nothing
    G = zeros(Int, n); L = zeros(Int, n); tmp = Vector{Int}(undef, n)
    G[1] = 1; L[1] = 1
    wshift = mod(w, n)
    for z in z1-1:-1:z0
        for f in s.fac
            x = _arg(f, z); f.a == 1 && (x += 1)
            if (f.a == 1) == (f.c > 0)
                _binomial_mul!(G, tmp, x)
            else
                0 < x < h || return nothing           # a vanishing denominator: not a regular segment
                _binomial_mul!(L, tmp, x)
            end
        end
        # G ← L ± ζ^w G
        @inbounds for j in 1:n
            a = j - wshift; a < 1 && (a += n)
            tmp[j] = s.alternating ? L[j] - G[a] : L[j] + G[a]
        end
        copyto!(G, tmp)
    end
    # ζ^h = −1 folds the coefficients to degree below h; then the remainder modulo Φ_{2h}
    @inbounds for j in 1:h
        G[j] -= G[j+h]
    end
    phi = _level_phi(k); d = length(phi) - 1
    @inbounds for i in h:-1:d+1
        c = G[i]
        iszero(c) && continue
        for t in 0:d
            p, o1 = Base.Checked.mul_with_overflow(c, phi[t+1])
            r, o2 = Base.Checked.sub_with_overflow(G[i-d+t], p)
            (o1 | o2) && return nothing
            G[i-d+t] = r
        end
    end
    @inbounds for i in 1:d
        iszero(G[i]) || return false
    end
    return true
end

"Is the value exactly zero at level k (empty, all terms vanishing, or cancelling)?"
is_zero_at_level(s::FactorialSum, k::Int) = is_zero_at_level(s, k, level_zero_table(k))

function is_zero_at_level(s::FactorialSum, k::Int, tab::LevelZeroTable)
    result = _level_zero_screen(s, k, tab)
    result === nothing || return result
    _, segs = classify_at_level(s, k)
    return _level_exact_zero(s, segs, k)
end

"Proved zero/nonzero, or `nothing` when exact confirmation is needed. No exact-field work here."
function _level_zero_screen(s::FactorialSum, k::Int, tab::LevelZeroTable)
    st, segs = classify_at_level(s, k)
    (st === :empty || st === :zero) && return true
    st === :pole && return false
    pairwise_zero(s) && return true
    reflection_zero(s, segs, k) && return true          # a proof, and cheaper than the screen
    return is_cancellation_zero(s, segs, k, tab) === false ? false : nothing
end

"""
    reflection_zero(s, segs, k) -> Bool

Does the sum vanish at level `k` by a level reflection? A proof, not a screen. At q = e^{iπ/h},
[m]!·[h−1−m]! = [h−1]!, which gives each term a canonical form. When T(z) and T(c − z) agree, with c = z₀ + z₁
odd, z ↦ c − z pairs the terms with opposite signs. For the 6j this is {β₁, β₂, β₃, k} = c − {α₁, …, α₄};
it accounts for 71% of the zeros among all 6j symbols with k ≤ 22.
"""
function reflection_zero(s::FactorialSum, segs, k::Int)
    (s.alternating && length(segs) == 1) || return false
    z0 = first(segs[1]); z1 = last(segs[1])
    c = z0 + z1
    isodd(c) || return false
    h = k + 2
    κo = 0; κi = 0
    @inbounds for f in s.fac
        a = Int(f.a); b = Int(f.b)
        abs(a) == 1 || return false
        lo, hi = minmax(a * z0 + b, a * z1 + b)
        (lo >= 0 && hi <= h - 1) || return false            # where the reflection identity holds
        a == 1 ? (κi += Int(f.c)) : (κo += Int(f.c))
    end
    κo == κi || return false
    # canonical offsets: original  a = 1: (b, e);        a = −1: (h−1−b, −e)
    #                    reflected a = 1: (h−1−c−b, −e);  a = −1: (b−c, e)
    @inbounds for g in s.fac
        for which in 1:2
            x = which == 1 ? (Int(g.a) == 1 ? Int(g.b) : h - 1 - Int(g.b)) :
                             (Int(g.a) == 1 ? h - 1 - c - Int(g.b) : Int(g.b) - c)
            mo = 0; mi = 0
            for f in s.fac
                a = Int(f.a); b = Int(f.b); e = Int(f.c)
                (a == 1 ? b : h - 1 - b) == x && (mo += a == 1 ? e : -e)
                (a == 1 ? h - 1 - c - b : b - c) == x && (mi += a == 1 ? -e : e)
            end
            mo == mi || return false
        end
    end
    return true
end
reflection_zero(s::FactorialSum, k::Int) = (st_segs = classify_at_level(s, k);
    st_segs[1] === :finite && reflection_zero(s, st_segs[2], k))


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
    ff::Vector{T}     # fm · 2^fe and gm · 2^ge as plain numbers, for products short enough to need no
    gf::Vector{T}     # separate exponent (`_small_pass`); `Float64` tables only
end

function QIntTables{T}(q, qi, ql, qil, fm, fe, fml, gm, ge, gml, classical::Bool) where {T}
    ff, gf = T === Float64 ? (ldexp.(fm, fe), ldexp.(gm, ge)) : (T[], T[])
    return QIntTables{T}(q, qi, ql, qil, fm, fe, fml, gm, ge, gml, classical, vcat(q, qi), vcat(ql, qil), ff, gf)
end

"""
    split_factorials(T, qhi) -> (fm, fe)

Split-exponent factorials from q-integers in a wider `BigFloat` precision: the running product stays in
`BigFloat` and each [n]! is rounded to `T` once, so products keep ~20 ε at any size.
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
    lo(x) = T(x - T(x))     # the wide value minus its rounding, then rounded
    return QIntTables{T}(hi.(qw), hi.(qiw), lo.(qw), lo.(qiw),
                         hi.(ms), fe, lo.(ms), hi.(gs), ge, lo.(gs), classical)
end

"[1], …, [N] at q = e^{iπ/h} in the current BigFloat precision, by the Chebyshev recurrence
[n+1] = 2cos(π/h)[n] - [n-1] (errors grow linearly, absorbed by the 64 guard bits)."
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
const QINT_TABLES_BYTES = Ref(0)            # size of the tables in QINT_TABLES, kept under its lock
const QINT_TABLES_LOCK = ReentrantLock()

qint_tables(::Type{Float64}, k::Int) = get_level!(() -> QIntTables(Float64, k), QINT_F64_TABLES, k)

function qint_tables(::Type{T}, k::Int) where {T}
    key = (T, k, T === BigFloat ? precision(BigFloat) : 0)
    tab = @lock QINT_TABLES_LOCK begin
        t = get(QINT_TABLES, key, nothing)
        if t === nothing
            t = QIntTables(T, k)
            size = Base.summarysize(t)
            # the same limit as the level caches: wide tables are large, and one is built per level and precision
            if QINT_TABLES_BYTES[] + size > LEVEL_CACHE_LIMIT[]
                empty!(QINT_TABLES); QINT_TABLES_BYTES[] = 0
            end
            QINT_TABLES[key] = t; QINT_TABLES_BYTES[] += size
        end
        t
    end
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

Classically each term ratio is a/b with integers a, b. Below 2^53 both are exact in floating point, so a
ratio costs one division instead of K table multiplications, and its low part is exact from one fma.
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
    # m^max(kn, kd) < 2^53. Products of integers below 2^53 are exact, and one that reaches 2^53 stays there.
    x = 1.0
    for _ in 1:max(kn, kd)
        x *= m
    end
    return x < 2.0^53
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

The ratio t_{z+1}/t_z as one load per factor. For [a z + b]!^{±1} with a = ±1 the factor is entry
`qq[a z + o]` of the combined table `[q; qi]`; `plan` holds the pairs (a, o), and `ok` is false for any other
slope or exponent (the general loop then runs). Bitwise the same result, without data-dependent branches.
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

"γ_m = m·u/(1 - m·u), the standard bound on m relative roundings (Inf when m·u ≥ 1/2)."
@inline function _gamma(::Type{T}, m::Integer) where {T}
    mu = m * _unit(T)
    return mu < 0.5 ? T(mu / (1 - mu)) : T(Inf)
end
