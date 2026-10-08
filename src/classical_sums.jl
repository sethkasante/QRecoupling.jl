# ---------------------------------------------------------------------------------
#  The classical limit (q → 1) from the factorial rule
#
#  At q = 1, the same ratio loop runs over tables of n and split-exponent n!, 
#  with the same cancellation estimate and precision escalation.
#  Cancellation of zeros also exist. Modular residues exclude nonzeros cheaply.
#---------------------------------------------------------------------------------

"Ordinary integers, their inverses and split-exponent factorials up to N, in the layout of the level tables."
function ClassicalTables(::Type{T}, N::Int) where {T}
    N = max(N, 1)
    return _qint_tables_from(T, () -> BigFloat[BigFloat(n) for n in 1:N]; classical = true)
end

# Tables are sized in powers of two, cached by the exponent, so one table serves every symbol whose
# largest factorial fits.
const CLASSICAL_F64_TABLES = LevelCache{QIntTables{Float64}}()

_classical_bucket(N::Int) = N <= 64 ? 6 : 8 * sizeof(Int) - leading_zeros(N - 1)    # max(6, ⌈log₂ N⌉)

function classical_tables(::Type{Float64}, N::Int)
    b = _classical_bucket(N)
    return get_level!(() -> ClassicalTables(Float64, 1 << b), CLASSICAL_F64_TABLES, b)
end

const CLASSICAL_TABLES = Dict{Tuple{DataType,Int,Int},Any}()
const CLASSICAL_TABLES_LOCK = ReentrantLock()

function classical_tables(::Type{T}, N::Int) where {T}
    b = _classical_bucket(N)
    key = (T, b, T === BigFloat ? precision(BigFloat) : 0)
    tab = @lock CLASSICAL_TABLES_LOCK get!(() -> ClassicalTables(T, 1 << b), CLASSICAL_TABLES, key)
    return tab::QIntTables{T}
end

"Largest factorial argument the rule can reach."
function max_argument(s::FactorialSum)
    m = 0
    for (n, _) in s.pre
        m = max(m, Int(n))
    end
    for f in s.fac
        m = max(m, _arg(f, s.zlo), _arg(f, s.zhi))
    end
    return m
end

# ---- exact classical zeros ----

"n! and 1/n! modulo a prime p > N (Montgomery form)."
struct ClassicalModTable
    m::Montgomery
    fact::Vector{UInt64}
    invfact::Vector{UInt64}
end

function ClassicalModTable(p::UInt64, N::Int)
    m = Montgomery(p)
    fact = Vector{UInt64}(undef, N + 1)
    fact[1] = m.one
    for n in 1:N
        fact[n+1] = mont_mul(m, fact[n], to_mont(m, n))
    end
    invf = similar(fact)
    invf[N+1] = mont_inv(m, fact[N+1])
    for n in N:-1:1
        invf[n] = mont_mul(m, invf[n+1], to_mont(m, n))
    end
    return ClassicalModTable(m, fact, invf)
end

# Two fixed word-size primes filter out nonzeros. Vanishing residues still require exact confirmation.
const CLASSICAL_PRIMES = (UInt64(4611686018427387847), UInt64(4611686018427387817))
const CLASSICAL_MOD_TABLES = LevelCache{Tuple{ClassicalModTable,ClassicalModTable}}()

function classical_mod_tables(N::Int)
    b = _classical_bucket(N)
    return get_level!(CLASSICAL_MOD_TABLES, b) do
        (ClassicalModTable(CLASSICAL_PRIMES[1], 1 << b), ClassicalModTable(CLASSICAL_PRIMES[2], 1 << b))
    end
end

"Is the classical sum exactly zero? The prefactor is a nonzero square root, so only the sum matters."
function is_classical_zero(s::FactorialSum, tabs = nothing)
    is_empty_sum(s) && return true
    tabs === nothing && (tabs = classical_mod_tables(max_argument(s)))
    _sum_vanishes_mod(s, tabs) || return false
    P, _ = try
        _horner_sum(s)
    catch e
        e isa OverflowError || rethrow()
        _horner_sum_big(s)
    end
    return iszero(P)
end

"Does the rule's sum vanish modulo every table in `tabs`? Each table holds [n]! and 1/[n]! at one point."
function _sum_vanishes_mod(s::FactorialSum, tabs)
    for tab in tabs
        m = tab.m
        acc = UInt64(0)
        @inbounds for z in s.zlo:s.zhi
            t = m.one
            for f in s.fac
                n = _arg(f, z)
                b = f.c > 0 ? tab.fact[n+1] : tab.invfact[n+1]
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

"""
    pairwise_zero(s) -> Bool

Does the sum cancel term by term, for every q? The reflection `z ↦ zlo + zhi - z` maps the factor
`[a z + b]!^c` to `[-a z + (a(zlo+zhi) + b)]!^c`. When it maps the rule's multiset of factors onto itself,
the terms at `z` and at its mirror are the same product of q-factorials; when the sum alternates and
`zlo + zhi` is odd they have opposite signs and no term is its own mirror, so every pair cancels.

This is a proof, not a screen, and it holds at every q and every level (a term that vanishes at a level
has a mirror that vanishes with it). It is the column-exchange selection rule of the 3j symbol —
`(j j j; m₁ m₂ m₁)` with `3j` odd, `(j₁ j₂ j₃; 0 0 0)` with `j₁ + j₂ + j₃` odd — read off the rule rather
than enumerated, so it covers every symbol whose sum has that shape. A cancellation zero whose sum is not
pairwise can be screened modularly and then confirmed by exact arithmetic.
"""
function pairwise_zero(s::FactorialSum)
    (is_empty_sum(s) || !s.alternating) && return false
    c = s.zlo + s.zhi
    isodd(c) || return false
    @inbounds for f in s.fac
        a, b, e = Int(f.a), Int(f.b), Int(f.c)
        ma, mb = -a, a * c + b
        nf = 0; nm = 0
        for g in s.fac
            ga, gb, ge = Int(g.a), Int(g.b), Int(g.c)
            ge == e || continue
            (ga == a && gb == b) && (nf += 1)
            (ga == ma && gb == mb) && (nm += 1)
        end
        nf == nm || return false
    end
    return true
end

# ---- identically zero at generic q: modular filtering and polynomial confirmation ----

const _GENERIC_ZERO_SEEDS = (UInt64(0x1f3a9c2d7e4b5a61), UInt64(0x2b7e151628aed2a6))
const GENERIC_MOD_TABLES = LevelCache{Tuple{ClassicalModTable,ClassicalModTable}}()

"[n]! and 1/[n]! at a point r of F_p where no [n], n ≤ N, vanishes."
function _generic_mod_table(p::UInt64, N::Int, seed::UInt64)
    m = Montgomery(p)
    r0 = seed % p
    qv = Vector{UInt64}(undef, N)
    while true
        r = to_mont(m, r0)
        ri = iszero(r0) ? r : mont_inv(m, r)
        d = mont_sub(m, r, ri)
        good = !iszero(r0) && !iszero(d)
        fact = Vector{UInt64}(undef, N + 1)
        if good
            dinv = mont_inv(m, d)
            fact[1] = m.one
            rn = m.one; rin = m.one
            for n in 1:N
                rn = mont_mul(m, rn, r); rin = mont_mul(m, rin, ri)
                qv[n] = mont_mul(m, mont_sub(m, rn, rin), dinv)   # [n] = (r^n − r^−n)/(r − r^−1)
                iszero(qv[n]) && (good = false; break)
                fact[n+1] = mont_mul(m, fact[n], qv[n])
            end
        end
        if good
            invf = similar(fact)
            invf[N+1] = mont_inv(m, fact[N+1])
            for n in N:-1:1
                invf[n] = mont_mul(m, invf[n+1], qv[n])
            end
            return ClassicalModTable(m, fact, invf)
        end
        r0 = UInt64((UInt128(r0) * 6364136223846793005 + 1442695040888963407) % p)
    end
end

function generic_mod_tables(N::Int)
    b = _classical_bucket(N)
    return get_level!(GENERIC_MOD_TABLES, b) do
        (_generic_mod_table(CLASSICAL_PRIMES[1], 1 << b, _GENERIC_ZERO_SEEDS[1]),
         _generic_mod_table(CLASSICAL_PRIMES[2], 1 << b, _GENERIC_ZERO_SEEDS[2]))
    end
end

"Is the sum the zero function of q? Nonzero residues prove it is not; see the section comment above."
function is_generic_zero(s::FactorialSum, tabs = nothing)
    is_empty_sum(s) && return true
    pairwise_zero(s) && return true
    tabs === nothing && (tabs = generic_mod_tables(max_argument(s)))
    return _sum_vanishes_mod(s, tabs) && iszero(generic_value(s))
end

# ---- exact evaluation by Horner nesting (replaces BigFloat escalation for Float64 results) ----

using Base.GMP: MPZ

function _split_product(tab::QIntTables{Float64}, pairs)
    m = 1.0; e = 0
    for (n, c) in pairs
        m, e = _split_mul(m, e, tab, Int(n), c)
    end
    return _renorm(m, e)
end

"""
Upper bound on the bits the Horner accumulators reach, Σ_z log₂(|a_z| + b_z), so that they can be allocated once.
"""
function _horner_bits(s::FactorialSum)
    bits = 64
    for z in (s.zhi - 1):-1:s.zlo
        for f in s.fac
            lo, hi, c = _factor_step(f, z)
            lo > hi && continue
            bits += abs(c) * (hi - lo + 1) * (64 - leading_zeros(max(hi, 1)))
        end
        bits += 1
    end
    return bits
end

"Exact S / t_lo = P / Q by Horner nesting; throws OverflowError if a ratio exceeds a machine word."
function _horner_sum(s::FactorialSum)
    nb = _horner_bits(s)
    P = BigInt(; nbits = nb); Q = BigInt(; nbits = nb); T1 = BigInt(; nbits = nb); T2 = BigInt(; nbits = nb)
    MPZ.set_si!(P, 1); MPZ.set_si!(Q, 1)
    for z in (s.zhi - 1):-1:s.zlo
        a = s.alternating ? -1 : 1
        b = 1
        for f in s.fac
            lo,hi,c = _factor_step(f,z)
            if lo == hi && abs(c) == 1                  # the 6j and 3j case: one integer per factor
                c > 0 ? (a = Base.checked_mul(a, lo)) : (b = Base.checked_mul(b, lo))
                continue
            end
            for x in lo:hi, _ in 1:abs(c)
                if c > 0
                    a = Base.checked_mul(a,x)
                else
                    b = Base.checked_mul(b,x)
                end
            end
        end
        MPZ.mul_si!(T1, Q, b)
        MPZ.mul_si!(T2, P, a)
        MPZ.add!(P, T1, T2)
        MPZ.set!(Q, T1)
    end
    return P, Q
end

"The classical value in Float64 with no cancellation anywhere: exact sum, split-exponent prefactor."
function classical_exact_float(s::FactorialSum)
    is_empty_sum(s) && return 0.0
    PQ = try
        _horner_sum(s)
    catch e
        e isa OverflowError ? nothing : rethrow()
    end
    PQ === nothing && return nothing
    P, Q = PQ
    iszero(P) && return 0.0
    tab = classical_tables(Float64, max_argument(s))
    mp, ep = _split_product(tab, ((Int(n), Int(c)) for (n, c) in s.pre))
    if s.sqrt_pre
        isodd(ep) && (mp *= 2; ep -= 1)
        mp = sqrt(mp); ep ÷= 2
    end
    mt, et = _split_product(tab, ((_arg(f, s.zlo), Int(f.c)) for f in s.fac))
    sP = max(ndigits(P; base = 2) - 64, 0); sQ = max(ndigits(Q; base = 2) - 64, 0)
    mr, er = frexp(Float64(P >> sP) / Float64(Q >> sQ))
    sg = (s.alternating && isodd(s.zlo)) ? -1 : 1
    return s.sign0 * sg * ldexp(mp * mt * mr, ep + et + er + sP - sQ)
end

# ---- short sums with small arguments ----
#
# For small n, n! and 1/n! are `Float64` values far inside the exponent range, so a product of a few of them
# needs no separate exponent. `_small_classical` runs the operations of the plain pass on these values, with
# each exponent folded into its factor. Floating-point products, quotients and square roots do not depend on
# the scaling while nothing overflows or underflows, so the value is bit for bit that of `_sum_at_level`, at
# a fraction of the fixed cost. A budget of binary exponents keeps every partial product in range.

"Largest factorial argument of the short path."
const SMALL_NMAX = 64
"Exponent budget of one product: factors of one sign may add up to at most this many bits."
const SMALL_BITS = 900
"Largest term, relative to the first, that the short path keeps; the plain pass rescales far above it."
const SMALL_TERM_MAX = 2.0^64
"Results outside [2^-800, 2^800] go to the plain pass, which carries its exponent separately."
const SMALL_RANGE = 2.0^800

# n! and 1/n! as the plain pass sees them (mantissa · 2^exponent of the classical table), and the number of
# bits e with 2^-e < 1/n! and n! < 2^e.
const SMALL_FACT, SMALL_INVFACT, SMALL_EXP = let tab = ClassicalTables(Float64, SMALL_NMAX)
    ntuple(i -> ldexp(tab.fm[i], tab.fe[i]), SMALL_NMAX + 1), ntuple(i -> ldexp(tab.gm[i], tab.ge[i]), SMALL_NMAX + 1),
    ntuple(i -> max(tab.fe[i], 1 - tab.ge[i]), SMALL_NMAX + 1)
end

"""
    _small_classical(s, rtol) -> Float64 or nothing

The plain pass for a `Float64` classical sum whose factorial arguments are at most `SMALL_NMAX` and whose
products fit the exponent budget, accepted at `rtol`. `nothing` when the rule is outside that class or the
estimate does not pass; the caller then runs the general evaluation. A returned value is the one the plain
pass of the general evaluation computes.
"""
@inline function _small_classical(s::FactorialSum, rtol::Float64)
    T = Float64
    u = _unit(T)
    P = one(T); np = 0; bpos = 0; bneg = 0
    @inbounds for (n, c) in s.pre
        0 <= n <= SMALL_NMAX || return nothing
        if c == 1
            P *= SMALL_FACT[n+1]; bpos += SMALL_EXP[n+1]
        elseif c == -1
            P *= SMALL_INVFACT[n+1]; bneg += SMALL_EXP[n+1]
        else
            return nothing
        end
        np += 1
    end
    (bpos <= SMALL_BITS && bneg <= SMALL_BITS) || return nothing
    relpre = _gamma(T, 2np)
    if s.sqrt_pre
        P = sqrt(P)
        relpre = (relpre / 2) + u
    end
    z0 = s.zlo; z1 = s.zhi
    M0 = one(T); nf = 0; kn = 0; kd = 0; m = 1; bpos = 0; bneg = 0
    @inbounds for f in s.fac
        (abs(f.a) == 1 && abs(f.c) == 1) || return nothing
        n0 = _arg(f, z0); n1 = _arg(f, z1)
        (0 <= n0 <= SMALL_NMAX && 0 <= n1 <= SMALL_NMAX) || return nothing
        if f.c == 1
            M0 *= SMALL_FACT[n0+1]; bpos += SMALL_EXP[n0+1]
        else
            M0 *= SMALL_INVFACT[n0+1]; bneg += SMALL_EXP[n0+1]
        end
        (f.a == 1) == (f.c > 0) ? (kn += 1) : (kd += 1)
        m = max(m, n0 + 1, n1 + 1)
        nf += 1
    end
    (bpos <= SMALL_BITS && bneg <= SMALL_BITS) || return nothing
    pw = one(T)                                  # the test of `_integer_ratios`: exact integer quotients
    for _ in 1:max(kn, kd)
        pw *= m
    end
    pw < 2.0^53 || return nothing
    t = one(T); ssum = one(T); W = zero(T); SS = zero(T); j = 0
    for z in z0:z1-1
        j += 1
        a, b = _int_ratio(s, z)
        t *= T(a) / T(b)
        ssum += t
        at = abs(t)
        W = fma(T(j), at, W); SS += abs(ssum)
        at > SMALL_TERM_MAX && return nothing
    end
    iszero(ssum) && return nothing
    g = _gamma(T, 2 * (z1 - z0))
    relfirst = _gamma(T, 2nf)
    sgn = (s.alternating && isodd(z0)) ? -one(T) : one(T)
    x = sgn * ssum * M0
    E = abs(M0) * (2 * u * W / (1 - g)^2 + u / (1 - u) * SS) + abs(ssum * M0) * (relfirst + u)
    xp = x * P
    bound = (E * abs(P) * (1 + relpre) + abs(xp) * (relpre + u)) * BOUND_SLACK(T)
    v = s.sign0 * xp
    (inv(SMALL_RANGE) < abs(v) < SMALL_RANGE && _certifies(v, bound, rtol)) || return nothing
    return v
end

"""
    classical_value(s, T) -> T

The symbol described by rule `s` at q = 1. Cancellation is measured, exact zeros 
come back as zero, and when too few digits survive a `Float64` result is  
recomputed exactly (Horner nesting in integers); other types escalate in `BigFloat`.
"""
function classical_value(s::FactorialSum, ::Type{T}; labels = nothing, workspace=nothing) where {T}
    is_empty_sum(s) && return zero(T)
    if T === Float64
        mode = POLICY[]
        if mode !== :compensated_only        # that policy has no plain pass to reproduce
            vs = _small_classical(s, mode === :strict || mode === :strict_lazy ? RTOL_CERTIFIED : RTOL_PLAIN)
            vs === nothing || return vs
        end
    end
    N = max_argument(s)
    segs = (s.zlo:s.zhi,)
    v, st = _certified_value(s, segs, 0, classical_tables(T, N), workspace)
    st === :done && return v
    target = _target_digits(T)
    pess = _bound_pessimism(s, segs, classical_tables(Float64, N))
    if T === Float64 && labels !== nothing       # the symbol as one entry of its column: O(distance), no κ
        vr = sixj_entry(labels, ClassicalQ(), _column_workspace(workspace))
        vr === nothing || return T(vr)           # a returned value is certified far from zero
    end
    (pairwise_zero(s) || is_classical_zero(s)) && return zero(T)
    if T === Float64
        v = classical_exact_float(s)
        v === nothing || return v
    end
    bits = 64 * cld(ceil(Int, (target + 26) * log2(10)), 64)
    while true
        w, B = setprecision(BigFloat, bits) do
            _sum_at_level(s, segs, 0, classical_tables(BigFloat, N))
        end
        _surviving_digits(w, B, pess) >= target + 2 && return T(w)
        bits *= 2
    end
end

"Is `q` the classical point q = 1?"
_is_classical(q) = !isnothing(q) && (q == 1 || q == 1.0)
