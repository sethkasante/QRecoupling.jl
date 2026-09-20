# ---------------------------------------------------------------------------------
#  The classical limit q → 1 from the same factorial rule
#
#  At q = 1 every q-integer is an ordinary integer, so the rule that describes a symbol at a level also
#  describes it classically: the same ratio loop runs over tables of n and log n!, with the same
#  cancellation estimate and precision escalation. There are no valuations to track (Φ_h never divides
#  an ordinary factorial), but cancellation zeros exist classically too — the "non-trivial zeros" of
#  6j symbols — and are settled the same way, by evaluating the rational sum modulo large primes.
#
#  This replaces a floating-point projection of the expanded form that lost all accuracy at large spin
#  (12% off at j = 100 for {j j j; j j j}).
# ---------------------------------------------------------------------------------

"Ordinary integers, their inverses and split-exponent factorials up to N, in the layout of the level tables."
function ClassicalTables(::Type{T}, N::Int) where {T}
    N = max(N, 1)
    return _qint_tables_from(T, () -> BigFloat[BigFloat(n) for n in 1:N]; classical = true)
end

# Tables are sized in powers of two, cached by the exponent, so one table serves every symbol whose
# largest factorial fits.
const CLASSICAL_F64_TABLES = LevelCache{QIntTables{Float64}}()

_classical_bucket(N::Int) = max(6, ceil(Int, log2(max(N, 2))))

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

# Two independent word-size primes: a nonzero rational sum vanishes modulo both with probability ~2⁻¹²⁴.
const CLASSICAL_PRIMES = (UInt64(4611686018427387847), UInt64(4611686018427387817))
const CLASSICAL_MOD_TABLES = LevelCache{Tuple{ClassicalModTable,ClassicalModTable}}()

function classical_mod_tables(N::Int)
    b = _classical_bucket(N)
    return get_level!(CLASSICAL_MOD_TABLES, b) do
        (ClassicalModTable(CLASSICAL_PRIMES[1], 1 << b), ClassicalModTable(CLASSICAL_PRIMES[2], 1 << b))
    end
end

"Is the classical sum exactly zero? The prefactor is a nonzero square root, so only the sum matters."
function is_classical_zero(s::FactorialSum)
    is_empty_sum(s) && return true
    for tab in classical_mod_tables(max_argument(s))
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

# ---- exact evaluation by Horner nesting (replaces BigFloat escalation for Float64 results) ----
#
# Consecutive terms differ by a ratio of small integers, t_{z+1}/t_z = a_z/b_z, so the sum nests as
# t_lo (1 + r_lo (1 + r_{lo+1} (1 + ...))) and, carried as P/Q from the top, P <- b Q + a P, Q <- b Q.
# Every step is big integer x machine word: exact, linear cost, no prime tables. The prefactor and first
# term come from split-exponent factorial tables (correctly rounded mantissa and binary exponent), which
# keep ~20 eps of accuracy at any size where log tables lose |log n!| eps. One rounding at the end.
# Measured against WignerSymbols.jl's exact prime-factorised route: 4-13x faster, error <= 6e-16.

using Base.GMP: MPZ

function _split_product(tab::QIntTables{Float64}, pairs)
    m = 1.0; e = 0
    for (n, c) in pairs
        m, e = _split_mul(m, e, tab, Int(n), c)
    end
    return _renorm(m, e)
end

"Exact S / t_lo = P / Q by Horner nesting; `nothing` if a ratio would overflow a machine word."
function _horner_sum(s::FactorialSum)
    P = big(1); Q = big(1); T1 = big(0); T2 = big(0)
    for z in (s.zhi - 1):-1:s.zlo
        a = s.alternating ? -1 : 1
        b = 1
        for f in s.fac
            lo,hi,c = _factor_step(f,z)
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

"""
    classical_value(s, T) -> T

The symbol described by rule `s` at q = 1. Same contract as the level path: the sum runs in `T`,
cancellation is measured, exact zeros come back as zero, and when too few digits survive a `Float64`
result is recomputed exactly (Horner nesting in integers); other types escalate in `BigFloat`.
"""
function classical_value(s::FactorialSum, ::Type{T}; labels = nothing, workspace=nothing) where {T}
    is_empty_sum(s) && return zero(T)
    N = max_argument(s)
    segs = (s.zlo:s.zhi,)
    v, st = _certified_value(s, segs, 0, classical_tables(T, N), workspace)
    st === :done && return v
    target = _target_digits(T)
    _, loss = _sum_at_level(s, segs, 0, classical_tables(T, N))
    if T === Float64 && labels !== nothing       # the symbol as one entry of its column: O(distance), no κ
        vr = sixj_entry(labels, ClassicalQ(), _column_workspace(workspace))
        vr === nothing || return T(vr)           # a returned value is certified far from zero
    end
    is_classical_zero(s) && return zero(T)
    if T === Float64
        v = classical_exact_float(s)
        v === nothing || return v
    end
    bits = 64 * cld(ceil(Int, (target + min(loss, 1e6) + 10) * log2(10)), 64)
    while true
        w, lossw = setprecision(BigFloat, bits) do
            _sum_at_level(s, segs, 0, classical_tables(BigFloat, N))
        end
        bits * log10(2.0) - lossw >= target + 2 && return T(w)
        bits *= 2
    end
end

"Is `q` the classical point q = 1?"
_is_classical(q) = !isnothing(q) && (q == 1 || q == 1.0)
