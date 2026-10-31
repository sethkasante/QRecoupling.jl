# Exact level values in quantum integers, in plain Julia.
#
# At h = k + 2 the field ℚ(2cos(π/h)) has degree d = φ(2h)/2. The q-integers [1], …, [d] are a basis of it,
# because [n] is monic of degree n − 1 in x, and x[n] = [n+1] + [n−1] keeps a value in this basis. A rule is
# evaluated at the d conjugates q ↦ q^a modulo word primes p ≡ 1 (mod 2h), which split completely, and the
# coordinates are recovered by interpolation, CRT and rational reconstruction.
#
# The number of primes comes from a proof, not from the digits settling. Clear denominators, so that
# E = D·J·C·(value − w) is an algebraic integer, with J, C and D integers. E vanishes at every place above every
# prime used, so their product M divides E. A nonzero algebraic integer has some conjugate of size at least 1, so
# |σE| < M for every σ forces E = 0, and σE is bounded term by term from the sizes of the conjugates, computed in
# floating point with a margin.
#
# C can be an integer because a q-integer is nearly a unit: [n] = unit · (1 − ζ_o)/(1 − ζ_h) with o = h/gcd(n, h),
# and 1 − ζ_o is a unit unless o is a power of a prime ℓ, where its valuation is 1/φ(o) with ℓ normalized to 1.

"Tables at the d places above one prime."
struct LevelPlaces
    m::Montgomery
    as::Vector{Int}             # the conjugates q ↦ q^a
    pow::Vector{UInt64}         # powers of the root of unity
    den::Vector{UInt64}         # 1/(q^a − q^(−a)) at each place
    fact::Matrix{UInt64}        # row n + 1: [n]! with each [jh] replaced by its limit [jh]/[h] = (−1)^(j−1) j
    invf::Matrix{UInt64}
    X::Vector{UInt64}           # x at each place
    nd::Vector{Vector{UInt64}}  # inverse differences of X, for interpolation
end

"[n] at place `i`, or the limit [n]/[h] when h divides n."
@inline function _place_label(t::LevelPlaces, n::Int, i::Int, h::Int)
    if rem(n, h) == 0
        j = n ÷ h
        return to_mont(t.m, iseven(j) ? -j : j)
    end
    an = mod(t.as[i] * n, 2h)
    return mont_mul(t.m, mont_sub(t.m, t.pow[an+1], t.pow[mod(-an, 2h)+1]), t.den[i])
end

function _fill_factorials!(fact, invf, t::LevelPlaces, h::Int)
    m = t.m; N = size(fact, 1) - 1
    @inbounds for i in axes(fact, 2)
        fact[1, i] = m.one
        for n in 1:N
            fact[n+1, i] = mont_mul(m, fact[n, i], _place_label(t, n, i, h))
        end
        invf[N+1, i] = mont_inv(m, fact[N+1, i])
        for n in N:-1:1
            invf[n, i] = mont_mul(m, invf[n+1, i], _place_label(t, n, i, h))
        end
    end
end

function LevelPlaces(p::UInt64, h::Int, N::Int)
    m = Montgomery(p); r = to_mont(m, root_of_unity(p, 2h))
    pow = Vector{UInt64}(undef, 2h); pow[1] = m.one
    for j in 2:2h
        pow[j] = mont_mul(m, pow[j-1], r)
    end
    as = [a for a in 1:h-1 if gcd(a, 2h) == 1]; d = length(as)
    X = [mont_add(m, pow[a+1], pow[2h-a+1]) for a in as]
    den = [mont_inv(m, mont_sub(m, pow[a+1], pow[2h-a+1])) for a in as]
    nd = [[mont_inv(m, mont_sub(m, X[i], X[i-j+1])) for i in j:d] for j in 2:d]
    fact = Matrix{UInt64}(undef, N + 1, d); invf = similar(fact)
    t = LevelPlaces(m, as, pow, den, fact, invf, X, nd)
    _fill_factorials!(fact, invf, t, h)
    return t
end

function _grow(t::LevelPlaces, h::Int, N::Int)
    fact = Matrix{UInt64}(undef, N + 1, length(t.as)); invf = similar(fact)
    g = LevelPlaces(t.m, t.as, t.pow, t.den, fact, invf, t.X, t.nd)
    _fill_factorials!(fact, invf, g, h)
    return g
end

"The primes in use at one level and their tables."
mutable struct LevelField
    h::Int
    d::Int
    places::Vector{LevelPlaces}
    cursor::UInt64
    logs::Matrix{Float64}       # row n + 1: log2 of the size of [n]!, as in `LevelPlaces`, at each conjugate
    valuation::Vector{Tuple{Int,Float64}}   # for each [n], n ≤ h/2: the prime ℓ dividing it and its valuation
end
LevelField(h::Int) = LevelField(h, euler_phi(2h) ÷ 2, LevelPlaces[], ((UInt64(1) << 62) - 1) ÷ UInt64(2h),
                                Matrix{Float64}(undef, 0, 0), _label_valuations(h))

"For each [n] with n ≤ h/2: `(ℓ, v)` when it is not a unit, with ℓ normalized to valuation 1, and `(0, 0.0)` when it is."
function _label_valuations(h::Int)
    out = fill((0, 0.0), h ÷ 2)
    power_of(o) = (fs = _prime_factors(o); length(fs) == 1 ? fs[1] : 0)
    top = h > 1 ? power_of(h) : 0
    for n in 2:h÷2
        o = h ÷ gcd(n, h)
        l = o > 1 ? power_of(o) : 0
        l == 0 && continue
        v = 1 / euler_phi(o) - (top == 0 ? 0.0 : 1 / euler_phi(h))
        v > 0 && (out[n] = (l, v))
    end
    return out
end

const LEVEL_FIELDS = LRU{Int,LevelField}(maxsize = 64)
const LEVEL_FIELD_LOCK = ReentrantLock()

_level_field(h::Int) = get!(() -> LevelField(h), LEVEL_FIELDS, h)

"Tables of the i-th prime at this level, with factorials up to at least N."
function _places(F::LevelField, i::Int, N::Int)
    lock(LEVEL_FIELD_LOCK) do
        while length(F.places) < i
            p = 1 + UInt64(2F.h) * F.cursor
            F.cursor -= 1
            is_prime_u64(p) && push!(F.places, LevelPlaces(p, F.h, max(N, 2F.h)))
        end
        t = F.places[i]
        if size(t.fact, 1) <= N
            t = _grow(t, F.h, max(N, 2 * (size(t.fact, 1) - 1)))
            F.places[i] = t
        end
        return t
    end
end

"Sizes of the factorials at the conjugates of this level, up to at least N."
function _place_logs(F::LevelField, N::Int)
    lock(LEVEL_FIELD_LOCK) do
        size(F.logs, 1) > N && return F.logs
        h = F.h; N = max(N, 2h)
        as = [a for a in 1:h-1 if gcd(a, 2h) == 1]
        logs = zeros(Float64, N + 1, length(as))
        for (i, a) in enumerate(as), n in 1:N
            logs[n+1, i] = logs[n, i] + (rem(n, h) == 0 ? log2(n ÷ h) :
                                         log2(abs(sinpi(a * n / h) / sinpi(a / h))))
        end
        F.logs = logs
        return logs
    end
end

"What is needed to evaluate a rule at the places of a level: its factorial arguments, and the labels in front."
struct LevelSum
    args::Matrix{Int}               # one column per term
    powers::Vector{Int}
    negative::Vector{Bool}
    front::Vector{Pair{Int,Int}}    # ∏[n]^e multiplying the sum
    top::Int
end

"""
Bits of J·C·front·S at its largest conjugate, of J and of C, each rounded up. C clears the q-integers below the
line: for each prime ℓ, the smallest power of ℓ whose valuation covers theirs. The factorial table already holds
the rational in each term, so J appears once.
"""
function _level_sum_bound(L::LevelSum, front::Vector{Int}, es::Vector{Vector{Int}}, cs, h::Int)
    F = _level_field(h); logs = _place_logs(F, L.top)
    m = zeros(Int, length(front))
    for e in es, n in eachindex(m)
        m[n] = max(m[n], -(e[n] + front[n]))
    end
    J = big(1)
    for c in cs
        J = lcm(J, denominator(c))
    end
    need = Dict{Int,Float64}()
    for n in 2:length(m)
        l, v = F.valuation[n]
        l != 0 && m[n] > 0 && (need[l] = get(need, l, 0.0) + m[n] * v)
    end
    Cbits = 0.0
    for (l, v) in need
        Cbits += ceil(v + 1e-6) * log2(l)
    end
    top = -Inf; scale = 1.0
    @inbounds for a in axes(logs, 2)
        extra = 0.0
        for n in 2:length(front)
            iszero(front[n]) || (extra += front[n] * (logs[n+1, a] - logs[n, a]))
        end
        most = -Inf
        for z in axes(L.args, 2)
            b = 0.0
            for j in axes(L.args, 1)
                b += L.powers[j] * logs[L.args[j, z]+1, a]
            end
            most = max(most, b)
        end
        scale = max(scale, abs(most), abs(extra))
        top = max(top, most + extra)
    end
    margin = 8 + 1e-9 * scale
    Jbits = ndigits(J; base = 2); C = ceil(Int, Cbits + 1)
    return max(ceil(Int, top + log2(size(L.args, 2)) + margin), 0) + Jbits + C, Jbits, C
end

function _level_sum(s::FactorialSum, zs::Vector{Int}, front::Vector{Int}, h::Int)
    fs = collect(s.fac)
    args = [_arg(f, z) for f in fs, z in zs]
    return LevelSum(args, [Int(f.c) for f in fs], [s.alternating && isodd(z) for z in zs],
                    _nonzero_pairs(front), max(maximum(args; init = 1), h ÷ 2))
end

"The value at place `i`, in Montgomery form."
function _place_value(t::LevelPlaces, L::LevelSum, i::Int, h::Int)
    m = t.m; acc = zero(UInt64)
    @inbounds for z in axes(L.args, 2)
        u = m.one
        for j in axes(L.args, 1)
            c = L.powers[j]
            b = c > 0 ? t.fact[L.args[j, z]+1, i] : t.invf[L.args[j, z]+1, i]
            for _ in 1:abs(c)
                u = mont_mul(m, u, b)
            end
        end
        acc = L.negative[z] ? mont_sub(m, acc, u) : mont_add(m, acc, u)
    end
    iszero(acc) && return acc
    for (n, e) in L.front
        b = _place_label(t, n, i, h)
        e < 0 && (b = mont_inv(m, b))
        acc = mont_mul(m, acc, mont_pow(m, b, abs(e)))
    end
    return acc
end

"Coordinates on [1], …, [d] of the polynomial taking the values `c` at the places, out of Montgomery form."
function _place_coordinates!(c::Vector{UInt64}, p::Vector{UInt64}, t::LevelPlaces)
    m = t.m; d = length(c)
    @inbounds for j in 2:d
        w = t.nd[j-1]
        for i in d:-1:j
            c[i] = mont_mul(m, mont_sub(m, c[i], c[i-1]), w[i-j+1])
        end
    end
    # Newton's form, from the top, with (x − X)[n] = [n+1] + [n−1] − X[n]
    fill!(p, zero(UInt64)); p[1] = c[d]
    @inbounds for i in d-1:-1:1
        X = t.X[i]; below = zero(UInt64)
        for n in 1:d-i+1
            here = p[n]
            v = mont_add(m, below, n < d ? p[n+1] : zero(UInt64))
            p[n] = mont_sub(m, v, mont_mul(m, X, here))
            below = here
        end
        p[1] = mont_add(m, p[1], c[i])
    end
    @inbounds for n in 1:d
        p[n] = from_mont(m, p[n])
    end
    return p
end

"""
Coordinates over a common denominator from their residues modulo M, or `nothing` while M is too small. Whatever
this returns is only a candidate: the bound in `_level_coordinates` decides.
"""
function _rational_coordinates(acc::Vector{BigInt}, M::BigInt)
    D = big(1); half = ndigits(M; base = 2) ÷ 2
    out = similar(acc)
    for (i, a) in enumerate(acc)
        u = mod(a * D, M)
        2u > M && (u -= M)
        if ndigits(u; base = 2) + ndigits(D; base = 2) >= half
            r0 = M; r1 = mod(a * D, M); s0 = big(0); s1 = big(1)
            while ndigits(r1; base = 2) >= half
                q = r0 ÷ r1
                r0, r1 = r1, r0 - q * r1
                s0, s1 = s1, s0 - q * s1
            end
            (iszero(s1) || !isone(gcd(r1, s1))) && return nothing
            D *= abs(s1)
        end
    end
    for (i, a) in enumerate(acc)
        u = mod(a * D, M)
        2u > M && (u -= M)
        ndigits(u; base = 2) + ndigits(D; base = 2) >= half && return nothing
        out[i] = u
    end
    return out, D
end

"""
The coordinates of the cleared front·S on [1], …, [d], over a common denominator: `(numerators, D)`, proved
by the bound described at the head of this file.
"""
function _level_coordinates(L::LevelSum, bits::Int, Jbits::Int, Cbits::Int, h::Int)
    F = _level_field(h); d = F.d
    acc = [big(0) for _ in 1:d]; M = big(1)
    c = Vector{UInt64}(undef, d); p = Vector{UInt64}(undef, d)
    S = inv(sinpi(1 / h)); need = bits
    i = 0
    while true
        i += 1
        t = _places(F, i, L.top)
        for a in 1:d
            c[a] = _place_value(t, L, a, h)
        end
        _place_coordinates!(c, p, t)
        q = big(t.m.p); Minv = invmod(mod(M, q), q)
        for n in 1:d
            acc[n] += M * mod((big(p[n]) - mod(acc[n], q)) * Minv, q)
        end
        M *= q
        ndigits(M; base = 2) > need || continue
        rec = _rational_coordinates(acc, M)
        rec === nothing && continue
        w, D = rec
        isone(gcd(D, M)) || continue
        wbits = 0.0
        for n in 1:d
            iszero(w[n]) && continue
            wbits = max(wbits, ndigits(w[n]; base = 2) + (n < 2 ? 0.0 : log2(min(n, S))))
        end
        total = max(ndigits(D; base = 2) + bits, Jbits + Cbits + ceil(Int, wbits + log2(d)) + 1) + 2
        ndigits(M; base = 2) > total && return w, D, i
        need = total
    end
end

"""
Do the terms of a rule cancel exactly? One nonzero residue settles a nonzero sum, usually at the first place of
the first prime, and the bound settles a vanishing one.
"""
function _level_sum_zero(s::FactorialSum, zs::Vector{Int}, es, cs, h::Int)
    F = _level_field(h)
    front = zeros(Int, length(first(es)))
    L = _level_sum(s, zs, front, h)
    iszero(_place_value(_places(F, 1, L.top), L, 1, h)) || return false
    need = first(_level_sum_bound(L, front, es, cs, h))
    bits = 0; i = 0
    while bits <= need
        i += 1
        t = _places(F, i, L.top)
        for a in 1:F.d
            iszero(_place_value(t, L, a, h)) || return false
        end
        bits += ndigits(t.m.p; base = 2) - 1
    end
    return true
end

"""
A value at level k: (Σ a_n [n]) √(r ∏[n]^f), with n ≤ d in the sum, exponents f = 1 under the root, and the
square part of the radicand already taken into the sum.
"""
struct LevelValue
    k::Int
    coeffs::Vector{Rational{BigInt}}
    radcoef::Rational{BigInt}
    rad::Vector{Pair{Int,Int}}
end

Base.iszero(v::LevelValue) = all(iszero, v.coeffs)

"The value of a rule at a level, in q-integers."
function level_value(s::FactorialSum, k::Integer)
    try
        return _level_value(s, Int(k), Rational{Int})
    catch e
        (e isa OverflowError || e isa InexactError) || rethrow()
        return _level_value(s, Int(k), Rational{BigInt})
    end
end

function _level_value(s::FactorialSum, k::Int, ::Type{R}) where {R}
    h = k + 2; d = euler_phi(2h) ÷ 2
    zero_value = LevelValue(k, fill(0 // big(1), d), 1 // big(1), Pair{Int,Int}[])
    data = _level_exponents(s, k, R)
    data === nothing && return zero_value
    zs, pre, prerat, es, cs = data
    front = [fld(a, 2) for a in pre]
    rad = [n => 1 for n in eachindex(pre) if isodd(pre[n])]
    L = _level_sum(s, zs, front, h)
    w, D, _ = _level_coordinates(L, _level_sum_bound(L, front, es, cs, h)..., h)
    all(iszero, w) && return zero_value
    outside, under = _split_square(prerat)
    scale = Int(s.sign0) * Rational{BigInt}(outside) / D
    return LevelValue(k, [scale * a for a in w], Rational{BigInt}(under), rad)
end

function numeric_value(v::LevelValue; bits::Integer = precision(BigFloat))
    w = setprecision(BigFloat, Int(bits) + 32) do
        h = v.k + 2
        qn(n) = sinpi(BigFloat(n) / h) / sinpi(inv(BigFloat(h)))
        total = sum(BigFloat(c) * qn(n) for (n, c) in enumerate(v.coeffs))
        total * sqrt(BigFloat(v.radcoef) * prod(qn(n)^e for (n, e) in v.rad; init = BigFloat(1)))
    end
    return setprecision(() -> BigFloat(w), BigFloat, Int(bits))
end
Base.float(v::LevelValue) = Float64(numeric_value(v; bits = 128))
