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
# |σE| < M for every σ forces E = 0, and σE is bounded term by term from the sizes of the conjugates. Those sizes
# are enclosed once per level and rounded outward to fixed-point integers, so the bound itself is integer addition.
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

"Invert a vector of nonzero residues with one modular inverse and linear work."
function _invert_places!(v::Vector{UInt64}, m::Montgomery)
    prefix = similar(v); run = m.one
    for i in eachindex(v)
        prefix[i] = run
        run = mont_mul(m, run, v[i])
    end
    iszero(run) && throw(DivideError())
    run = mont_inv(m, run)
    for i in reverse(eachindex(v))
        old = v[i]
        v[i] = mont_mul(m, run, prefix[i])
        run = mont_mul(m, run, old)
    end
    return v
end

function LevelPlaces(p::UInt64, h::Int, N::Int)
    m = Montgomery(p); r = to_mont(m, root_of_unity(p, 2h))
    pow = Vector{UInt64}(undef, 2h); pow[1] = m.one
    for j in 2:2h
        pow[j] = mont_mul(m, pow[j-1], r)
    end
    as = [a for a in 1:h-1 if gcd(a, 2h) == 1]; d = length(as)
    X = [mont_add(m, pow[a+1], pow[2h-a+1]) for a in as]
    den = _invert_places!([mont_sub(m, pow[a+1], pow[2h-a+1]) for a in as], m)
    nd = [_invert_places!([mont_sub(m, X[i], X[i-j+1]) for i in j:d], m) for j in 2:d]
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

"""
A relation ∏[n]^e · G² = c used to rewrite square classes, with G a field element given by its coordinates, or left
empty for 1. `pivot > 0` is a label whose exponent is odd in `e`; `pivot < 0` is minus a prime whose exponent is odd
in `c`, with every exponent in `e` even.
"""
struct ClassRow
    pivot::Int
    e::Vector{Int}
    c::Rational{BigInt}
    root::Vector{BigInt}
end

"The primes in use at one level and their tables."
mutable struct LevelField
    h::Int
    d::Int
    places::Vector{LevelPlaces}
    cursor::UInt64
    valuation::Vector{Tuple{Int,Rational{Int}}}   # for each [n], n ≤ h/2: the prime ℓ dividing it and its valuation
    basis::Vector{Vector{Pair{Int,BigInt}}} # coordinates on [1], …, [d] of [d+1], …, [h÷2]; empty until asked for
    basisbits::Int              # bits of the largest of those coordinates
    classes::Union{Nothing,Vector{ClassRow}}
    inverses::Dict{Int,Tuple{Vector{BigInt},BigInt}}    # coordinates of 1/[n], as they are asked for
    logsines::Union{Nothing,Tuple{Vector{Int},Vector{Int}}} # `_log_sine_table(h)`, once it is asked for
end
LevelField(h::Int) = LevelField(h, _totient(2h) ÷ 2, LevelPlaces[], ((UInt64(1) << 62) - 1) ÷ UInt64(2h),
                                _label_valuations(h),
                                Vector{Pair{Int,BigInt}}[], 0, nothing, Dict{Int,Tuple{Vector{BigInt},BigInt}}(), nothing)

"For each [n] with n ≤ h/2: `(ℓ, v)` when it is not a unit, with ℓ normalized to valuation 1, and `(0, 0//1)` when it is."
function _label_valuations(h::Int)
    out = fill((0, 0//1), h ÷ 2)
    power_of(o) = (fs = _prime_factors(o); length(fs) == 1 ? fs[1] : 0)
    top = h > 1 ? power_of(h) : 0
    for n in 2:h÷2
        o = h ÷ gcd(n, h)
        l = o > 1 ? power_of(o) : 0
        l == 0 && continue
        v = 1 // _totient(o) - (top == 0 ? 0//1 : 1 // _totient(h))
        v > 0 && (out[n] = (l, v))
    end
    return out
end

const LEVEL_FIELDS = LRU{Int,LevelField}(maxsize = 64)
const LEVEL_FIELD_LOCK = ReentrantLock()

_level_field(h::Int) = lock(LEVEL_FIELD_LOCK) do
    get!(() -> LevelField(h), LEVEL_FIELDS, h)
end

# Mutable tables grow after insertion, so an entry-count limit alone cannot bound retained memory.
const LEVEL_FIELD_CACHE_LIMIT = Ref(64 * 1024^2)
function _trim_level_fields!(current::LevelField)
    sum(Base.summarysize(F) for F in values(LEVEL_FIELDS); init = 0) <= LEVEL_FIELD_CACHE_LIMIT[] && return
    empty!(LEVEL_FIELDS)
    Base.summarysize(current) <= LEVEL_FIELD_CACHE_LIMIT[] && (LEVEL_FIELDS[current.h] = current)
    return nothing
end

"Tables of the i-th prime at this level, with factorials up to at least N."
function _places(F::LevelField, i::Int, N::Int)
    lock(LEVEL_FIELD_LOCK) do
        grew = false
        while length(F.places) < i
            p = 1 + UInt64(2F.h) * F.cursor
            F.cursor -= 1
            if is_prime_u64(p)
                push!(F.places, LevelPlaces(p, F.h, N))
                grew = true
            end
        end
        t = F.places[i]
        if size(t.fact, 1) <= N
            t = _grow(t, F.h, max(N, 2 * (size(t.fact, 1) - 1)))
            F.places[i] = t
            grew = true
        end
        grew && _trim_level_fields!(F)
        return t
    end
end

"The fixed-point logarithms of the sines of this level, kept with its other tables."
function _log_sines(F::LevelField)
    F.logsines === nothing || return F.logsines
    lock(LEVEL_FIELD_LOCK) do
        F.logsines === nothing || return F.logsines
        F.logsines = _log_sine_table(F.h)
        _trim_level_fields!(F)
        return F.logsines
    end
end

"What is needed to evaluate a rule at the places of a level: its factorial arguments, and the labels in front."
struct LevelRule
    args::Matrix{Int}               # one column per term
    powers::Vector{Int}
    negative::Vector{Bool}
    front::Vector{Pair{Int,Int}}    # ∏[n]^e multiplying the sum
    top::Int
end

"""
The sizes of Σ_t c_t^scale ∏[n]^(base[n] + scale·e_t[n]) at its conjugates, term by term: for each conjugate the largest
term in fixed point, then the bits of J and of C. C clears the q-integers below the line: for each prime ℓ, the smallest
power of ℓ whose valuation covers theirs. J clears the rationals, once for each power of `scale`.
"""
function _level_sizes(base::Vector{Int}, scale::Int, es::Vector{Vector{Int}}, cs, h::Int)
    F = _level_field(h)
    m = zeros(Int, length(base))
    for e in es, n in eachindex(m)
        m[n] = max(m[n], -(scale * e[n] + base[n]))
    end
    J = big(1)
    for c in cs
        J = lcm(J, denominator(c))
    end
    need = Dict{Int,Rational{BigInt}}()
    for n in 2:length(m)
        l, v = F.valuation[n]
        l != 0 && m[n] > 0 && (need[l] = get(need, l, big(0)//big(1)) + big(m[n]) * v)
    end
    # C clears every negative prime-ideal valuation. Integer ceilings, no floating threshold.
    Cbits = sum(Int(cld(numerator(v), denominator(v))) * ndigits(l; base = 2) for (l, v) in need; init = 0)
    Jbits = ndigits(J; base = 2)
    # Each log size in fixed point: a positive power takes its upper bound, a negative one its lower bound.
    up, lo = _log_sines(F)
    N = length(base)
    cbits = [scale * (ndigits(abs(numerator(c)); base = 2) - ndigits(denominator(c); base = 2) + 1) for c in cs]
    # Nothing below can wrap around if the largest possible total fits, which keeps the integers a proof.
    weight = maximum(e -> sum(abs(scale * e[n] + base[n]) for n in 2:N; init = 0), es)
    Base.Checked.checked_add(Base.Checked.checked_mul(2weight, max(maximum(abs, up), maximum(abs, lo))),
                             Base.Checked.checked_mul(maximum(abs, cbits), 1 << LOG_FRACTION_BITS))
    above = Vector{Int}(undef, N); below = Vector{Int}(undef, N)
    tops = Int[]
    for a in 1:h-1
        gcd(a, 2h) == 1 || continue
        b = min(a, h - a)
        @inbounds for n in 2:N
            r = rem(a * n, h); r = min(r, h - r)
            above[n] = up[r] - lo[b]; below[n] = lo[r] - up[b]
        end
        top = typemin(Int)
        for t in eachindex(es)
            iszero(cs[t]) && continue
            e = es[t]; total = cbits[t] << LOG_FRACTION_BITS
            @inbounds for n in 2:N
                power = scale * e[n] + base[n]
                total += power * (power > 0 ? above[n] : below[n])
            end
            top = max(top, total)
        end
        push!(tops, top)
    end
    return tops, Jbits, Cbits
end

"""
Bits of J·C·front·S at its largest conjugate, of J and of C, each rounded up. The factorial table already holds the
rational in each term, so J appears once.
"""
function _level_sum_bound(L::LevelRule, front::Vector{Int}, es::Vector{Vector{Int}}, cs, h::Int)
    tops, Jbits, Cbits = _level_sizes(front, 1, es, cs, h)
    bits = iszero(count(!iszero, cs)) ? 0 : cld(maximum(tops), 1 << LOG_FRACTION_BITS) + ndigits(length(es); base = 2)
    return max(bits, 0) + Jbits + Cbits, Jbits, Cbits
end

function _level_rule(s::FactorialSum, zs::Vector{Int}, front::Vector{Int}, h::Int)
    fs = collect(s.fac)
    args = [_arg(f, z) for f in fs, z in zs]
    return LevelRule(args, [Int(f.c) for f in fs], [s.alternating && isodd(z) for z in zs],
                    _nonzero_pairs(front), maximum(args; init = 1))
end

"The value at place `i`, in Montgomery form."
function _place_value(t::LevelPlaces, L::LevelRule, i::Int, h::Int)
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
function _level_coordinates(L::LevelRule, bits::Int, Jbits::Int, Cbits::Int, h::Int)
    F = _level_field(h); d = F.d
    acc = [big(0) for _ in 1:d]; M = big(1)
    c = Vector{UInt64}(undef, d); p = Vector{UInt64}(undef, d)
    need = bits
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
        # Every conjugate satisfies |[n]| <= n. This integer bound also covers the candidate.
        wbound = sum(abs(w[n]) * n for n in 1:d)
        total = max(ndigits(D; base = 2) + bits, Jbits + Cbits + ndigits(wbound; base = 2)) + 2
        ndigits(M; base = 2) > total && return w, D, i
        need = total
    end
end

# ---- deciding without expanding ----
#
# A nonzero algebraic integer E has a nonzero norm, and each place at which E vanishes is a prime ideal of norm p dividing
# E, so such places have a product of at most ∏_σ |σE|. Once the primes of the places read pass that, E = 0. The sizes of
# the conjugates add up here, where a coordinate's bound needs the largest of them: a unit is large at one conjugate and
# small at another, and in the norm it counts for nothing.

"""
Do the terms of a rule cancel exactly? One nonzero residue settles a nonzero sum, usually at the first place of
the first prime, and the bound settles a vanishing one.
"""
function _level_sum_zero(s::FactorialSum, zs::Vector{Int}, es, cs, h::Int)
    F = _level_field(h)
    front = zeros(Int, length(first(es)))
    L = _level_rule(s, zs, front, h)
    iszero(_place_value(_places(F, 1, L.top), L, 1, h)) || return false
    tops, Jbits, Cbits = _level_sizes(front, 1, es, cs, h)
    need = _norm_bits(tops, ndigits(length(es); base = 2) + Jbits + Cbits)
    bits = 0; i = 0
    while bits <= need
        i += 1
        t = _places(F, i, L.top)
        for a in 1:F.d
            iszero(_place_value(t, L, a, h)) || return false
            bits += ndigits(t.m.p; base = 2) - 1
        end
    end
    return true
end

"Bits of the norm: the fixed-point sizes `tops` at the conjugates, each with `extra` bits more."
_norm_bits(tops::Vector{Int}, extra::Int) =
    Int(cld(sum(Int128, tops; init = Int128(0)), Int128(1) << LOG_FRACTION_BITS)) + length(tops) * extra

"""
A value read through its rule, c·s with the weights (w, e2), for what is decided at a level by squares: its contributing
terms and its rational factor. Its square lies in the field of q, unlike the value.
"""
struct LevelSquare
    s::FactorialSum
    zs::Vector{Int}
    w::Int
    e2::Int
    c::Rational{BigInt}
end
LevelSquare(s::FactorialSum, k::Int, weight::Tuple{Int,Int} = (0, 0), c = 1) = LevelSquare(s,
    Int[z for seg in last(classify_at_level(s, k)) for z in seg], weight..., Rational{BigInt}(c))

"The square at the place q ↦ q^a, from factorials alone: nothing here is folded, merged or reduced."
_place_square(t::LevelPlaces, V::LevelSquare, i::Int, a::Int, h::Int) =
    _rule_square(t.m, V, V.s, t.fact, t.invf, i, x -> t.pow[mod(a * x, 2h)+1])

"The square of `V` from tables of factorials in column `i`, with `power(x)` the residue of q^x."
function _rule_square(m::Montgomery, V::LevelSquare, s::FactorialSum, fact, invf, i::Int, power)
    acc = zero(UInt64)
    @inbounds for z in V.zs
        u = m.one
        for f in s.fac
            n = _arg(f, z); c = Int(f.c)
            b = c > 0 ? fact[n+1, i] : invf[n+1, i]
            for _ in 1:abs(c)
                u = mont_mul(m, u, b)
            end
        end
        iszero(V.w) || (u = mont_mul(m, u, power(V.w * z)))
        acc = s.alternating && isodd(z) ? mont_sub(m, acc, u) : mont_add(m, acc, u)
    end
    iszero(acc) && return acc
    acc = mont_mul(m, acc, acc)
    @inbounds for (n, c) in s.pre
        e = s.sqrt_pre ? Int(c) : 2Int(c)
        b = e > 0 ? fact[Int(n)+1, i] : invf[Int(n)+1, i]
        for _ in 1:abs(e)
            acc = mont_mul(m, acc, b)
        end
    end
    iszero(V.e2) || (acc = mont_mul(m, acc, power(V.e2)))
    if !isone(V.c)
        p = big(m.p)
        c = mont_mul(m, to_mont(m, mod(numerator(V.c), p)), mont_inv(m, to_mont(m, mod(denominator(V.c), p))))
        acc = mont_mul(m, acc, mont_mul(m, c, c))
    end
    return acc
end

"""
Do two generic values have different squares? Their squares are rational functions of q, read here at one point modulo
one prime, so a different residue proves them different. Equal residues say nothing, and neither does a point where some
[n] vanishes.
"""
function _generic_squares_differ(f, g)
    V = LevelSquare(f.rule, collect(f.rule.zlo:f.rule.zhi), f.weight..., big(1)//big(1))
    W = LevelSquare(g.rule, collect(g.rule.zlo:g.rule.zhi), g.weight..., big(1)//big(1))
    N = max(_rule_top(f.rule, [f.rule.zlo, f.rule.zhi]), _rule_top(g.rule, [g.rule.zlo, g.rule.zhi]))
    m = Montgomery(UInt64(0x7fffffffffffffe7))
    Q = to_mont(m, 0x1234567); Qinv = mont_inv(m, Q); x = mont_add(m, Q, Qinv)
    labels = Vector{UInt64}(undef, N)
    last = zero(UInt64); now = m.one
    for n in 1:N
        iszero(now) && return false
        labels[n] = now
        last, now = now, mont_sub(m, mont_mul(m, x, now), last)
    end
    fact = Matrix{UInt64}(undef, N + 1, 1); invf = similar(fact)
    fact[1, 1] = m.one
    for n in 1:N
        fact[n+1, 1] = mont_mul(m, fact[n, 1], labels[n])
    end
    invf[N+1, 1] = mont_inv(m, fact[N+1, 1])
    for n in N:-1:1
        invf[n, 1] = mont_mul(m, invf[n+1, 1], labels[n])
    end
    power(a) = a >= 0 ? mont_pow(m, Q, a) : mont_pow(m, Qinv, -a)
    return _rule_square(m, V, V.s, fact, invf, 1, power) != _rule_square(m, W, W.s, fact, invf, 1, power)
end

"""
Bits an integer K needs for K·v² to be an algebraic integer, and the fixed-point sizes of K·v² at the conjugates. Only
called once two squares agree at a first place, since this is the one step here that reads the folded rule.
"""
function _square_sizes(V::LevelSquare, k::Int)
    h = k + 2
    data = try
        _level_exponents(V.s, k, Rational{Int})
    catch e
        (e isa OverflowError || e isa InexactError) || rethrow()
        _level_exponents(V.s, k, Rational{BigInt})
    end
    zs, pre, prerat, es, cs = data
    tops, Jbits, Cbits = _level_sizes(pre, 2, es, cs, h)
    r = Rational{BigInt}(prerat) * V.c^2
    Kbits = 2Jbits + Cbits + ndigits(denominator(r); base = 2)
    extra = 2ndigits(length(es); base = 2) + 2Jbits + Cbits + ndigits(numerator(r); base = 2)
    return Kbits, tops .+ (extra << LOG_FRACTION_BITS)
end

"Do two squares differ at the first place? That is a proof that they differ, and usually the whole of the work."
function _first_place_differs(V::LevelSquare, W::LevelSquare, k::Int)
    h = k + 2
    t = _places(_level_field(h), 1, _squares_top(V, W))
    return _place_square(t, V, 1, t.as[1], h) != _place_square(t, W, 1, t.as[1], h)
end
_squares_top(V::LevelSquare, W::LevelSquare) =
    max(_rule_top(V.s, isempty(V.zs) ? [V.s.zlo] : V.zs), _rule_top(W.s, isempty(W.zs) ? [W.s.zlo] : W.zs))

"""
    _level_squares_equal(V, W, k) -> Bool

Whether two values have equal squares at a level. With K, K′ clearing them, E = K′K(v² − w²) is an algebraic integer, of
degree φ(2h) when powers of q occur and φ(2h)/2 when not, so the places needed are those of the norm bound above.
"""
function _level_squares_equal(V::LevelSquare, W::LevelSquare, k::Int)
    h = k + 2; F = _level_field(h)
    N = _squares_top(V, W)
    KV, a = _square_sizes(V, k); KW, b = _square_sizes(W, k)
    tops = [max(a[i] + (KW << LOG_FRACTION_BITS), b[i] + (KV << LOG_FRACTION_BITS)) for i in eachindex(a)]
    weighted = !iszero(V.w) || !iszero(W.w) || !iszero(V.e2) || !iszero(W.e2)
    need = _norm_bits(tops, 1) * (weighted ? 2 : 1)
    bits = 0; i = 0
    while bits <= need
        i += 1
        t = _places(F, i, N)
        for j in eachindex(t.as)
            for a in (weighted ? (t.as[j], 2h - t.as[j]) : (t.as[j],))
                _place_square(t, V, j, a, h) == _place_square(t, W, j, a, h) || return false
                bits += ndigits(t.m.p; base = 2) - 1
            end
        end
    end
    return true
end

"""
    LevelValue

A value at level k: (Σ_n num[n] [n]) / den · √(r ∏_{n ∈ rad} [n]), with n ≤ d = φ(2h)/2 in the sum. The root is
positive, because r > 0 and [n] > 0 for n ≤ h/2, and (r, rad) is reduced by the relations of the level, so that
equal roots carry equal keys as far as those relations see.
"""
struct LevelValue
    k::Int
    num::Vector{BigInt}
    den::BigInt
    r::BigInt
    rad::Vector{Int}
end

Base.iszero(v::LevelValue) = all(iszero, v.num)
coordinates(v::LevelValue) = [a // v.den for a in v.num]

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
    h = k + 2
    data = _level_exponents(s, k, R)
    data === nothing && return zero(LevelValue, k)
    zs, pre, prerat, es, cs = data
    front, outside, r, rad, roots = _canonical_class(copy(pre), Rational{BigInt}(prerat), h)
    L = _level_rule(s, zs, front, h)
    w, D, _ = _level_coordinates(L, _level_sum_bound(L, front, es, cs, h)..., h)
    all(iszero, w) && return zero(LevelValue, k)
    F = _level_field(h)
    for G in roots
        w = _label_mul(w, G, F)
    end
    c = Int(s.sign0) * outside
    return _level_make(k, w .* numerator(c), D * denominator(c), r, rad)
end
