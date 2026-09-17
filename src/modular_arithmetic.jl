
# ---------------------------------------------------------------------------------
#   Word-size modular arithmetic for exact zero tests
#
#   Montgomery multiplication modulo odd primes p < 2^63, primes p ≡ 1 (mod n) found with a
#   deterministic Miller–Rabin test, and primitive n-th roots of unity modulo p. No dependencies.
# ---------------------------------------------------------------------------------


"Montgomery arithmetic modulo an odd prime p < 2^63 (values kept in Montgomery form a·2^64 mod p)."
struct Montgomery
    p::UInt64
    pinv::UInt64      # -p⁻¹ mod 2^64
    r2::UInt64        # 2^128 mod p
    one::UInt64       # 2^64 mod p
end

function Montgomery(p::Integer)
    p = UInt64(p)
    (isodd(p) && p < (UInt64(1) << 63)) || throw(ArgumentError("Montgomery arithmetic needs an odd modulus below 2^63"))
    inv = p                           # p·p ≡ 1 (mod 8); each Newton step doubles the correct bits
    for _ in 1:6
        inv *= 2 - p * inv
    end
    r1 = UInt64((UInt128(1) << 64) % p)
    return Montgomery(p, -inv, UInt64((UInt128(r1) * r1) % p), r1)
end

@inline function mont_mul(m::Montgomery, a::UInt64, b::UInt64)
    t = UInt128(a) * b
    q = (t % UInt64) * m.pinv
    u = UInt64((t + UInt128(q) * m.p) >> 64)
    return u >= m.p ? u - m.p : u
end

@inline mont_add(m::Montgomery, a::UInt64, b::UInt64) = (s = a + b; s >= m.p ? s - m.p : s)
@inline mont_sub(m::Montgomery, a::UInt64, b::UInt64) = a >= b ? a - b : a + m.p - b

to_mont(m::Montgomery, a::Integer) = mont_mul(m, UInt64(mod(Int128(a), Int128(m.p))), m.r2)
from_mont(m::Montgomery, a::UInt64) = mont_mul(m, a, UInt64(1))

function mont_pow(m::Montgomery, a::UInt64, e::Integer)
    r = m.one
    while e > 0
        isodd(e) && (r = mont_mul(m, r, a))
        a = mont_mul(m, a, a)
        e >>= 1
    end
    return r
end

mont_inv(m::Montgomery, a::UInt64) = mont_pow(m, a, m.p - 2)


# ---- primes ----

@inline _mulmod(a::UInt64, b::UInt64, n::UInt64) = UInt64((UInt128(a) * b) % n)

function _powmod(a::UInt64, e::UInt64, n::UInt64)
    r = UInt64(1)
    a %= n
    while e > 0
        isodd(e) && (r = _mulmod(r, a, n))
        a = _mulmod(a, a, n)
        e >>= 1
    end
    return r
end

const _MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

"Deterministic primality test for 64-bit integers (Miller–Rabin with the first twelve prime bases)."
function is_prime_u64(n::UInt64)
    n < 2 && return false
    for b in _MR_BASES
        n == b && return true
        n % b == 0 && return false
    end
    d = n - 1
    s = 0
    while iseven(d)
        d >>= 1
        s += 1
    end
    for b in _MR_BASES
        x = _powmod(UInt64(b), d, n)
        (x == 1 || x == n - 1) && continue
        witness = true
        for _ in 1:s-1
            x = _mulmod(x, x, n)
            if x == n - 1
                witness = false
                break
            end
        end
        witness && return false
    end
    return true
end

"The `count` largest primes p ≡ 1 (mod n) below `below`, largest first."
function primes_1mod(n::Int, count::Int; below::UInt64 = UInt64(1) << 62)
    ps = UInt64[]
    t = (below - 1) ÷ UInt64(n)
    while length(ps) < count
        p = 1 + UInt64(n) * t
        is_prime_u64(p) && push!(ps, p)
        t -= 1
    end
    return ps
end

"Distinct prime factors of a small positive integer."
function _prime_factors(n::Int)
    fs = Int[]
    m = n
    d = 2
    while d * d <= m
        if m % d == 0
            push!(fs, d)
            while m % d == 0
                m ÷= d
            end
        end
        d += 1
    end
    m > 1 && push!(fs, m)
    return fs
end

"A primitive n-th root of unity modulo a prime p ≡ 1 (mod n)."
function root_of_unity(p::UInt64, n::Int)
    (p - 1) % UInt64(n) == 0 || throw(ArgumentError("$p is not ≡ 1 (mod $n)"))
    fs = _prime_factors(n)
    e = (p - 1) ÷ UInt64(n)
    for x in UInt64(2):UInt64(1_000_000)
        g = _powmod(x, e, p)
        all(q -> _powmod(g, UInt64(n ÷ q), p) != 1, fs) && return g
    end
    error("no primitive $n-th root of unity found modulo $p")
end
