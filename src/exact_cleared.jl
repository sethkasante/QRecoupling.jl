# ---------------------------------------------------------------------------------
#  Exact level values without field arithmetic
#
#  A ratio of consecutive Racah terms is a ratio of products of q-integers, and for the standard symbols
#  it has as many factors up as down, so the (q − q⁻¹) powers cancel identically:
#
#      t_{z+1}/t_z  =  ∏ (q^{a_i} − q^{-a_i}) / ∏ (q^{b_i} − q^{-b_i}).
#
#  Work in ℤ[q]/(q^{2h} − 1) rather than in ℚ(ζ_{2h}). There, multiplying by (qⁿ − q⁻ⁿ) is
#
#      out[i] = c[i−n] − c[i+n],
#
#  two shifted subtractions: no multiplication, no rational coefficients, and no reduction modulo the
#  cyclotomic polynomial. The cleared Horner P ← B·Q + A·P, Q ← B·Q therefore runs in machine words modulo
#  a few primes, and the CRT is *exact* — P and Q are integer polynomials, so there is no denominator to
#  reconstruct and the number of primes follows from a closed-form coefficient bound.
#
#  Both are palindromic: each ratio multiplies by an even total number of antisymmetric binomials when the
#  factor counts match, so P(1/q) = ±P(q) with a parity that the recursion tracks. Only the coefficients of
#  q⁰ … q^h are stored.
#
#  In a benchmark of {30 30 30 30 30 30} at k = 200, the nesting cost
#  1.4 ms here against 14.3 ms for the same Horner over field elements and 286 ms for the term-by-term
#  projector. What remains expensive is the *final division* of two field elements (19 ms), which only
#  canonicalises; a cleared result that is divided lazily keeps the whole factor.
# ---------------------------------------------------------------------------------

"""
Distinct primes just below 2^62 for the multimodular Horner. Unlike the zero test these need no root of
unity — only primality, because the reconstruction is of *integer* polynomials — but they must really be
prime: Garner inverts each modulo the others. Generated with the package's own Miller–Rabin rather than
written down.
"""
const CLEARED_PRIMES = let ps = UInt64[], n = (UInt64(1) << 62) - UInt64(1)
    while length(ps) < 16
        is_prime_u64(n) && push!(ps, n)
        n -= UInt64(2)
    end
    ps
end

"Largest number of q-integer factors on one side of a ratio that the fixed scratch holds."
const CLEARED_MAX_FACTORS = 32

@inline _cl_sub(a::UInt64, b::UInt64, p::UInt64) = a >= b ? a - b : a + p - b
@inline _cl_add(a::UInt64, b::UInt64, p::UInt64) = (c = a + b; c >= p ? c - p : c)

"""
    _ratio_qints!(num, den, s, z) -> (nnum, nden, sign) or nothing

The q-integer indices of `t_{z+1}/t_z`, split by side, written into the caller's buffers. `nothing` when a
side needs more than `CLEARED_MAX_FACTORS` entries.
"""
function _ratio_qints!(num::Vector{Int}, den::Vector{Int}, s::FactorialSum, z::Int)
    nn = 0; nd = 0
    @inbounds for f in s.fac
        lo, hi, c = _factor_step(f, z)
        (c == 0 || lo > hi) && continue
        for n in lo:hi, _ in 1:abs(c)
            if c > 0
                nn += 1; nn > CLEARED_MAX_FACTORS && return nothing
                num[nn] = n
            else
                nd += 1; nd > CLEARED_MAX_FACTORS && return nothing
                den[nd] = n
            end
        end
    end
    return (nn, nd, s.alternating ? -1 : 1)
end

"""
    _cleared_plan(s, h) -> (bits, nsteps) or nothing

Checks the preconditions and returns the coefficient bound. The ratios must be *balanced* (equally many
factors up and down, so the (q - q⁻¹) powers cancel) and no denominator index may be a multiple of `h`,
where the q-integer vanishes: the cleared Horner would then carry a zero denominator. A vanishing
*numerator* is fine — it truncates the sum, exactly as the level does.

`bits` bounds log₂ of any coefficient of P or Q: each step multiplies by at most `max(nnum, nden)`
binomials with unit coefficients and then adds once.
"""
function _cleared_plan(s::FactorialSum, h::Int)
    num = Vector{Int}(undef, CLEARED_MAX_FACTORS)
    den = Vector{Int}(undef, CLEARED_MAX_FACTORS)
    bits = 2
    nsteps = 0
    for z in s.zlo:(s.zhi - 1)
        r = _ratio_qints!(num, den, s, z)
        r === nothing && return nothing
        nn, nd, _ = r
        nn == nd || return nothing                       # unbalanced: (q − 1/q) powers would survive
        @inbounds for i in 1:nd
            iszero(mod(den[i], h)) && return nothing     # a zero denominator cannot be cleared
        end
        bits += max(nn, nd) + 1
        nsteps += 1
    end
    return (bits, nsteps)
end

"""
Coefficient of `q^x` of a folded vector: entries 0…h are stored, and the reflection carries the parity
`sigma` (+1 palindromic, −1 antipalindromic), so index m−i mirrors index i up to sign.
"""
@inline function _cl_at(c::Vector{UInt64}, x::Int, sigma::Int, m::Int, h::Int, p::UInt64)
    r = mod(x, m)
    @inbounds v = r <= h ? c[r+1] : c[m-r+1]
    return (r <= h || sigma > 0) ? v : (iszero(v) ? v : p - v)
end

"out ← (qⁿ − q⁻ⁿ)·c on the folded range; the parity flips."
@inline function _cl_mul_qint!(out::Vector{UInt64}, c::Vector{UInt64}, n::Int,
                               sigma::Int, m::Int, h::Int, p::UInt64)
    @inbounds for i in 0:h
        out[i+1] = _cl_sub(_cl_at(c, i - n, sigma, m, h, p),
                           _cl_at(c, i + n, sigma, m, h, p), p)
    end
    return -sigma
end

"""
    _cleared_horner_mod(s, m, h, p) -> (P, Q, sigma)

Σ/t_first = P/Q in (ℤ/p)[q]/(q^m − 1), folded to indices 0…h, by shifted subtraction only.
"""
function _cleared_horner_mod(s::FactorialSum, m::Int, h::Int, p::UInt64)
    L = h + 1
    P = zeros(UInt64, L); Q = zeros(UInt64, L)
    A = zeros(UInt64, L); B = zeros(UInt64, L); S = zeros(UInt64, L)
    num = Vector{Int}(undef, CLEARED_MAX_FACTORS)
    den = Vector{Int}(undef, CLEARED_MAX_FACTORS)
    @inbounds P[1] = UInt64(1); Q[1] = UInt64(1)
    sigma = 1
    for z in (s.zhi - 1):-1:s.zlo
        nn, nd, sg = _ratio_qints!(num, den, s, z)
        copyto!(A, P); sa = sigma
        @inbounds for i in 1:nn
            sa = _cl_mul_qint!(S, A, num[i], sa, m, h, p); copyto!(A, S)
        end
        copyto!(B, Q); sb = sigma
        @inbounds for i in 1:nd
            sb = _cl_mul_qint!(S, B, den[i], sb, m, h, p); copyto!(B, S)
        end
        # balanced counts make the two parities equal, so the sum has one parity
        @inbounds if sg < 0
            for i in 1:L; P[i] = _cl_sub(B[i], A[i], p); end
        else
            for i in 1:L; P[i] = _cl_add(B[i], A[i], p); end
        end
        copyto!(Q, B)
        sigma = sb
    end
    return P, Q, sigma
end

@inline function _cl_powmod(a::UInt64, e::UInt64, p::UInt64)
    r = UInt64(1); b = a % p; ee = e
    while !iszero(ee)
        isodd(ee) && (r = UInt64(widemul(r, b) % p))
        b = UInt64(widemul(b, b) % p); ee >>= 1
    end
    return r
end
@inline _cl_invmod(a::UInt64, p::UInt64) = _cl_powmod(a, p - UInt64(2), p)

"""
    _cleared_garner!(out, residues, ps, digits, M, half)

Mixed-radix (Garner) CRT of one coefficient position across primes: the digits are machine words and only
the final assembly touches a `BigInt`, in place.
"""
function _cleared_garner!(out::Vector{BigInt}, residues, ps::Vector{UInt64}, inv,
                          digits::Vector{UInt64}, M::BigInt, half::BigInt, which::Int, L::Int)
    np = length(ps)
    acc = BigInt()
    @inbounds for i in 1:L
        for a in 1:np
            r = residues[a][which][i] % ps[a]
            for b in 1:a-1
                r = _cl_sub(r, digits[b] % ps[a], ps[a])
                r = UInt64(widemul(r, inv[a][b]) % ps[a])
            end
            digits[a] = r
        end
        MPZ.set_ui!(acc, digits[np])
        for a in (np-1):-1:1
            MPZ.mul_ui!(acc, ps[a])
            MPZ.add_ui!(acc, digits[a])
        end
        # `MPZ.set` copies: `BigInt(acc)` is the identity, which would alias the accumulator
        out[i] = acc > half ? acc - M : MPZ.set(acc)
    end
    return out
end

"""
    _cleared_nesting(s, k) -> (P, Q, sigma) or nothing

Σ/t_first as a pair of exact integer polynomials in ℤ[q]/(q^{2h} − 1), folded to indices 0…h with parity
`sigma`. The number of primes comes from the coefficient bound, so the reconstruction is exact rather than
probabilistic.
"""
function _cleared_nesting(s::FactorialSum, k::Int)
    h = k + 2
    plan = _cleared_plan(s, h)
    plan === nothing && return nothing
    bits, nsteps = plan
    nsteps == 0 && return nothing                        # a one-term sum has nothing to nest
    # Route by length against field size. The cleared Horner costs O(primes · steps · 2h) word operations
    # whatever the coefficients do, while the term-by-term walk's field elements grow denser with every
    # step; so the cleared form wins once the sum is long relative to the field. Measurements over
    # 20 cases placed the crossover at nsteps ≈ √(2h/6); this heuristic selected the faster route in
    # all of those cases.
    6 * nsteps^2 >= 2h || return nothing
    np = cld(bits + 2, 61)
    np <= length(CLEARED_PRIMES) || return nothing
    m = 2h; L = h + 1
    ps = CLEARED_PRIMES[1:np]
    res = Vector{Tuple{Vector{UInt64},Vector{UInt64},Int}}(undef, np)
    for a in 1:np
        res[a] = _cleared_horner_mod(s, m, h, ps[a])
    end
    sigma = res[1][3]
    if np == 1                                           # one prime: lift to the symmetric range directly
        p = ps[1]; halfp = p >> 1
        lift(v) = BigInt(v > halfp ? Int128(v) - Int128(p) : Int128(v))
        return ([lift(res[1][1][i]) for i in 1:L], [lift(res[1][2][i]) for i in 1:L], sigma)
    end
    inv = [[_cl_invmod(ps[b] % ps[a], ps[a]) for b in 1:a-1] for a in 1:np]
    M = prod(BigInt.(ps)); half = M >> 1
    digits = Vector{UInt64}(undef, np)
    P = Vector{BigInt}(undef, L); Q = Vector{BigInt}(undef, L)
    _cleared_garner!(P, res, ps, inv, digits, M, half, 1, L)
    _cleared_garner!(Q, res, ps, inv, digits, M, half, 2, L)
    return (P, Q, sigma)
end

"Integer content of a coefficient vector."
function _cleared_content(c::Vector{BigInt})
    g = BigInt(0)
    for v in c
        g = gcd(g, v)
        isone(g) && break
    end
    return g
end

"""
Unfold a parity-`sigma` vector on 0…h to the full 2h coefficients, dividing by `g`.

`g` must be the content **common to the numerator and the denominator**: scaling them by different
constants would change the quotient they represent.
"""
function _cleared_unfold(c::Vector{BigInt}, sigma::Int, m::Int, h::Int, g::BigInt)
    out = zeros(BigInt, m)
    @inbounds for i in 0:h
        v = (iszero(g) || isone(g)) ? c[i+1] : div(c[i+1], g)
        out[i+1] = v
        j = mod(-i, m)
        j == i && continue
        out[j+1] = sigma > 0 ? v : -v
    end
    return out
end

"""
    _cleared_exact_sum(s, dcr, k, V_exact, V_inv, ζ, h) -> field element or nothing

The exact Racah sum Σ_z project(root·base·∏ratios) with the ratio walk replaced by the cleared Horner: the
first monomial is projected once, and the nesting comes from two integer polynomials. `nothing` when the
preconditions of `_cleared_nesting` do not hold, so the caller keeps its term-by-term projector.
"""
function _cleared_exact_sum(s::FactorialSum, dcr::DCR, k::Int, V_exact, V_inv, ζ, h::Int)
    nest = _cleared_nesting(s, k)
    nest === nothing && return nothing
    P, Q, sigma = nest
    m = 2h
    K = parent(ζ)
    R, _ = polynomial_ring(ZZ, "x")
    # one content for both sides: the quotient is what represents the sum
    g = gcd(_cleared_content(P), _cleared_content(Q))
    qv = K(R(_cleared_unfold(Q, sigma, m, h, g)))
    iszero(qv) && return nothing
    pv = K(R(_cleared_unfold(P, sigma, m, h, g)))
    # the first term, projected exactly as the term-by-term walk does
    buf = CycloBuffer(dcr.max_d)
    mul!(buf, dcr.root, dcr.base)
    first_mono = snapshot(buf)
    e_rad = _phi_exponent(dcr.radical, h)
    v2 = e_rad + 2 * _phi_exponent(first_mono, h)
    v2 < 0 && throw(DomainError(k, "Topological pole at level k=$k."))
    # With no vanishing denominator the valuations only increase along the sum, so a finite first term
    # guarantees the rest: the nesting cannot introduce a pole.
    v2 > 0 && return zero(ζ)
    first_val = _project_monomial_nemo_internal(first_mono, V_exact, V_inv, ζ, h)
    return first_val * (pv // qv)
end
