# Exact q=1 specialization of a factorial rule. Immutable prime tables are cached.
const CLASSICAL_EXACT_PRIMES = LevelCache{Vector{Int}}()

function _classical_primes(N::Int)
    b = _classical_bucket(N)
    return get_level!(CLASSICAL_EXACT_PRIMES,b) do
        limit = 1 << b
        sieve = trues(limit)
        sieve[1] = false
        for p in 2:isqrt(limit)
            sieve[p] || continue
            for m in p*p:p:limit
                sieve[m] = false
            end
        end
        findall(sieve)
    end
end

@inline function _factorial_valuation(n::Int,p::Int)
    e = 0
    while n >= p
        n ÷= p
        e += n
    end
    return e
end

"Combine factorial powers in the squared prefactor times the squared initial term."
function _classical_square_factors(s::FactorialSum)
    pairs = Pair{Int,Int}[]
    sizehint!(pairs,length(s.pre)+length(s.fac))
    scale = s.sqrt_pre ? 1 : 2
    for (n,c) in s.pre
        n > 1 && c != 0 && push!(pairs,Int(n)=>scale*Int(c))
    end
    for f in s.fac
        n = _arg(f,s.zlo)
        n > 1 && f.c != 0 && push!(pairs,n=>2Int(f.c))
    end
    sort!(pairs;by=first)
    w = 0; i = 1
    while i <= length(pairs)
        n,c = pairs[i]; i += 1
        while i <= length(pairs) && pairs[i].first == n
            c = Base.checked_add(c,pairs[i].second)
            i += 1
        end
        if c != 0
            w += 1
            pairs[w] = n=>c
        end
    end
    resize!(pairs,w)
    return pairs
end

"""
Reduced factorial product via Legendre valuations; cancellation precedes BigInt products. Prime powers
are gathered in a machine word and multiplied into the BigInt only when the word is full, instead of one
`big(p)^e` per prime. `num` and `den` have disjoint prime supports, so they are coprime by construction.
"""
function _classical_factorial_ratio(pairs)
    num = big(1); den = big(1)
    isempty(pairs) && return num,den
    N = last(pairs).first
    wn = UInt64(1); wd = UInt64(1)
    for p in _classical_primes(N)
        p > N && break
        e = 0
        for (n,c) in pairs
            e = Base.checked_add(e,Base.checked_mul(c,_factorial_valuation(n,p)))
        end
        iszero(e) && continue
        pp = UInt64(p)
        for _ in 1:abs(e)
            if e > 0
                w, over = Base.mul_with_overflow(wn, pp)
                over ? (MPZ.mul_ui!(num, wn); wn = pp) : (wn = w)
            else
                w, over = Base.mul_with_overflow(wd, pp)
                over ? (MPZ.mul_ui!(den, wd); wd = pp) : (wd = w)
            end
        end
    end
    MPZ.mul_ui!(num, wn); MPZ.mul_ui!(den, wd)
    return num,den
end

"Exact Horner fallback when an adjacent ratio does not fit a machine integer."
function _horner_sum_big(s::FactorialSum)
    P = big(1); Q = big(1)
    for z in (s.zhi-1):-1:s.zlo
        a = big(s.alternating ? -1 : 1); b = big(1)
        for f in s.fac
            lo,hi,c = _factor_step(f,z)
            (c == 0 || lo > hi) && continue
            # Compute each interval product once, then raise it to its integer power.
            block = big(1)
            for n in lo:hi
                MPZ.mul_si!(block,block,n)
            end
            power = block^abs(c)
            c > 0 ? MPZ.mul!(a,a,power) : MPZ.mul!(b,b,power)
        end
        P = b*Q + a*P
        Q *= b
    end
    # One gcd reduction at the end. 
    g = gcd(P,Q)
    return div(P,g), div(Q,g)
end

"""
    classical_exact(s::FactorialSum) -> ClassicalResult

Evaluate a factorial rule exactly at q=1. The finite sum uses integer Horner nesting. 
The result is a sign and a reduced squared rational value.
"""
function classical_exact(s::FactorialSum)
    is_empty_sum(s) && return zero(ClassicalResult)
    P,Q = try
        _horner_sum(s)
    catch e
        e isa OverflowError || rethrow()
        _horner_sum_big(s)
    end
    iszero(P) && return zero(ClassicalResult)
    num,den = _classical_factorial_ratio(_classical_square_factors(s))
    # num/den and P/Q are each coprime once P/Q is reduced
    g = gcd(P,Q)
    isone(g) || (P = div(P,g); Q = div(Q,g))
    P2 = P*P; Q2 = Q*Q
    g1 = gcd(num,Q2); isone(g1) || (num = div(num,g1); Q2 = div(Q2,g1))
    g2 = gcd(den,P2); isone(g2) || (den = div(den,g2); P2 = div(P2,g2))
    sq = Rational{BigInt}(num*P2, den*Q2)
    sg = s.sign0 * (s.alternating && isodd(s.zlo) ? -1 : 1) * sign(P)
    return ClassicalResult(sg,sq)
end
