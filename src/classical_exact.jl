# Exact q=1 specialization of a factorial rule. No cyclotomic representation or
# floating-point zero screen is involved. Only immutable prime tables are cached.
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

"Reduced factorial product via Legendre valuations; cancellation precedes BigInt products."
function _classical_factorial_ratio(pairs)
    num = big(1); den = big(1)
    isempty(pairs) && return num,den
    for p in _classical_primes(last(pairs).first)
        p > last(pairs).first && break
        e = 0
        for (n,c) in pairs
            e = Base.checked_add(e,Base.checked_mul(c,_factorial_valuation(n,p)))
        end
        iszero(e) && continue
        dest = e > 0 ? num : den
        power = abs(e)
        if power == 1
            MPZ.mul_si!(dest,dest,p)
        elseif power == 2 && p <= isqrt(typemax(Int))
            MPZ.mul_si!(dest,dest,p*p)
        else
            MPZ.mul!(dest,dest,big(p)^power)
        end
    end
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
        g = gcd(P,Q)
        P = div(P,g); Q = div(Q,g)
    end
    return P,Q
end

"""
    classical_exact(s::FactorialSum) -> ClassicalResult

Evaluate a factorial rule exactly at q=1. The finite sum uses integer Horner nesting;
Legendre valuations cancel factorial powers in the squared prefactor and initial
term before materializing integers. The result is a sign and a reduced squared
rational value. No DCR, approximate arithmetic, or modular zero decision is used.
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
    r = P//Q
    sq = (num//den)*r^2
    sg = s.sign0 * (s.alternating && isodd(s.zlo) ? -1 : 1) * sign(P)
    return ClassicalResult(sg,sq)
end
