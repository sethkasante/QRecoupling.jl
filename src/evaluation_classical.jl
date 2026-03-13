
#evaluation_classical.jl


# --- Classical Limit computations ---


# --- Precomputations (Sieve) ---
const CLASSICAL_SIEVE = Int[]

struct ClassicalResult
    sign::Int
    sq_val::Rational{BigInt} # The value under the square root
end

function Base.show(io::IO, res::ClassicalResult)
    res.sign == 0 && return print(io, "0.0")
    
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    # Check if it's a perfect square (A purely rational 6j-symbol)
    if s_num^2 == num && s_den^2 == den
        print(io, res.sign < 0 ? "-" : "", Rational(s_num, s_den))
    else
        # It is NOT a perfect square, show as exact square root
        print(io, res.sign < 0 ? "-" : "", "√(", res.sq_val, ")")
    end
end

# 1. Conversion to Floats 
Base.Float64(res::ClassicalResult) = res.sign * sqrt(Float64(res.sq_val))
Base.BigFloat(res::ClassicalResult) = res.sign * sqrt(BigFloat(res.sq_val))

# 2. Conversion to Exact Rational (Only works if it's a perfect square)
function Base.Rational(res::ClassicalResult)
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    if s_num^2 == num && s_den^2 == den
        return res.sign * (s_num // s_den)
    else
        throw(DomainError(res.sq_val, "This symbol is an irrational square root and cannot be represented as a pure Rational."))
    end
end

# 3. Arithmetic Overloads for seamless integration with numerical scripts
Base.:+(res::ClassicalResult, x::Number) = Float64(res) + x
Base.:+(x::Number, res::ClassicalResult) = x + Float64(res)
Base.:*(res::ClassicalResult, x::Number) = Float64(res) * x
Base.:*(x::Number, res::ClassicalResult) = x * Float64(res)


# ============================================================
# Classical Evaluation (q = 1) Limit
# ============================================================

"""
    classical_sieve!(max_d::Int)

Maintains a Sieve of Eratosthenes-style cache where:
- sieve[d] = p  if d = p^k (prime power)
- sieve[d] = 1  otherwise
Used for O(1) evaluation of Φ_d(1).
"""
function classical_sieve!(max_d::Int)
    curr_len = length(CLASSICAL_SIEVE)
    if max_d > curr_len
        new_size = max(max_d, 2 * curr_len, 5000) # Pre-allocate aggressively
        resize!(CLASSICAL_SIEVE, new_size)
        for i in (curr_len+1):new_size
            if i < 2
                CLASSICAL_SIEVE[i] = 0
            else
                # Factorization determines the value of the cyclotomic polynomial at q=1
                facs = collect(Nemo.factor(i))
                # Φ_d(1) = p if d is a power of prime p, else 1
                CLASSICAL_SIEVE[i] = (length(facs) == 1) ? Int(facs[1][1]) : 1
            end
        end
    end
    return nothing
end

"""
    evaluate_classical(m::CycloMonomial)

Evaluates the monomial at q=1 to a floating point number using the prime power sieve.
"""
function evaluate_classical(m::CycloMonomial)
    m.sign == 0 && return 0.0
    classical_sieve!(length(m.exps))
    
    log_val = 0.0
    @inbounds for d in 2:length(m.exps)
        e = m.exps[d]
        e == 0 && continue
        
        p = CLASSICAL_SIEVE[d]
        if p > 1
            log_val += e * log(Float64(p))
        end
    end
    return m.sign * exp(log_val)
end

"""
    evaluate_to_rational(m::CycloMonomial)

Converts a CycloMonomial to an exact Rational{BigInt} at q=1.
Highly optimized to bypass Dict allocations and use fast GMP paths.
"""
function evaluate_to_rational(m::CycloMonomial)
    m.sign == 0 && return 0//1
    classical_sieve!(length(m.exps))
    
    num = BigInt(1)
    den = BigInt(1)
    
    @inbounds for d in 2:length(m.exps)
        e = m.exps[d]
        e == 0 && continue
        
        p = CLASSICAL_SIEVE[d]
        p == 1 && continue
        
        # Fast-paths for e = ±1 (Extremely common in hypergeometric ratios)
        # Avoids allocating a `BigInt` just to hold the prime `p`
        if e == 1
            num *= p
        elseif e == -1
            den *= p
        elseif e > 1
            num *= BigInt(p)^e
        else
            den *= BigInt(p)^abs(e)
        end
    end
    
    return m.sign > 0 ? (num // den) : -(num // den)
end

# ============================================================
# Main Classical Evaluators
# ============================================================

export qracah6j_classical_exact, qracah3j_classical_exact

function qracah6j_classical_exact(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res = q6j_generic(j1, j2, j3, j4, j5, j6)
    
    # Structural zero
    (res.pref_sq.sign == 0 || res.m_min.sign == 0) && return ClassicalResult(0, 0//1)

    # 1. Evaluate Δ² to Rational
    pref_sq_rat = evaluate_to_rational(res.pref_sq)
    
    # 2. Evaluate Sum to Rational (Hypergeometric Pattern)
    m_min_rat = evaluate_to_rational(res.m_min)
    
    sum_val = m_min_rat
    curr_term = m_min_rat
    
    for r in res.ratios
        curr_term *= evaluate_to_rational(r)
        sum_val += curr_term
    end
    
    # Final = sign(sum) * sqrt(pref_sq * sum^2)
    final_sign = sign(sum_val)
    final_sq_val = pref_sq_rat * (sum_val^2)
    
    return ClassicalResult(final_sign, final_sq_val)
end

function qracah3j_classical_exact(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    res = q3j_generic(j1, j2, j3, m1, m2, m3)
    
    # Structural zero
    (res.pref_sq.sign == 0 || res.m_min.sign == 0) && return ClassicalResult(0, 0//1)

    pref_sq_rat = evaluate_to_rational(res.pref_sq)
    m_min_rat = evaluate_to_rational(res.m_min)
    
    sum_val = m_min_rat
    curr_term = m_min_rat
    
    for r in res.ratios
        curr_term *= evaluate_to_rational(r)
        sum_val += curr_term
    end
    
    final_sign = sign(sum_val)
    final_sq_val = pref_sq_rat * (sum_val^2)
    
    return ClassicalResult(final_sign, final_sq_val)
end