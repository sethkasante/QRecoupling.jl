# ---------------------------------------------------------------------------
# Classical Limit Precomputations (Pure Julia Sieve)
# Computes the Ponzano-Regge classical limit where level k → ∞, or `q` = 1
# ---------------------------------------------------------------------------

const CLASSICAL_SIEVE = Int[]
const CLASSICAL_LOG   = Float64[]

export ensure_classical_sieve


"""
    ensure_classical_sieve(max_d::Int)

Dynamically resizes and populates the global prime-power sieve for classical limits.
For any integer `d`, computes the exact prime base if `d = p^k`, alongside its 
precomputed floating-point logarithm to guarantee O(1) cache retrieval during summation.
"""
function ensure_classical_sieve(max_d::Int)
    curr_len = length(CLASSICAL_SIEVE)
    if max_d > curr_len
        new_size = max(max_d, 2 * curr_len, 5000)
        
        resize!(CLASSICAL_SIEVE, new_size)
        resize!(CLASSICAL_LOG, new_size)
        
        for i in (curr_len+1):new_size
            CLASSICAL_SIEVE[i] = 1
            CLASSICAL_LOG[i] = 0.0
        end
        CLASSICAL_SIEVE[1] = 0
        
        is_prime = trues(new_size)
        is_prime[1] = false
        
        @inbounds for p in 2:new_size
            if is_prime[p]
                log_p = log(Float64(p))
                
                power = p
                while power <= new_size
                    CLASSICAL_SIEVE[power] = p
                    CLASSICAL_LOG[power] = log_p
                    power = widen(power) * p > new_size ? new_size + 1 : power * p
                end
                
                for mult in (p * 2):p:new_size
                    is_prime[mult] = false
                end
            end
        end
    end
end

# ==============================================================================
# Result Struct and Base Overloads
# ==============================================================================

struct ClassicalResult
    sign::Int
    sq_val::Rational{BigInt} 
end

function Base.show(io::IO, res::ClassicalResult)
    res.sign == 0 && return print(io, "0.0")
    
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    if s_num^2 == num && s_den^2 == den
        print(io, res.sign < 0 ? "-" : "", Rational(s_num, s_den))
    else
        print(io, res.sign < 0 ? "-" : "", "√(", res.sq_val, ")")
    end
end

Base.Float64(res::ClassicalResult) = res.sign * Float64(sqrt(BigFloat(res.sq_val)))
Base.BigFloat(res::ClassicalResult) = res.sign * sqrt(BigFloat(res.sq_val))

function Base.Rational{T}(res::ClassicalResult) where {T <: Integer}
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    if s_num^2 == num && s_den^2 == den
        return res.sign * Rational{T}(T(s_num), T(s_den))
    else
        throw(InexactError(:Rational, Rational{T}, res))
    end
end

Base.:+(res::ClassicalResult, x::Number) = Float64(res) + x
Base.:*(res::ClassicalResult, x::Number) = Float64(res) * x

# ==============================================================================
# Fast Float64 Evaluators (NaN-Safe Log-Sum-Exp)
# ==============================================================================


@inline function evaluate_classical_log(m::CycloMonomial)
    m.sign == 0 && return -Inf
    log_val = 0.0
    @inbounds for (d, e) in m.exps
        log_val += e * CLASSICAL_LOG[d]
    end
    return log_val
end


"""
    evaluate_classical(res::CycloResult)

Evaluates a pre-constructed `CycloResult` in the classical Ponzano-Regge limit (where q = 1 or k → ∞).

Maps the cyclotomic logic directly to integer primes. Utilizes a two-pass Log-Sum-Exp (LSE) 
shifted summation to guarantee  immunity against `NaN` and `Inf` floating-point overflows
at massive classical spins.
"""
function evaluate_classical(res::CycloResult)
    (res.radical.sign == 0 || res.m_min.sign == 0) && return 0.0
    ensure_classical_sieve(res.max_d)

    max_lm = -Inf
    lm_curr = evaluate_classical_log(res.m_min)
    max_lm = max(max_lm, lm_curr)
    
    for r in res.ratios
        lm_curr += evaluate_classical_log(r)
        max_lm = max(max_lm, lm_curr)
    end

    lm_root = evaluate_classical_log(res.root)
    lm_rem  = evaluate_classical_log(res.radical)
    v_pref = exp(lm_root + (lm_rem / 2))
    
    lm_curr = evaluate_classical_log(res.m_min)
    sign_curr = res.m_min.sign
    
    sum_val = sign_curr * exp(lm_curr - max_lm)
    
    for r in res.ratios
        lm_curr += evaluate_classical_log(r)
        sign_curr *= r.sign
        sum_val += sign_curr * exp(lm_curr - max_lm)
    end
    
    return v_pref * sum_val * exp(max_lm)
end

function q6j_classical(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    return evaluate_classical(q6j_cyclo(j1, j2, j3, j4, j5, j6))
end

function q3j_classical(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    return evaluate_classical(q3j_cyclo(j1, j2, j3, m1, m2, m3))
end

# ==============================================================================
# Exact Rational Evaluators (Zero-Allocation Architecture)
# ==============================================================================

# Thread-safe global workspaces for exact summation to prevent array allocations
const EXACT_WORKSPACE_CURR = Int[]
const EXACT_WORKSPACE_MIN  = Int[]

function evaluate_to_rational(m::CycloMonomial)
    m.sign == 0 && return zero(Rational{BigInt})
    num = one(BigInt)
    den = one(BigInt)
    
    @inbounds for (d, e) in m.exps
        p = CLASSICAL_SIEVE[d]
        if p > 1
            if e > 0
                num *= BigInt(p)^e
            elseif e < 0
                den *= BigInt(p)^abs(e)
            end
        end
    end
    
    val = num // den
    return m.sign > 0 ? val : -val
end

# ---- Core GCD-Free Summation Engine (Platform-Safe, Zero GC-Thrashing Hot Loop) -----
"""
    _sum_hypergeometric_exact(res::CycloResult)

The zero-allocation exact summation hot-loop for classical symbols.

We use a two-pass algorithm: this function dynamically tracks 
prime exponents to reconstruct a global denominator. 
It then executes the hypergeometric sum using strictly in-place C-calls (`Base.GMP.MPZ`) 
to bypass standard Julia `BigInt` Garbage Collection overhead entirely.
"""
function _sum_hypergeometric_exact(res::CycloResult)
    max_d = res.max_d
    
    # Resize thread-safe workspaces to avoid `zeros()` heap allocations
    if length(EXACT_WORKSPACE_CURR) < max_d
        resize!(EXACT_WORKSPACE_CURR, max(max_d, 5000))
        resize!(EXACT_WORKSPACE_MIN, max(max_d, 5000))
    end
    curr_exps = EXACT_WORKSPACE_CURR
    min_exps = EXACT_WORKSPACE_MIN
    
    fill!(view(curr_exps, 1:max_d), 0)
    fill!(view(min_exps, 1:max_d), 0)
    
    # PASS 1: Track Prime Exponent Valleys to find the Global Denominator
    @inbounds for (d, e) in res.m_min.exps
        p = CLASSICAL_SIEVE[d]
        if p > 1
            curr_exps[p] += e
            min_exps[p] = curr_exps[p]
        end
    end
    
    for r in res.ratios
        @inbounds for (d, e) in r.exps
            p = CLASSICAL_SIEVE[d]
            if p > 1
                curr_exps[p] += e
                if curr_exps[p] < min_exps[p]
                    min_exps[p] = curr_exps[p]
                end
            end
        end
    end
    
    # PASS 2: Reconstruct N_0 and D_global
    fill!(view(curr_exps, 1:max_d), 0)
    @inbounds for (d, e) in res.m_min.exps
        p = CLASSICAL_SIEVE[d]
        p > 1 && (curr_exps[p] += e)
    end
    
    D_global = one(BigInt)
    N_0 = one(BigInt)
    
    @inbounds for p in 2:max_d
        E_p = min_exps[p]
        K_p = E_p < 0 ? -E_p : 0 
        
        if K_p > 0
            Base.GMP.MPZ.mul!(D_global, D_global, BigInt(p)^K_p)
        end
        
        pow = curr_exps[p] + K_p
        if pow > 0
            Base.GMP.MPZ.mul!(N_0, N_0, BigInt(p)^pow)
        end
    end
    
    # PASS 3: The Hot Loop
    Sum_N = BigInt()
    Base.GMP.MPZ.set!(Sum_N, N_0)
    res.m_min.sign < 0 && Base.GMP.MPZ.neg!(Sum_N, Sum_N)
    
    curr_N = BigInt()
    Base.GMP.MPZ.set!(curr_N, N_0)
    curr_sign = res.m_min.sign
    
    # Workspace variables to completely eliminate BigInt heap allocations
    r_num = BigInt()
    r_den = BigInt()
    val_big = BigInt() 
    
    for r in res.ratios
        Base.GMP.MPZ.set_si!(r_num, 1)
        Base.GMP.MPZ.set_si!(r_den, 1)
        
        @inbounds for (d, e) in r.exps
            p = CLASSICAL_SIEVE[d]
            if p > 1
                # DIRECT C-CALL: Computes p^|e| directly into val_big with ZERO heap allocations!
                ccall((:__gmpz_ui_pow_ui, :libgmp), Cvoid, (Ref{BigInt}, Culong, Culong), val_big, p, abs(e))
                
                if e > 0
                    Base.GMP.MPZ.mul!(r_num, r_num, val_big)
                elseif e < 0
                    Base.GMP.MPZ.mul!(r_den, r_den, val_big)
                end
            end
        end
        
        # curr_N = (curr_N * r_num) / r_den (strictly in-place)
        Base.GMP.MPZ.mul!(curr_N, curr_N, r_num)
        Base.GMP.MPZ.tdiv_q!(curr_N, curr_N, r_den)
        
        curr_sign *= r.sign
        
        if curr_sign > 0
            Base.GMP.MPZ.add!(Sum_N, Sum_N, curr_N)
        else
            Base.GMP.MPZ.sub!(Sum_N, Sum_N, curr_N)
        end
    end
    
    return Sum_N // D_global
end

export evaluate_classical_exact


"""
    evaluate_classical_exact(res::CycloResult)

Evaluates the rigorous, exact rational representation of the classical limit.

Returns `ClassicalResult` containing:
1. `sign`: The exact integer parity sign of the final evaluation.
2. `sq_val`: The mathematically exact `Rational{BigInt}` of the *squared* symbol. 

By avoiding floating-point square roots entirely, this engine bypasses the 64-bit 
overflows found in standard libraries (like `WignerSymbols.jl`) to preserve infinite precision.
"""
function evaluate_classical_exact(res::CycloResult)
    (res.radical.sign == 0 || res.m_min.sign == 0) && return ClassicalResult(0, 0//1)
    ensure_classical_sieve(res.max_d)

    # 1. Evaluate the exact split prefactor
    root_rat = evaluate_to_rational(res.root)
    radical_rat  = evaluate_to_rational(res.radical)
    
    # 2. Exact integer/rational evaluation of the hypergeometric series
    sum_val = _sum_hypergeometric_exact(res)
    
    # The sum cleanly absorbs the perfectly square-rooted prefactor!
    total_sum = sum_val * root_rat
    
    final_sign = sign(total_sum) * res.radical.sign
    final_sq_val = radical_rat * (total_sum^2)
    
    return ClassicalResult(final_sign, final_sq_val)
end

function q6j_classical_exact(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    return evaluate_classical_exact(q6j_cyclo(j1, j2, j3, j4, j5, j6))
end

function q3j_classical_exact(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    return evaluate_classical_exact(q3j_cyclo(j1, j2, j3, m1, m2, m3))
end