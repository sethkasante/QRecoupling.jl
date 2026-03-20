
# ==============================================================================
# Classical Limit Precomputations (Pure Julia Sieve)
# ==============================================================================
# - CLASSICAL_SIEVE[d] = p  if d = p^k (prime power), else 1
# - CLASSICAL_LOG[d]   = log(p) if d = p^k, else 0.0
const CLASSICAL_SIEVE = Int[]
const CLASSICAL_LOG   = Float64[]

export ensure_classical_sieve

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
                
                # Mark all prime powers p^k
                power = p
                while power <= new_size
                    CLASSICAL_SIEVE[power] = p
                    CLASSICAL_LOG[power] = log_p
                    power = widen(power) * p > new_size ? new_size + 1 : power * p
                end
                
                # Mark multiples as composite
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
# Fast Float64 Evaluators (Log-Space State Tracking)
# ==============================================================================

@inline function evaluate_classical_log(m::CycloMonomial)
    m.sign == 0 && return -Inf
    log_val = 0.0
    @inbounds for (d, e) in m.exps
        log_val += e * CLASSICAL_LOG[d] # O(1) Cache retrieval, zero math!
    end
    return log_val
end


"""
    evaluate_classical(res::CycloResult)

Evaluates a pre-constructed CycloResult in the classical Ponzano-Regge limit (q -> 1).
Uses the O(1) CLASSICAL_LOG Sieve for blisteringly fast floating-point summation.
"""
function evaluate_classical(res::CycloResult)
    res.pref_sq.sign == 0 && return 0.0
    ensure_classical_sieve(res.max_d)

    # 1. Prefactor
    pref_lm = evaluate_classical_log(res.pref_sq)
    v_pref = exp(pref_lm / 2)
    
    # 2. Base Term
    lm_curr = evaluate_classical_log(res.m_min)
    sign_curr = res.m_min.sign
    sum_val = sign_curr * exp(lm_curr)
    
    # 3. Ratio Summation
    for r in res.ratios
        lm_curr += evaluate_classical_log(r)
        sign_curr *= r.sign
        sum_val += sign_curr * exp(lm_curr)
    end
    
    return v_pref * sum_val
end


export qracah6j_classical, qracah3j_classical, evaluate_classical




function qracah6j_classical(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res = q6j_cyclo(j1, j2, j3, j4, j5, j6)
    return evaluate_classical(res)
end

function qracah3j_classical(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    return return evaluate_classical(res)
end

# ==============================================================================
# Exact Rational Evaluators (GCD-Free Two-Pass Algorithm)
# ==============================================================================

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

# # The Core GCD-Free Summation Engine (GMP Optimized)
# function _sum_hypergeometric_exact(res::CycloResult)
#     max_d = res.max_d
#     curr_exps = zeros(Int, max_d)
#     min_exps = zeros(Int, max_d)
    
#     # PASS 1: Track Prime Exponent Valleys to find the Global Denominator
#     @inbounds for (d, e) in res.m_min.exps
#         p = CLASSICAL_SIEVE[d]
#         if p > 1
#             curr_exps[p] += e
#             min_exps[p] = curr_exps[p]
#         end
#     end
    
#     for r in res.ratios
#         @inbounds for (d, e) in r.exps
#             p = CLASSICAL_SIEVE[d]
#             if p > 1
#                 curr_exps[p] += e
#                 if curr_exps[p] < min_exps[p]
#                     min_exps[p] = curr_exps[p]
#                 end
#             end
#         end
#     end
    
#     # PASS 2: Reconstruct N_0 and D_global
#     fill!(curr_exps, 0)
#     @inbounds for (d, e) in res.m_min.exps
#         p = CLASSICAL_SIEVE[d]
#         p > 1 && (curr_exps[p] += e)
#     end
    
#     D_global = one(BigInt)
#     N_0 = one(BigInt)
    
#     # Pre-cache BigInt primes to prevent repeated BigInt(p) allocations
#     cached_primes = Dict{Int, BigInt}()
#     _get_p(p::Int) = get!(cached_primes, p) do; BigInt(p); end
    
#     @inbounds for p in 2:max_d
#         E_p = min_exps[p]
#         K_p = E_p < 0 ? -E_p : 0 
        
#         if K_p > 0
#             Base.GMP.MPZ.mul!(D_global, D_global, _get_p(p)^K_p)
#         end
        
#         pow = curr_exps[p] + K_p
#         if pow > 0
#             Base.GMP.MPZ.mul!(N_0, N_0, _get_p(p)^pow)
#         end
#     end
    
#     # PASS 3: The Hot Loop (In-Place GMP BigInt Arithmetic, ZERO Allocations!)
#     Sum_N = BigInt()
#     Base.GMP.MPZ.set!(Sum_N, N_0)
#     res.m_min.sign < 0 && Base.GMP.MPZ.neg!(Sum_N, Sum_N)
    
#     curr_N = BigInt()
#     Base.GMP.MPZ.set!(curr_N, N_0)
#     curr_sign = res.m_min.sign
    
#     # Workspace variables to prevent loop allocations
#     r_num = BigInt()
#     r_den = BigInt()
    
#     for r in res.ratios
#         Base.GMP.MPZ.set_si!(r_num, 1)
#         Base.GMP.MPZ.set_si!(r_den, 1)
        
#         @inbounds for (d, e) in r.exps
#             p = CLASSICAL_SIEVE[d]
#             if p > 1
#                 if e > 0
#                     Base.GMP.MPZ.mul!(r_num, r_num, _get_p(p)^e)
#                 elseif e < 0
#                     Base.GMP.MPZ.mul!(r_den, r_den, _get_p(p)^(-e))
#                 end
#             end
#         end
        
#         # curr_N = (curr_N * r_num) / r_den (in-place)
#         Base.GMP.MPZ.mul!(curr_N, curr_N, r_num)
#         Base.GMP.MPZ.tdiv_q!(curr_N, curr_N, r_den) # truncating exact division
        
#         curr_sign *= r.sign
        
#         if curr_sign > 0
#             Base.GMP.MPZ.add!(Sum_N, Sum_N, curr_N)
#         else
#             Base.GMP.MPZ.sub!(Sum_N, Sum_N, curr_N)
#         end
#     end
    
#     # Only ONE final fraction allocation
#     return Sum_N // D_global
# end

# The Core GCD-Free Summation Engine (Platform-Safe, ZERO-Allocation Hot Loop)
function _sum_hypergeometric_exact(res::CycloResult)
    max_d = res.max_d
    curr_exps = zeros(Int, max_d)
    min_exps = zeros(Int, max_d)
    
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
    fill!(curr_exps, 0)
    @inbounds for (d, e) in res.m_min.exps
        p = CLASSICAL_SIEVE[d]
        p > 1 && (curr_exps[p] += e)
    end
    
    D_global = one(BigInt)
    N_0 = one(BigInt)
    
    # We can afford BigInt allocations here because this loop executes exactly ONCE 
    # per 6j-symbol, rendering GC overhead practically non-existent.
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
    
    # Workspace variables for 100% ZERO allocation in the hot loop
    r_num = BigInt()
    r_den = BigInt()
    temp_val = BigInt()
    
    for r in res.ratios
        # Reset workspaces using native Julia Int
        Base.GMP.MPZ.set_si!(r_num, 1)
        Base.GMP.MPZ.set_si!(r_den, 1)
        
        @inbounds for (d, e) in r.exps
            p = CLASSICAL_SIEVE[d]
            if p > 1
                # For a single step ratio, p^|e| is very small. Native Int is perfectly safe.
                val = p^abs(e)
                Base.GMP.MPZ.set_si!(temp_val, val)
                
                if e > 0
                    Base.GMP.MPZ.mul!(r_num, r_num, temp_val)
                elseif e < 0
                    Base.GMP.MPZ.mul!(r_den, r_den, temp_val)
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
    
    # Only ONE final fraction allocation at the very end
    return Sum_N // D_global
end









# ==============================================================================
# Fast Float64 Evaluators (NaN-Safe Log-Sum-Exp)
# ==============================================================================

@inline function evaluate_classical_log(m::CycloMonomial)
    m.sign == 0 && return -Inf
    log_val = 0.0
    @inbounds for (d, e) in m.exps
        log_val += e * CLASSICAL_LOG[d] # O(1) Cache retrieval
    end
    return log_val
end

"""
    evaluate_classical(res::CycloResult)

Evaluates a pre-constructed CycloResult in the classical Ponzano-Regge limit.
Uses a two-pass Log-Sum-Exp (LSE) shift to guarantee immunity against `NaN` overflow.
Note: Susceptible to Float64 catastrophic cancellation at j > 80.
"""
function evaluate_classical(res::CycloResult)
    res.pref_sq.sign == 0 && return 0.0
    ensure_classical_sieve(res.max_d)

    # PASS 1: Find the global maximum log-magnitude to prevent Inf
    max_lm = -Inf
    lm_curr = evaluate_classical_log(res.m_min)
    max_lm = max(max_lm, lm_curr)
    
    for r in res.ratios
        lm_curr += evaluate_classical_log(r)
        max_lm = max(max_lm, lm_curr)
    end

    # PASS 2: Shifted Summation
    pref_lm = evaluate_classical_log(res.pref_sq)
    v_pref = exp(pref_lm / 2)
    
    lm_curr = evaluate_classical_log(res.m_min)
    sign_curr = res.m_min.sign
    
    # We subtract max_lm from the exponent. The max exp is now exactly exp(0) = 1.0.
    sum_val = sign_curr * exp(lm_curr - max_lm)
    
    for r in res.ratios
        lm_curr += evaluate_classical_log(r)
        sign_curr *= r.sign
        sum_val += sign_curr * exp(lm_curr - max_lm)
    end
    
    # Multiply the shifted max back in at the very end
    return v_pref * sum_val * exp(max_lm)
end

function qracah6j_classical(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res = q6j_cyclo(j1, j2, j3, j4, j5, j6)
    return evaluate_classical(res)
end

function qracah3j_classical(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    res = q3j_cyclo(j1, j2, j3, m1, m2, m3) # Assuming you have this constructor
    return evaluate_classical(res)
end


# ==============================================================================
# Exact Rational Evaluators (The WignerSymbols Killer)
# ==============================================================================

"""
    evaluate_classical_exact(res::CycloResult)

Evaluates the exact classical limit using the zero-allocation GMP GCD-bypass.
Returns a `ClassicalResult` containing the exact squared rational and parity sign.
"""
function evaluate_classical_exact(res::CycloResult)
    res.pref_sq.sign == 0 && return ClassicalResult(0, 0//1)
    ensure_classical_sieve(res.max_d)

    # 1. Exact fraction of the triangle prefactor squared
    pref_sq_val = evaluate_to_rational(res.pref_sq)
    
    # 2. Exact integer/rational evaluation of the hypergeometric series
    sum_val = _sum_hypergeometric_exact(res)
    
    # The actual exact 6j symbol squared is pref_sq * sum^2
    final_sign = sign(sum_val) * res.pref_sq.sign
    sq_val = pref_sq_val * sum_val^2
    
    return ClassicalResult(final_sign, sq_val)
end

function qracah6j_classical_exact(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res = q6j_cyclo(j1, j2, j3, j4, j5, j6)
    return evaluate_classical_exact(res)
end


export qracah6j_classical_exact


export qracah6j_classical_exact, qracah3j_classical_exact

function qracah6j_classical_exact(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res = q6j_cyclo(j1, j2, j3, j4, j5, j6)
    (res.pref_sq.sign == 0 || res.m_min.sign == 0) && return ClassicalResult(0, 0//1)
    ensure_classical_sieve(res.max_d)
    
    pref_sq_rat = evaluate_to_rational(res.pref_sq)
    val_sum = _sum_hypergeometric_exact(res)
    
    final_sign = sign(val_sum)
    final_sq_val = pref_sq_rat * (val_sum^2)
    
    return ClassicalResult(final_sign, final_sq_val)
end

function qracah3j_classical_exact(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    res = q3j_cyclo(j1, j2, j3, m1, m2, m3)
    (res.pref_sq.sign == 0 || res.m_min.sign == 0) && return ClassicalResult(0, 0//1)
    ensure_classical_sieve(res.max_d)
    
    pref_sq_rat = evaluate_to_rational(res.pref_sq)
    val_sum = _sum_hypergeometric_exact(res)
    
    final_sign = sign(val_sum)
    final_sq_val = pref_sq_rat * (val_sum^2)
    
    return ClassicalResult(final_sign, final_sq_val)
end













# The Core GCD-Free Summation Engine
# function _sum_hypergeometric_exact(res::CycloResult)
#     max_d = res.max_d
#     curr_exps = zeros(Int, max_d)
#     min_exps = zeros(Int, max_d)
    
#     # PASS 1: Track Prime Exponent Valleys to find the Global Denominator
#     @inbounds for (d, e) in res.m_min.exps
#         p = CLASSICAL_SIEVE[d]
#         if p > 1
#             curr_exps[p] += e
#             min_exps[p] = curr_exps[p]
#         end
#     end
    
#     for r in res.ratios
#         @inbounds for (d, e) in r.exps
#             p = CLASSICAL_SIEVE[d]
#             if p > 1
#                 curr_exps[p] += e
#                 if curr_exps[p] < min_exps[p]
#                     min_exps[p] = curr_exps[p]
#                 end
#             end
#         end
#     end
    
#     # PASS 2: Reconstruct N_0 and D_global
#     fill!(curr_exps, 0)
#     @inbounds for (d, e) in res.m_min.exps
#         p = CLASSICAL_SIEVE[d]
#         p > 1 && (curr_exps[p] += e)
#     end
    
#     D_global = one(BigInt)
#     N_0 = one(BigInt)
    
#     @inbounds for p in 2:max_d
#         E_p = min_exps[p]
#         K_p = E_p < 0 ? -E_p : 0 # K_p is the max negative depth of prime p
        
#         if K_p > 0
#             D_global *= BigInt(p)^K_p
#         end
        
#         pow = curr_exps[p] + K_p
#         if pow > 0
#             N_0 *= BigInt(p)^pow
#         end
#     end
    
#     # PASS 3: The Hot Loop (Pure Integer Arithmetic, ZERO GCDs!)
#     Sum_N = res.m_min.sign > 0 ? N_0 : -N_0
#     curr_N = N_0
#     curr_sign = res.m_min.sign
    
#     for r in res.ratios
#         r_num = one(BigInt)
#         r_den = one(BigInt)
#         @inbounds for (d, e) in r.exps
#             p = CLASSICAL_SIEVE[d]
#             if p > 1
#                 if e > 0
#                     r_num *= BigInt(p)^e
#                 elseif e < 0
#                     r_den *= BigInt(p)^(-e)
#                 end
#             end
#         end
        
#         # Because we scaled by D_global, curr_N * r_num is GUARANTEED divisible by r_den
#         # We use standard div() to completely bypass fractional GCD logic
#         curr_N = div(curr_N * r_num, r_den)
#         curr_sign *= r.sign
        
#         if curr_sign > 0
#             Sum_N += curr_N
#         else
#             Sum_N -= curr_N
#         end
#     end
    
#     # Only ONE GCD is performed at the very end to finalize the fraction
#     return Sum_N // D_global
# end



# # ==============================================================================
# # Classical Limit Precomputations (Pure Julia Sieve)
# # ==============================================================================
# # - CLASSICAL_SIEVE[d] = p  if d = p^k (prime power), else 1
# # - CLASSICAL_LOG[d]   = log(p) if d = p^k, else 0.0
# const CLASSICAL_SIEVE = Int[]
# const CLASSICAL_LOG   = Float64[]

# function ensure_classical_sieve(max_d::Int)
#     curr_len = length(CLASSICAL_SIEVE)
#     if max_d > curr_len
#         new_size = max(max_d, 2 * curr_len, 5000)
        
#         resize!(CLASSICAL_SIEVE, new_size)
#         resize!(CLASSICAL_LOG, new_size)
        
#         # Default states for non-prime-powers
#         for i in (curr_len+1):new_size
#             CLASSICAL_SIEVE[i] = 1
#             CLASSICAL_LOG[i] = 0.0
#         end
#         CLASSICAL_SIEVE[1] = 0
        
#         is_prime = trues(new_size)
#         is_prime[1] = false
        
#         @inbounds for p in 2:new_size
#             if is_prime[p]
#                 log_p = log(Float64(p))
                
#                 # 1. Mark all prime powers p^k
#                 power = p
#                 while power <= new_size
#                     CLASSICAL_SIEVE[power] = p
#                     CLASSICAL_LOG[power] = log_p
#                     # Prevent integer overflow during condition check
#                     power = widen(power) * p > new_size ? new_size + 1 : power * p
#                 end
                
#                 # 2. Mark multiples as composite
#                 for mult in (p * 2):p:new_size
#                     is_prime[mult] = false
#                 end
#             end
#         end
#     end
# end

# # ==============================================================================
# # Result Struct and Base Overloads
# # ==============================================================================

# struct ClassicalResult
#     sign::Int
#     sq_val::Rational{BigInt} 
# end

# function Base.show(io::IO, res::ClassicalResult)
#     res.sign == 0 && return print(io, "0.0")
    
#     num, den = numerator(res.sq_val), denominator(res.sq_val)
#     s_num, s_den = isqrt(num), isqrt(den)
    
#     if s_num^2 == num && s_den^2 == den
#         print(io, res.sign < 0 ? "-" : "", Rational(s_num, s_den))
#     else
#         print(io, res.sign < 0 ? "-" : "", "√(", res.sq_val, ")")
#     end
# end

# # Use BigFloat internally for sqrt to prevent Float64 overflow on massive BigInts
# Base.Float64(res::ClassicalResult) = res.sign * Float64(sqrt(BigFloat(res.sq_val)))
# Base.BigFloat(res::ClassicalResult) = res.sign * sqrt(BigFloat(res.sq_val))

# function Base.Rational{T}(res::ClassicalResult) where {T <: Integer}
#     num, den = numerator(res.sq_val), denominator(res.sq_val)
#     s_num, s_den = isqrt(num), isqrt(den)
    
#     if s_num^2 == num && s_den^2 == den
#         return res.sign * Rational{T}(T(s_num), T(s_den))
#     else
#         throw(InexactError(:Rational, Rational{T}, res))
#     end
# end

# Base.:+(res::ClassicalResult, x::Number) = Float64(res) + x
# Base.:*(res::ClassicalResult, x::Number) = Float64(res) * x

# # ==============================================================================
# # Fast Float64 Evaluators (Log-Space State Tracking)
# # ==============================================================================

# """
#     evaluate_classical_log(m::CycloMonomial)
# Evaluates the log-magnitude using the precomputed Sieve Logarithms.
# O(N_sparse) time, zero allocations, no trig/log calls!
# """
# @inline function evaluate_classical_log(m::CycloMonomial)
#     m.sign == 0 && return -Inf
#     log_val = 0.0
#     @inbounds for (d, e) in m.exps
#         log_val += e * CLASSICAL_LOG[d]
#     end
#     return log_val
# end

# export qracah6j_classical, qracah3j_classical

# function qracah6j_classical(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
#     res = q6j_cyclo(j1, j2, j3, j4, j5, j6)
#     res.pref_sq.sign == 0 && return 0.0
#     ensure_classical_sieve(res.max_d)

#     # 1. Prefactor
#     pref_lm = evaluate_classical_log(res.pref_sq)
#     v_pref = exp(pref_lm / 2)
    
#     # 2. Base Term
#     lm_curr = evaluate_classical_log(res.m_min)
#     sign_curr = res.m_min.sign
#     sum_val = sign_curr * exp(lm_curr)
    
#     # 3. Ratio Summation (Log-Space prevents under/overflow!)
#     for r in res.ratios
#         lm_curr += evaluate_classical_log(r)
#         sign_curr *= r.sign
#         sum_val += sign_curr * exp(lm_curr)
#     end
    
#     return v_pref * sum_val
# end

# function qracah3j_classical(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
#     res = q3j_cyclo(j1, j2, j3, m1, m2, m3)
#     res.pref_sq.sign == 0 && return 0.0
#     ensure_classical_sieve(res.max_d)

#     pref_lm = evaluate_classical_log(res.pref_sq)
#     v_pref = exp(pref_lm / 2)
    
#     lm_curr = evaluate_classical_log(res.m_min)
#     sign_curr = res.m_min.sign
#     sum_val = sign_curr * exp(lm_curr)
    
#     for r in res.ratios
#         lm_curr += evaluate_classical_log(r)
#         sign_curr *= r.sign
#         sum_val += sign_curr * exp(lm_curr)
#     end
    
#     return v_pref * sum_val
# end

# # ==============================================================================
# # Exact Rational Evaluators
# # ==============================================================================

# """
#     evaluate_to_rational(m::CycloMonomial)
# Converts a sparse CycloMonomial to an exact Rational{BigInt} at q=1.
# """
# function evaluate_to_rational(m::CycloMonomial)
#     m.sign == 0 && return zero(Rational{BigInt})
#     num = one(BigInt)
#     den = one(BigInt)
    
#     @inbounds for (d, e) in m.exps
#         p = CLASSICAL_SIEVE[d]
#         if p > 1
#             if e > 0
#                 num *= BigInt(p)^e
#             elseif e < 0
#                 den *= BigInt(p)^abs(e)
#             end
#         end
#     end
    
#     val = num // den
#     return m.sign > 0 ? val : -val
# end

# export qracah6j_classical_exact, qracah3j_classical_exact

# function qracah6j_classical_exact(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
#     res = q6j_cyclo(j1, j2, j3, j4, j5, j6)
#     (res.pref_sq.sign == 0 || res.m_min.sign == 0) && return ClassicalResult(0, 0//1)
#     ensure_classical_sieve(res.max_d)
    
#     pref_sq_rat = evaluate_to_rational(res.pref_sq)
#     m_min_rat   = evaluate_to_rational(res.m_min)
    
#     sum_val = one(Rational{BigInt})
#     curr_ratio_prod = one(Rational{BigInt})
    
#     # Recursive Ratio Multiplication (Fastest Path)
#     for r in res.ratios
#         curr_ratio_prod *= evaluate_to_rational(r)
#         sum_val += curr_ratio_prod
#     end
    
#     total_rational_sum = m_min_rat * sum_val
#     final_sign = sign(total_rational_sum)
#     final_sq_val = pref_sq_rat * (total_rational_sum^2)
    
#     return ClassicalResult(final_sign, final_sq_val)
# end

# function qracah3j_classical_exact(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
#     res = q3j_cyclo(j1, j2, j3, m1, m2, m3)
#     (res.pref_sq.sign == 0 || res.m_min.sign == 0) && return ClassicalResult(0, 0//1)
#     ensure_classical_sieve(res.max_d)
    
#     pref_sq_rat = evaluate_to_rational(res.pref_sq)
#     m_min_rat   = evaluate_to_rational(res.m_min)
    
#     sum_val = one(Rational{BigInt})
#     curr_ratio_prod = one(Rational{BigInt})
    
#     for r in res.ratios
#         curr_ratio_prod *= evaluate_to_rational(r)
#         sum_val += curr_ratio_prod
#     end
    
#     total_rational_sum = m_min_rat * sum_val
#     final_sign = sign(total_rational_sum)
#     final_sq_val = pref_sq_rat * (total_rational_sum^2)
    
#     return ClassicalResult(final_sign, final_sq_val)
# end






# # ==============================================================================
# # Classical Limit Precomputations (Pure Julia Sieve)
# # ==============================================================================
# # - CLASSICAL_SIEVE[d] = p  if d = p^k (prime power)
# # - CLASSICAL_SIEVE[d] = 1  otherwise
# const CLASSICAL_SIEVE = Int[]

# function ensure_classical_sieve(max_d::Int)
#     if max_d > length(CLASSICAL_SIEVE)
#         new_size = max(max_d, 5000)
#         resize!(CLASSICAL_SIEVE, new_size)
#         fill!(CLASSICAL_SIEVE, 1) # Default to 1 (not a prime power)
#         CLASSICAL_SIEVE[1] = 0
        
#         is_prime = trues(new_size)
#         is_prime[1] = false
        
#         @inbounds for p in 2:new_size
#             if is_prime[p]
#                 # 1. Mark all prime powers p^k
#                 power = p
#                 while power <= new_size
#                     CLASSICAL_SIEVE[power] = p
#                     # Prevent integer overflow during condition check
#                     power = widen(power) * p > new_size ? new_size + 1 : power * p
#                 end
                
#                 # 2. Mark multiples as composite
#                 for mult in (p * 2):p:new_size
#                     is_prime[mult] = false
#                 end
#             end
#         end
#     end
# end

# # ==============================================================================
# # Result Struct and Base Overloads
# # ==============================================================================

# struct ClassicalResult
#     sign::Int
#     sq_val::Rational{BigInt} # The exact value under the square root
# end

# function Base.show(io::IO, res::ClassicalResult)
#     res.sign == 0 && return print(io, "0.0")
    
#     num, den = numerator(res.sq_val), denominator(res.sq_val)
#     s_num, s_den = isqrt(num), isqrt(den)
    
#     if s_num^2 == num && s_den^2 == den
#         print(io, res.sign < 0 ? "-" : "", Rational(s_num, s_den))
#     else
#         print(io, res.sign < 0 ? "-" : "", "√(", res.sq_val, ")")
#     end
# end

# Base.Float64(res::ClassicalResult) = res.sign * sqrt(Float64(res.sq_val))
# Base.BigFloat(res::ClassicalResult) = res.sign * sqrt(BigFloat(res.sq_val))

# function Base.Rational{T}(res::ClassicalResult) where {T <: Integer}
#     num, den = numerator(res.sq_val), denominator(res.sq_val)
#     s_num, s_den = isqrt(num), isqrt(den)
    
#     if s_num^2 == num && s_den^2 == den
#         return res.sign * Rational{T}(T(s_num), T(s_den))
#     else
#         throw(InexactError(:Rational, Rational{T}, res))
#     end
# end

# Base.:+(res::ClassicalResult, x::Number) = Float64(res) + x
# Base.:*(res::ClassicalResult, x::Number) = Float64(res) * x

# # ==============================================================================
# # Fast Float64 Evaluators (Log-Space State Tracking)
# # ==============================================================================

# """
#     evaluate_classical_log(m::CycloMonomial)
# Evaluates the log-magnitude of a sparse cyclomonomial at q=1.
# """
# @inline function evaluate_classical_log(m::CycloMonomial)
#     m.sign == 0 && return -Inf
#     log_val = 0.0
#     @inbounds for (d, e) in m.exps
#         p = CLASSICAL_SIEVE[d]
#         if p > 1
#             log_val += e * log(Float64(p))
#         end
#     end
#     return log_val
# end

# export qracah6j_classical, qracah3j_classical

# function qracah6j_classical(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
#     res = q6j_cyclo(j1, j2, j3, j4, j5, j6)
#     res.pref_sq.sign == 0 && return 0.0
#     ensure_classical_sieve(res.max_d)

#     # 1. Prefactor
#     pref_lm = evaluate_classical_log(res.pref_sq)
#     v_pref = exp(pref_lm / 2)
    
#     # 2. Base Term
#     lm_curr = evaluate_classical_log(res.m_min)
#     sign_curr = res.m_min.sign
#     sum_val = sign_curr * exp(lm_curr)
    
#     # 3. Ratio Summation (Log-Space prevents under/overflow!)
#     for r in res.ratios
#         lm_curr += evaluate_classical_log(r)
#         sign_curr *= r.sign
#         sum_val += sign_curr * exp(lm_curr)
#     end
    
#     return v_pref * sum_val
# end

# function qracah3j_classical(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
#     res = q3j_cyclo(j1, j2, j3, m1, m2, m3)
#     res.pref_sq.sign == 0 && return 0.0
#     ensure_classical_sieve(res.max_d)

#     pref_lm = evaluate_classical_log(res.pref_sq)
#     v_pref = exp(pref_lm / 2)
    
#     lm_curr = evaluate_classical_log(res.m_min)
#     sign_curr = res.m_min.sign
#     sum_val = sign_curr * exp(lm_curr)
    
#     for r in res.ratios
#         lm_curr += evaluate_classical_log(r)
#         sign_curr *= r.sign
#         sum_val += sign_curr * exp(lm_curr)
#     end
    
#     return v_pref * sum_val
# end

# # ==============================================================================
# # Exact Rational Evaluators (Prime Exponent Tracker)
# # ==============================================================================

# # Helper to rapidly convert a prime exponent array into a BigInt fraction
# function _exps_to_rational(p_exps::Vector{Int})
#     num = one(BigInt)
#     den = one(BigInt)
#     @inbounds for p in 2:length(p_exps)
#         e = p_exps[p]
#         if e > 0
#             num *= BigInt(p)^e
#         elseif e < 0
#             den *= BigInt(p)^abs(e)
#         end
#     end
#     return num // den
# end

# export qracah6j_classical_exact, qracah3j_classical_exact

# function qracah6j_classical_exact(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
#     res = q6j_cyclo(j1, j2, j3, j4, j5, j6)
#     (res.pref_sq.sign == 0 || res.m_min.sign == 0) && return ClassicalResult(0, 0//1)
#     ensure_classical_sieve(res.max_d)
    
#     # 1. Prefactor
#     p_exps_pref = zeros(Int, res.max_d)
#     @inbounds for (d, e) in res.pref_sq.exps
#         p = CLASSICAL_SIEVE[d]; p > 1 && (p_exps_pref[p] += e)
#     end
#     pref_sq_rat = _exps_to_rational(p_exps_pref)
    
#     # 2. Base Term (Running Prime Exponent Tracker)
#     p_exps_curr = zeros(Int, res.max_d)
#     @inbounds for (d, e) in res.m_min.exps
#         p = CLASSICAL_SIEVE[d]; p > 1 && (p_exps_curr[p] += e)
#     end
    
#     sign_curr = res.m_min.sign
#     sum_val = sign_curr * _exps_to_rational(p_exps_curr)
    
#     # 3. Ratio Summation (Pure integer addition!)
#     for r in res.ratios
#         sign_curr *= r.sign
#         @inbounds for (d, e) in r.exps
#             p = CLASSICAL_SIEVE[d]; p > 1 && (p_exps_curr[p] += e)
#         end
#         sum_val += sign_curr * _exps_to_rational(p_exps_curr)
#     end
    
#     # Final Result: sqrt(pref * sum^2)
#     final_sign = sign(sum_val)
#     final_sq_val = pref_sq_rat * (sum_val^2)
    
#     return ClassicalResult(final_sign, final_sq_val)
# end

# function qracah3j_classical_exact(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
#     res = q3j_cyclo(j1, j2, j3, m1, m2, m3)
#     (res.pref_sq.sign == 0 || res.m_min.sign == 0) && return ClassicalResult(0, 0//1)
#     ensure_classical_sieve(res.max_d)
    
#     p_exps_pref = zeros(Int, res.max_d)
#     @inbounds for (d, e) in res.pref_sq.exps
#         p = CLASSICAL_SIEVE[d]; p > 1 && (p_exps_pref[p] += e)
#     end
#     pref_sq_rat = _exps_to_rational(p_exps_pref)
    
#     p_exps_curr = zeros(Int, res.max_d)
#     @inbounds for (d, e) in res.m_min.exps
#         p = CLASSICAL_SIEVE[d]; p > 1 && (p_exps_curr[p] += e)
#     end
    
#     sign_curr = res.m_min.sign
#     sum_val = sign_curr * _exps_to_rational(p_exps_curr)
    
#     for r in res.ratios
#         sign_curr *= r.sign
#         @inbounds for (d, e) in r.exps
#             p = CLASSICAL_SIEVE[d]; p > 1 && (p_exps_curr[p] += e)
#         end
#         sum_val += sign_curr * _exps_to_rational(p_exps_curr)
#     end
    
#     final_sign = sign(sum_val)
#     final_sq_val = pref_sq_rat * (sum_val^2)
    
#     return ClassicalResult(final_sign, final_sq_val)
# end