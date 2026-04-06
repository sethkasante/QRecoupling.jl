
# ----------------------------------------------------------------
#           --- Projection to classical  ---- 
# Classical Limit (q -> 1): Numerical & Exact Projectors
# Computes the Ponzano-Regge classical limit where level k → ∞
# ----------------------------------------------------------------

using Base.GMP.MPZ

# --- Prime Power Sieve (Φ_d(1) Logic) ---
const CLASSICAL_SIEVE = Int[]
const CLASSICAL_LOG   = Float64[]

"""
    ensure_classical_sieve(max_d::Int)
Populates the sieve where Φ_d(1) = p if d = pᵏ, and 1 otherwise.
"""
function ensure_classical_sieve(max_d::Int)
    curr_len = length(CLASSICAL_SIEVE)
    max_d <= curr_len && return
    
    new_size = max(max_d, 2 * curr_len, 8192)
    resize!(CLASSICAL_SIEVE, new_size)
    resize!(CLASSICAL_LOG, new_size)
    
    fill!(view(CLASSICAL_SIEVE, curr_len+1:new_size), 1)
    fill!(view(CLASSICAL_LOG, curr_len+1:new_size), 0.0)
    
    is_prime = trues(new_size)
    @inbounds for p in 2:new_size
        if is_prime[p]
            lp = log(Float64(p))
            power = p
            while power <= new_size
                CLASSICAL_SIEVE[power] = p
                CLASSICAL_LOG[power] = lp
                power > (new_size ÷ p) ? break : (power *= p)
            end
            for mult in (p * 2):p:new_size; is_prime[mult] = false; end
        end
    end
end

# --- Exact Rational Workspaces ---
const CLASSICAL_WS_CURR = Int[]
const CLASSICAL_WS_MIN  = Int[]

function _ensure_workspaces(n::Int)
    if length(CLASSICAL_WS_CURR) < n
        resize!(CLASSICAL_WS_CURR, max(n, 8192))
        resize!(CLASSICAL_WS_MIN, max(n, 8192))
    end
end

# ---------------------------------------------------------
# NUMERICAL PROJECTORS (Zero-Allocation Streaming LSE)
# ---------------------------------------------------------

"""
    project_classical(m::CyclotomicMonomial, ::Type{T}=Float64)
Evaluates a monomial as q -> 1. Returns a real value of type T.
"""
function project_classical(m::CyclotomicMonomial, ::Type{T}=Float64) where T
    m.sign == 0 && return zero(T)
    ensure_classical_sieve(m.max_d)
    lmag = CLASSICAL_LOG
    
    lm = zero(T)
    pe = m.phi_exps
    @inbounds for j in eachindex(pe)
        d, e = pe[j]
        lm += e * lmag[d]
    end
    
    (isnan(lm) || lm == -T(Inf)) && return zero(T)
    return m.sign * exp(lm)
end

"""
    project_classical(dcr::DCR, ::Type{T}=Float64)
Evaluates the full DCR series as q -> 1. 
Uses a single-pass Streaming Log-Sum-Exp algorithm to prevent overflow.
"""
function project_classical(dcr::DCR, ::Type{T}=Float64) where T
    ensure_classical_sieve(dcr.max_d)
    lmag = CLASSICAL_LOG
    
    @inline _eval_log(m) = begin
        m.sign == 0 && return (-T(Inf), 0)
        lm = zero(T)
        pe = m.phi_exps
        @inbounds for j in eachindex(pe)
            d, e = pe[j]; lm += e * lmag[d]
        end
        return (lm, m.sign)
    end

    (lr, sr), (lrad, srad), (lb, sb) = _eval_log(dcr.root), _eval_log(dcr.radical), _eval_log(dcr.base)
    
    curr_l = lr + (0.5 * lrad) + lb
    curr_s = sr * sb 
    
    (isinf(curr_l) || curr_s == 0 || isnan(curr_l)) && return zero(T)

    max_l = curr_l
    fsum = T(curr_s)

    @inbounds for i in 1:length(dcr.ratios)
        (rl, rs) = _eval_log(dcr.ratios[i])
        curr_l += rl
        curr_s *= rs
        
        if curr_l > max_l
            fsum = fsum * exp(max_l - curr_l) + curr_s
            max_l = curr_l
        else
            fsum += curr_s * exp(curr_l - max_l)
        end
    end

    final_val = exp(max_l) * fsum
    return isnan(final_val) ? zero(T) : final_val
end


# ---------------------------------------------------------
# EXACT PROJECTORS (Rational BigInt / Zero-GCD)
# ---------------------------------------------------------

function _fast_rat_inplace(m::CyclotomicMonomial, sieve::Vector{Int}, buf_big::BigInt)
    m.sign == 0 && return 0//1
    num, den = one(BigInt), one(BigInt)
    pe = m.phi_exps
    @inbounds for j in eachindex(pe)
        d, e = pe[j]
        p = sieve[d]
        p <= 1 && continue
        
        ccall((:__gmpz_ui_pow_ui, :libgmp), Cvoid, (Ref{BigInt}, Culong, Culong), buf_big, Culong(p), Culong(abs(e)))
        e > 0 ? MPZ.mul!(num, num, buf_big) : MPZ.mul!(den, den, buf_big)
    end
    return m.sign * (num // den)
end

"""
    project_classical_exact(m::CyclotomicMonomial)
Returns σ * √(val²) exactly. Result is a ClassicalResult.
"""
function project_classical_exact(m::CyclotomicMonomial)
    m.sign == 0 && return ClassicalResult(0, 0//1)
    ensure_classical_sieve(m.max_d)
    val_big = BigInt()
    rat = _fast_rat_inplace(m, CLASSICAL_SIEVE, val_big)
    return ClassicalResult(sign(rat), rat^2)
end

"""
    project_classical_exact(dcr::DCR)
Evaluates the DCR series exactly as a Rational BigInt (Zero-GCD).
"""
function project_classical_exact(dcr::DCR)
    ensure_classical_sieve(dcr.max_d)
    _ensure_workspaces(dcr.max_d)
    
    sieve = CLASSICAL_SIEVE
    curr_exps = CLASSICAL_WS_CURR
    min_exps  = CLASSICAL_WS_MIN
    
    @inbounds fill!(view(curr_exps, 1:dcr.max_d), 0)
    @inbounds fill!(view(min_exps, 1:dcr.max_d), 0)

    # 1. Pass 1: Valley Tracking (Denominator Identification)
    for (d, e) in dcr.base.phi_exps
        p = sieve[d]; p <= 1 && continue
        curr_exps[p] += e; min_exps[p] = min(0, curr_exps[p])
    end
    for r in dcr.ratios
        for (d, e) in r.phi_exps
            p = sieve[d]; p <= 1 && continue
            curr_exps[p] += e
            (curr_exps[p] < min_exps[p]) && (min_exps[p] = curr_exps[p])
        end
    end

    val_big = BigInt()
    D_glob, N_0 = BigInt(1), BigInt(1)

    # 2. Pass 2: Reconstruct N_0 and D_global
    @inbounds fill!(view(curr_exps, 1:dcr.max_d), 0)
    for (d, e) in dcr.base.phi_exps
        p = sieve[d]; p > 1 && (curr_exps[p] += e)
    end

    @inbounds for p in 2:dcr.max_d
        K_p = min_exps[p] < 0 ? -min_exps[p] : 0 
        if K_p > 0
            ccall((:__gmpz_ui_pow_ui, :libgmp), Cvoid, (Ref{BigInt}, Culong, Culong), val_big, Culong(p), Culong(K_p))
            MPZ.mul!(D_glob, D_glob, val_big)
        end
        pow = curr_exps[p] + K_p
        if pow > 0
            ccall((:__gmpz_ui_pow_ui, :libgmp), Cvoid, (Ref{BigInt}, Culong, Culong), val_big, Culong(p), Culong(pow))
            MPZ.mul!(N_0, N_0, val_big)
        end
    end

    # 3. Pass 3: Integer Sum Hot Loop
    Sum_N = BigInt(); MPZ.set!(Sum_N, N_0)
    (dcr.base.sign < 0) && MPZ.neg!(Sum_N, Sum_N)
    curr_N = BigInt(); MPZ.set!(curr_N, N_0)
    curr_sign = dcr.base.sign
    
    r_num, r_den = BigInt(), BigInt()
    for r in dcr.ratios
        MPZ.set_si!(r_num, 1); MPZ.set_si!(r_den, 1)
        pe = r.phi_exps
        @inbounds for j in eachindex(pe)
            d, e = pe[j]; p = sieve[d]; p <= 1 && continue
            ccall((:__gmpz_ui_pow_ui, :libgmp), Cvoid, (Ref{BigInt}, Culong, Culong), val_big, Culong(p), Culong(abs(e)))
            e > 0 ? MPZ.mul!(r_num, r_num, val_big) : MPZ.mul!(r_den, r_den, val_big)
        end
        MPZ.mul!(curr_N, curr_N, r_num); MPZ.tdiv_q!(curr_N, curr_N, r_den)
        curr_sign *= r.sign
        curr_sign > 0 ? MPZ.add!(Sum_N, Sum_N, curr_N) : MPZ.sub!(Sum_N, Sum_N, curr_N)
    end

    # 4. Pass 4: Final Prefactor Fusion
    root_rat = _fast_rat_inplace(dcr.root, sieve, val_big)
    rad_rat  = _fast_rat_inplace(dcr.radical, sieve, val_big)
    
    total_sum_rat = (Sum_N // D_glob) * root_rat
    return ClassicalResult(sign(total_sum_rat) * dcr.radical.sign, rad_rat * (total_sum_rat^2))
end

