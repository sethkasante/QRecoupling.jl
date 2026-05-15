

# ------------------------------------------------------------------------
#           --- Projection to Classical Limit ----
#  Evaluates unprojected CyclotomicMonomial and DCR at q = 1
#  Classical Limit (q → 1): Numerical & Exact Projectors
#  Computes the Ponzano-Regge classical limit where level k → ∞
#
#  Thread-safe model:
#    - CLASSICAL_SIEVE / CLASSICAL_LOG : grow-only after init, one lock
# ------------------------------------------------------------------------

using Base.GMP.MPZ

# --------------------------------
#  Prime-Power Sieve   
# Φ_d(1) = p if d = pᵏ, else 1
# --------------------------------

const CLASSICAL_SIEVE = Int[]
const CLASSICAL_LOG = Float64[]
const CLASSICAL_SIEVE_LOCK = ReentrantLock()

"""
    ensure_classical_sieve(max_d::Int)

Grow the shared prime-power sieve to cover indices `1:max_d`.

**Threading note:** The outer length check is a lock-free read on a
grow-only `Vector`.  Because we never shrink `CLASSICAL_SIEVE`, a stale
read of a too-small length is the only possible race outcome, and that
causes a redundant (but correct) lock acquisition.  The authoritative
check is always the one inside the lock.
"""
function ensure_classical_sieve(max_d::Int)
    # Fast path — safe: grow-only vector, stale read -> redundant lock at worst
    length(CLASSICAL_SIEVE) >= max_d && return

    @lock CLASSICAL_SIEVE_LOCK begin
        length(CLASSICAL_SIEVE) >= max_d && return   # authoritative check

        old_size = length(CLASSICAL_SIEVE)
        new_size = max(max_d, 2 * old_size, 8192)

        resize!(CLASSICAL_SIEVE, new_size)
        resize!(CLASSICAL_LOG,   new_size)

        # Initialise new region: Φ_d(1) = 1 for non-prime-powers
        fill!(view(CLASSICAL_SIEVE, old_size+1:new_size), 1)
        fill!(view(CLASSICAL_LOG,   old_size+1:new_size), 0.0)

        # Sieve of Eratosthenes — mark prime powers only
        is_prime = trues(new_size)
        @inbounds for p in 2:new_size
            is_prime[p] || continue
            lp = log(Float64(p))

            # All powers pᵏ ≤ new_size
            power = p
            while power <= new_size
                CLASSICAL_SIEVE[power] = p
                CLASSICAL_LOG[power]   = lp
                # Overflow-safe: power ≤ new_size < typemax(Int)/p guaranteed
                power * p > new_size && break
                power *= p
            end

            # Strike composites
            for mult in (2p):p:new_size
                is_prime[mult] = false
            end
        end
    end
    return
end


# ----------------------------------
#  Internal log-magnitude helpers
# ----------------------------------

"""
    _log_mag_classical(m::CyclotomicMonomial, lmag::Vector{Float64}, ::Type{T})

Returns `(log_magnitude::T, sign::Int)` for monomial `m` in the q→1 limit.

`lmag` must be the `CLASSICAL_LOG` vector, passed by the caller after
`ensure_classical_sieve` has been called — this avoids repeated global
lookups in hot loops.
"""
@inline function _log_mag_classical(m::CyclotomicMonomial,lmag::Vector{Float64},::Type{T}) where T
    m.sign == 0 && return (-T(Inf), 0)

    lm = zero(T)
    @inbounds for j in eachindex(m.phi_exps)
        d, e = m.phi_exps[j]
        lm += T(e) * T(lmag[d])
    end
    return (lm, m.sign)
end

# --------------------------------------------------
#  Numerical Projectors  (Zero-Allocation LSE)
# --------------------------------------------------

"""
    project_classical(m::CyclotomicMonomial, ::Type{T}=Float64) -> T

Evaluate monomial `m` at q -> 1.  Returns a real value of type `T`.
"""
function project_classical(m::CyclotomicMonomial, ::Type{T}=Float64) where T
    m.sign == 0 && return zero(T)

    ensure_classical_sieve(m.max_d)
    # lmag alias is safe: sieve is already large enough, no resize can follow
    lmag = CLASSICAL_LOG

    lm = zero(T)
    @inbounds for j in eachindex(m.phi_exps)
        d, e = m.phi_exps[j]
        lm += T(e) * T(lmag[d])
    end

    isinf(lm) && lm < 0 && return zero(T)
    # NaN = 0*Inf in the exponents 
    @assert !isnan(lm) "NaN in log-magnitude: degenerate monomial (sign=$(m.sign))"

    return m.sign * exp(lm)
end


"""
    project_classical(dcr::DCR, ::Type{T}=Float64) -> T

Evaluate the full DCR series at q -> 1 using a single-pass Streaming
Log-Sum-Exp algorithm to prevent overflow/underflow.
"""
function project_classical(dcr::DCR, ::Type{T}=Float64) where T
    dcr.base.sign == 0 && return zero(T)

    ensure_classical_sieve(dcr.max_d)
    lmag = CLASSICAL_LOG   # alias safe: sieve is fully initialised above

    (lr,  sr) = _log_mag_classical(dcr.root, lmag, T)
    (lrad, _) = _log_mag_classical(dcr.radical, lmag, T)
    (lb,  sb) = _log_mag_classical(dcr.base, lmag, T)

    curr_l = lr + T(0.5) * lrad + lb
    curr_s = sr * sb

    (isinf(curr_l) || curr_s == 0) && return zero(T)
    @assert !isnan(curr_l) "NaN in prefactor log-magnitude: degenerate DCR"

    max_l = curr_l
    fsum  = T(curr_s)

    @inbounds for i in eachindex(dcr.ratios)
        (rl, rs) = _log_mag_classical(dcr.ratios[i], lmag, T)
        curr_l  += rl
        curr_s  *= rs

        if curr_l > max_l
            fsum   = fsum * exp(max_l - curr_l) + T(curr_s)
            max_l  = curr_l
        else
            fsum  += T(curr_s) * exp(curr_l - max_l)
        end
    end

    final_val = exp(max_l) * fsum
    return isnan(final_val) ? zero(T) : final_val
end


# -------------------------------------------------
#  Exact Result (Rational BigInt / Zero-GCD)
# -------------------------------------------------

"""
    _fast_rat_inplace(m, sieve, buf) → Rational{BigInt}

Compute `m` as an exact rational in the q->1 limit using the prime-power
sieve. `buf` is a reusable `BigInt` scratch buffer (mutated, not returned).
"""
function _fast_rat_inplace(m::CyclotomicMonomial,sieve::Vector{Int},buf::BigInt)
    m.sign == 0 && return 0 // one(BigInt)

    num = one(BigInt)
    den = one(BigInt)

    @inbounds for j in eachindex(m.phi_exps)
        d, e = m.phi_exps[j]
        p    = sieve[d]
        p <= 1 && continue  # Φ_d(1) = 1 

        ae = Culong(abs(e))
        ccall((:__gmpz_ui_pow_ui, :libgmp), Cvoid,
              (Ref{BigInt}, Culong, Culong), buf, Culong(p), ae)

        e > 0 ? MPZ.mul!(num, num, buf) : MPZ.mul!(den, den, buf)
    end

    return m.sign * (num // den)
end


"""
    project_classical_exact(m::CyclotomicMonomial) -> ClassicalResult

Evaluate monomial `m` exactly at q → 1.  Returns a `ClassicalResult`
encodes `σ · √(val)` where `val` is a `Rational{BigInt}`.
"""
function project_classical_exact(m::CyclotomicMonomial)
    m.sign == 0 && return ClassicalResult(0, 0 // one(BigInt))

    ensure_classical_sieve(m.max_d)
    buf = BigInt()
    rat = _fast_rat_inplace(m, CLASSICAL_SIEVE, buf)
    return ClassicalResult(sign(rat), rat^2)
end


"""
    project_classical_exact(dcr::DCR) -> ClassicalResult

Evaluate the full DCR series exactly at q -> 1 as a `Rational{BigInt}`.

Algorithm: (four passes, no GCD on every ratio):
  1. Find the global denominator exponents
  2. Denominator + N_0 reconstruction
  3. Sum hot loop: accumulate numerators over `D_glob`
  4. Fuse prefactor: multiply root and radical

Function is fully thread-safe.
"""
function project_classical_exact(dcr::DCR)
    dcr.base.sign == 0 && return ClassicalResult(0, 0 // one(BigInt))

    ensure_classical_sieve(dcr.max_d)
    sieve = CLASSICAL_SIEVE

    max_p = dcr.max_d  # prime index never exceeds max_d
    curr_exps = zeros(Int, max_p)
    min_exps = zeros(Int, max_p)

    # denominator identification 
    @inbounds for (d, e) in dcr.base.phi_exps
        p = sieve[d];  p <= 1 && continue
        curr_exps[p] += e
        curr_exps[p] < min_exps[p] && (min_exps[p] = curr_exps[p])
    end

    @inbounds for r in dcr.ratios
        for (d, e) in r.phi_exps
            p = sieve[d];  p <= 1 && continue
            curr_exps[p] += e
            curr_exps[p] < min_exps[p] && (min_exps[p] = curr_exps[p])
        end
    end

    # reconstruct D_glob and N_0 
    fill!(curr_exps, 0)
    @inbounds for (d, e) in dcr.base.phi_exps
        p = sieve[d];  p > 1 && (curr_exps[p] += e)
    end

    buf    = BigInt()
    D_glob = one(BigInt)
    N_0    = one(BigInt)

    @inbounds for p in 2:max_p
        K_p = min_exps[p] < 0 ? -min_exps[p] : 0

        if K_p > 0
            ccall((:__gmpz_ui_pow_ui, :libgmp), Cvoid,
                  (Ref{BigInt}, Culong, Culong), buf, Culong(p), Culong(K_p))
            MPZ.mul!(D_glob, D_glob, buf)
        end

        pow = curr_exps[p] + K_p
        if pow > 0
            ccall((:__gmpz_ui_pow_ui, :libgmp), Cvoid,
                  (Ref{BigInt}, Culong, Culong), buf, Culong(p), Culong(pow))
            MPZ.mul!(N_0, N_0, buf)
        end
    end

    # Integer Sum: loop 
    Sum_N  = MPZ.set(N_0)
    dcr.base.sign < 0 && MPZ.neg!(Sum_N, Sum_N)

    curr_N    = MPZ.set(N_0)
    curr_sign = dcr.base.sign
    r_num     = BigInt()
    r_den     = BigInt()

    @inbounds for r in dcr.ratios
        MPZ.set_si!(r_num, 1)
        MPZ.set_si!(r_den, 1)

        for (d, e) in r.phi_exps
            p = sieve[d];  p <= 1 && continue
            ccall((:__gmpz_ui_pow_ui, :libgmp), Cvoid,
                  (Ref{BigInt}, Culong, Culong), buf, Culong(p), Culong(abs(e)))
            e > 0 ? MPZ.mul!(r_num, r_num, buf) : MPZ.mul!(r_den, r_den, buf)
        end

        MPZ.mul!(curr_N, curr_N, r_num)
        MPZ.tdiv_q!(curr_N, curr_N, r_den)
        curr_sign *= r.sign

        curr_sign > 0 ? MPZ.add!(Sum_N, Sum_N, curr_N) : MPZ.sub!(Sum_N, Sum_N, curr_N)
    end

    # fuse prefactor 
    root_rat = _fast_rat_inplace(dcr.root,     sieve, buf)
    rad_rat  = _fast_rat_inplace(dcr.radical,  sieve, buf)

    # Radical magnitude must be unsigned. Sign is tracked via dcr.radical.sign
    rad_unsigned = abs(rad_rat)

    total_sum_rat  = (Sum_N // D_glob) * root_rat
    result_sign    = sign(total_sum_rat) * dcr.radical.sign

    return ClassicalResult(result_sign, rad_unsigned * total_sum_rat^2)
end

# ---------------------------------------- 
#  QPhase Classical Projectors (q -> 1)
# ----------------------------------------
project_classical(p::QPhase, ::Type{T}=Float64) where T = T(p.sign)
project_classical_exact(p::QPhase) = ClassicalResult(sign(p), p.sign == 0 ? 0//1 : 1//1)