
#generic_engine.jl

# Compute cyclomonomial representation of Racah Wigner 3j, 6j Symbols ---------



# ---- Allocation-Free Operations ------ 

@inline function _ensure_capacity!(buf::SymbolicBuffer, n::Int)
    old_len = length(buf.exps)
    if n > old_len
        resize!(buf.exps, n)
        @inbounds for i in (old_len+1):n
            buf.exps[i] = 0
        end
    end
    return nothing
end

function add_qfact!(buf::SymbolicBuffer, n::Int, power::Int=1)
    n <= 1 && return nothing
    buf.z_pow += power * (-(n * (n - 1)) ÷ 2)
    _ensure_capacity!(buf, n)
    @inbounds for d in 2:n
        buf.exps[d] += power * div(n, d)
    end
    return nothing
end

function add_qint!(buf::SymbolicBuffer, n::Int, power::Int=1)
    n <= 1 && return nothing
    buf.z_pow += power * (1 - n)
    _ensure_capacity!(buf, n)
    @inbounds for d in 2:n
        if n % d == 0
            buf.exps[d] += power
        end
    end
    return nothing
end


#---- Single allocation setup ------- 


function qdelta2_symb!(buf::SymbolicBuffer, j1::Spin, j2::Spin, j3::Spin)
    add_qfact!(buf, Int(j1+j2-j3), 1)
    add_qfact!(buf, Int(j1-j2+j3), 1)
    add_qfact!(buf, Int(-j1+j2+j3), 1)
    add_qfact!(buf, Int(j1+j2+j3+1), -1)
    return nothing
end

function qtricoeff2_symb!(buf::SymbolicBuffer, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    add_qfact!(buf, Int(j1+j2-j3), 1); add_qfact!(buf, Int(j1-j2+j3), 1); add_qfact!(buf, Int(-j1+j2+j3), 1); add_qfact!(buf, Int(j1+j2+j3+1), -1)
    add_qfact!(buf, Int(j1+j5-j6), 1); add_qfact!(buf, Int(j1-j5+j6), 1); add_qfact!(buf, Int(-j1+j5+j6), 1); add_qfact!(buf, Int(j1+j5+j6+1), -1)
    add_qfact!(buf, Int(j2+j4-j6), 1); add_qfact!(buf, Int(j2-j4+j6), 1); add_qfact!(buf, Int(-j2+j4+j6), 1); add_qfact!(buf, Int(j2+j4+j6+1), -1)
    add_qfact!(buf, Int(j3+j4-j5), 1); add_qfact!(buf, Int(j3-j4+j5), 1); add_qfact!(buf, Int(-j3+j4+j5), 1); add_qfact!(buf, Int(j3+j4+j5+1), -1)
    return nothing
end


#----- qRacah Series (Hypergeometric Ratio Method) ----- 

function q6jseries_symb(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::Vector{CycloMonomial}
    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    S_z = CycloMonomial[]
    z_min = max(α1, α2, α3, α4)
    z_max = min(β1, β2, β3)
    
    if z_min <= z_max
        # Allocate a single buffer for the entire series evaluation
        cap = max(z_max + 1, β1 + 1, β2 + 1, β3 + 1)
        buf = SymbolicBuffer(cap)
        
        # Initialize the first term
        buf.sign = iseven(z_min) ? 1 : -1
        add_qfact!(buf, z_min+1, 1)
        add_qfact!(buf, z_min-α1, -1); add_qfact!(buf, z_min-α2, -1); add_qfact!(buf, z_min-α3, -1); add_qfact!(buf, z_min-α4, -1)
        add_qfact!(buf, β1-z_min, -1); add_qfact!(buf, β2-z_min, -1); add_qfact!(buf, β3-z_min, -1)
        
        push!(S_z, snapshot(buf))

        # Iteratively update using the ratio method
        for z in z_min+1 : z_max
            buf.sign = -buf.sign # Alternating sum
            
            add_qint!(buf, z + 1, 1)
            add_qint!(buf, β1 - z + 1, 1)
            add_qint!(buf, β2 - z + 1, 1)
            add_qint!(buf, β3 - z + 1, 1)
            
            add_qint!(buf, z - α1, -1)
            add_qint!(buf, z - α2, -1)
            add_qint!(buf, z - α3, -1)
            add_qint!(buf, z - α4, -1)
            
            push!(S_z, snapshot(buf))
        end
    end
    return S_z
end

function qracah6j_generic(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    # Estimate max factorial needed for prefactor to avoid resizes
    cap = Int(2 * max(j1+j2+j3, j1+j5+j6, j2+j4+j6, j3+j4+j5) + 2)
    buf = SymbolicBuffer(cap)
    
    qtricoeff2_symb!(buf, j1, j2, j3, j4, j5, j6)
    Tc2 = snapshot(buf)
    
    series = q6jseries_symb(j1, j2, j3, j4, j5, j6)
    return GenericResult(Tc2, series)
end




export q6j_recursive


# src/generic_engine.jl

"""
    q6j_recursive(j1, j2, j3, j4, j5, j6)
Constructs the Generic6j struct using the ratio method.
"""
function q6j_recursive(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    α = (Int(j1+j2+j3), Int(j1+j5+j6), Int(j2+j4+j6), Int(j3+j4+j5))
    β = (Int(j1+j2+j4+j5), Int(j1+j3+j4+j6), Int(j2+j3+j5+j6))
    
    z_min, z_max = max(α...), min(β...)
    cap = max(z_max + 1, β...)
    buf = SymbolicBuffer(cap)
    
    # 1. Prefactor
    qtricoeff2_symb!(buf, j1, j2, j3, j4, j5, j6)
    pref_sq = snapshot(buf)
    
    # 2. Base Term M(z_min)
    # Reset buffer
    buf.sign = iseven(z_min) ? 1 : -1; buf.z_pow = 0; fill!(buf.exps, 0)
    add_qfact!(buf, z_min+1, 1)
    for a in α; add_qfact!(buf, z_min-a, -1); end
    for b in β; add_qfact!(buf, b-z_min, -1); end
    m_min = snapshot(buf)
    
    # 3. Ratios R_z
    ratios = CycloMonomial[]
    for z in z_min : z_max-1
        # Ratio = - [z+2] * Π[β-z] / Π[z+1-α]
        buf.sign = -1; buf.z_pow = 0; fill!(buf.exps, 0)
        add_qint!(buf, z + 2, 1)
        for b in β; add_qint!(buf, b - z, 1); end
        for a in α; add_qint!(buf, z + 1 - a, -1); end
        push!(ratios, snapshot(buf))
    end
    
    return Generic6j(pref_sq, m_min, ratios, z_min:z_max)
end







# ----- quantum 3j symbols (cyclomonomial)

function q3jseries_symb(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin)
    α1 = Int(j3 - j2 + m1) 
    α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 - j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    S_z = CycloMonomial[]
    z_min = max(-α1, -α2, 0)
    z_max = min(β1, β2, β3)

    if z_min <= z_max
        # Allocate single buffer
        cap = max(z_max + max(α1, α2), β1, β2, β3) + 1
        buf = SymbolicBuffer(cap)
        
        # Initialize first term
        buf.sign = isodd(Int(z_min + α1 - α2)) ? -1 : 1
        add_qfact!(buf, z_min, -1)
        add_qfact!(buf, α1 + z_min, -1); add_qfact!(buf, α2 + z_min, -1)
        add_qfact!(buf, β1 - z_min, -1); add_qfact!(buf, β2 - z_min, -1); add_qfact!(buf, β3 - z_min, -1)
        
        push!(S_z, snapshot(buf))

        # Iteratively update
        for z in z_min+1 : z_max
            buf.sign = -buf.sign
            
            add_qint!(buf, β1 - z + 1, 1)
            add_qint!(buf, β2 - z + 1, 1)
            add_qint!(buf, β3 - z + 1, 1)
            
            add_qint!(buf, z, -1)
            add_qint!(buf, α1 + z, -1)
            add_qint!(buf, α2 + z, -1)
            
            push!(S_z, snapshot(buf))
        end
    end
    return S_z
end

function qracah3j_generic(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    cap = Int(j1 + j2 + j3 + abs(m1) + abs(m2) + 2)
    buf = SymbolicBuffer(cap)
    
    qdelta2_symb!(buf, j1, j2, j3) 
    add_qfact!(buf, Int(j1+m1), 1); add_qfact!(buf, Int(j1-m1), 1)
    add_qfact!(buf, Int(j2+m2), 1); add_qfact!(buf, Int(j2-m2), 1)
    add_qfact!(buf, Int(j3-m1-m2), 1); add_qfact!(buf, Int(j3+m1+m2), 1)
              
    pref_sq = snapshot(buf)
    series = q3jseries_symb(j1, j2, j3, m1, m2)
    return GenericResult(pref_sq, series)
end

export qracah3j_recursive

function qracah3j_recursive(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    # --- 1. Top-Level Prefactor (Squared) ---
    cap_pref = Int(j1 + j2 + j3 + abs(m1) + abs(m2) + 2)
    buf_pref = SymbolicBuffer(cap_pref)
    
    qdelta2_symb!(buf_pref, j1, j2, j3) 
    add_qfact!(buf_pref, Int(j1+m1), 1); add_qfact!(buf_pref, Int(j1-m1), 1)
    add_qfact!(buf_pref, Int(j2+m2), 1); add_qfact!(buf_pref, Int(j2-m2), 1)
    add_qfact!(buf_pref, Int(j3-m1-m2), 1); add_qfact!(buf_pref, Int(j3+m1+m2), 1)
              
    pref_sq = snapshot(buf_pref)

    # --- 2. Hypergeometric Bounds ---
    α1 = Int(j3 - j2 + m1) 
    α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 - j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    z_min = max(-α1, -α2, 0)
    z_max = min(β1, β2, β3)

    # Handle structural topological zeros (empty sum)
    if z_min > z_max
        empty_m = snapshot(SymbolicBuffer(1)) # Monomial with sign = 0
        return Generic6j(pref_sq, empty_m, CycloMonomial[], 0:-1)
    end

    # --- 3. Initial Summand M(z_min) ---
    cap_m = max(z_max + max(α1, α2), β1, β2, β3) + 1
    buf_m = SymbolicBuffer(cap_m)
    
    buf_m.sign = isodd(Int(z_min + α1 - α2)) ? -1 : 1
    add_qfact!(buf_m, z_min, -1)
    add_qfact!(buf_m, α1 + z_min, -1); add_qfact!(buf_m, α2 + z_min, -1)
    add_qfact!(buf_m, β1 - z_min, -1); add_qfact!(buf_m, β2 - z_min, -1); add_qfact!(buf_m, β3 - z_min, -1)
    
    m_min = snapshot(buf_m)

    # --- 4. Recursive Ratios R_z ---
    ratios = CycloMonomial[]
    sizehint!(ratios, z_max - z_min)
    
    # R_z is the ratio of term(z) / term(z-1)
    for z in z_min+1 : z_max
        buf_r = SymbolicBuffer(cap_m)
        buf_r.sign = -1 # Alternating sum gives a -1 factor between terms
        
        # Numerator of ratio (from the shifted denominator factorials of the 3j formula)
        add_qint!(buf_r, β1 - z + 1, 1)
        add_qint!(buf_r, β2 - z + 1, 1)
        add_qint!(buf_r, β3 - z + 1, 1)
        
        # Denominator of ratio (from the shifted numerator factorials)
        add_qint!(buf_r, z, -1)
        add_qint!(buf_r, α1 + z, -1)
        add_qint!(buf_r, α2 + z, -1)
        
        push!(ratios, snapshot(buf_r))
    end

    return Generic6j(pref_sq, m_min, ratios, z_min:z_max)
end

# ============================================================
# Classical Evaluation (q = 1) Limit
# ============================================================

# Fast prime power check
function is_prime_power(n::Int)
    # Checks if n is a prime power p^k and returns p.
    # Assumes Nemo is exported in the main module.
    facs = Nemo.factor(n)
    # Collect the factorization into an array of (prime, exponent) tuples
    pairs = collect(facs)
    
    # If there is exactly one prime base, it's a prime power
    if length(pairs) == 1
        return Int(pairs[1][1]) # Return the prime base
    end
    return 0
end

"""
    evaluate_classical(m::CycloMonomial)
Evaluates the symbolic monomial at the limit q -> 1.
"""
function evaluate_classical(m::CycloMonomial)
    m.sign == 0 && return 0.0
    
    # The value of \\Phi_d(1) is:
    # - p if d = p^k (a prime power)
    # - 1 if d has two or more distinct prime factors
    # - 0 if d = 1 (but our d starts at 2)
    
    log_val = 0.0
    for (d, exp) in enumerate(m.exps)
        (d < 2 || exp == 0) && continue
        
        p = is_prime_power(d) 
        if p > 0
            log_val += exp * log(Float64(p))
        end
    end
    
    return m.sign * exp(log_val)
end

function qracah6j_classical(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    q6j = qracah6j_generic(j1, j2, j3, j4, j5, j6)
    q6j.pref_sq.sign == 0 && return 0.0
    
    sumz = sum(evaluate_classical.(q6j.series))
    pref_sq_val = evaluate_classical(q6j.pref_sq)
    
    return Float64(sqrt(BigFloat(pref_sq_val)) * BigFloat(sumz))
end

function qracah3j_classical(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    q3j = qracah3j_generic(j1, j2, j3, m1, m2, m3)
    
    if q3j.pref_sq.sign == 0
        return 0.0
    end
    
    sumz = sum(evaluate_classical.(q3j.series))
    pref_sq_val = evaluate_classical(q3j.pref_sq)
    
    return Float64(sqrt(BigFloat(pref_sq_val)) * BigFloat(sumz))
end

#-------- Evaluations -------- 

# export evaluate_phi_stable

# """
#     evaluate_phi_stable(d::Int, k::Int)

# Computes the Log-Magnitude and Phase of Φ_d(q) at q = exp(i * 2π / (k+2)).
# Uses the trigonometric product formula for extreme numerical stability and caches the result.
# """
# # Efficiently evaluate Φ_d(x) using the product of (x^k - 1) for x root of unity
# # function evaluate_phi_stable(d::Int, k::Int)
# #     cache_key = (d, k)
# #     if haskey(PHI_EVAL_CACHE, cache_key)
# #         return PHI_EVAL_CACHE[cache_key]
# #     end

# #     # q = exp(iϕ) where ϕ = 2π / (k+2)
# #     ϕ = big(2) * big(π) / (k + 2)
    
# #     total_log_mag = big(0.0)
# #     total_phase = big(0.0)
    
# #     for n in 1:d
# #         if d % n == 0
# #             m = Nemo.moebius_mu(d ÷ n)
# #             m == 0 && continue

# #             # Φ_d(x) = Π_{n|d} (x^n - 1)^μ(d/n)
# #             # We evaluate |x^n - 1| and angle(x^n - 1)
# #             # |e^{inϕ} - 1| = |2 * sin(nϕ/2)|
# #             # angle(e^{inϕ} - 1) = (nϕ + π)/2
# #             val_sin = 2 * sin(n * ϕ / 2)
# #             mag_n = abs(val_sin)
# #             phase_n = (n * ϕ + big(π)) / 2
            
# #             # If sin() is negative, the abs() flipped the sign. 
# #             # We must add π to the phase to compensate.
# #             if val_sin < 0
# #                 phase_n += big(π)
# #             end
            
# #             total_log_mag += m * log(mag_n)
# #             total_phase += m * phase_n
# #         end
# #     end
    
# #     res = (total_log_mag, total_phase)
# #     PHI_EVAL_CACHE[cache_key] = res
# #     return res
# # end

# """
#     build_phi_table(D_max::Int, k::Int)

# Computes the log-magnitude and phase for ALL cyclotomic polynomials Φ_d(q) 
# from d = 1 to D_max simultaneously using the Fast Möbius Transform.
# """

# const GLOBAL_SIEVE_CACHE = Dict{Int, Tuple{Vector{BigFloat}, Vector{BigFloat}}}()

# """
#     get_phi_table(k::Int)
# Fetches or builds the Sieve table for level k. 
# Only takes 2.2ms the first time it is called for a new k!
# """
# function get_phi_table(k::Int)
#     haskey(GLOBAL_SIEVE_CACHE, k) && return GLOBAL_SIEVE_CACHE[k]
    
#     # The maximum d we will ever encounter for SU(2)_k is k+2
#     D_max = k + 2
#     lmag, lphs = build_phi_table(D_max, k)
    
#     GLOBAL_SIEVE_CACHE[k] = (lmag, lphs)
#     return (lmag, lphs)
# end

# export build_phi_table

# function build_phi_table(D_max::Int, k::Int)
#     h = k + 2
    
#     # 1. Allocate arrays for the results V(d)
#     V_mag = zeros(BigFloat, D_max)
#     V_phs = zeros(BigFloat, D_max)
    
#     # 2. Base Step: Compute F(n) exactly ONCE per n
#     # We do exactly D_max sine and log calls. Not a single one more.
#     for n in 1:D_max
#         # Handle the topological zero case safely
#         if n % h == 0
#             V_mag[n] = big(-Inf)
#             V_phs[n] = big(0.0)
#             continue
#         end

#         val_sin = 2 * sinpi(big(n) / h)
#         V_mag[n] = log(abs(val_sin))
        
#         p_n = (big(n) * big(π) / h + big(π)) / 2
#         if val_sin < 0
#             p_n += big(π)
#         end
#         V_phs[n] = p_n
#     end
    
#     # 3. The Fast Möbius Transform (Dirichlet Sieve)
#     # We solve F[n] = sum_{d|n} V[d] for V[d] in-place!
#     # No Nemo.divisors, no moebius_mu, just pure addition.
#     @inbounds for d in 1:D_max
#         # Subtract V[d] from all its multiples
#         for m in (2 * d):d:D_max
#             V_mag[m] -= V_mag[d]
#             V_phs[m] -= V_phs[d]
#         end
#     end
    
#     return V_mag, V_phs
# end



# function evaluate_phi_stable(d::Int, k::Int)
#     cache_key = (d, k)
#     if haskey(PHI_EVAL_CACHE, cache_key)
#         return PHI_EVAL_CACHE[cache_key]
#     end

#     h = k + 2
#     # Use sinpi for better accuracy and speed
#     # sin(n * pi / h) == sinpi(n / h)
    
#     total_log_mag = big(0.0)
#     total_phase = big(0.0)
    
#     # Φ_d(x) = Π_{n|d} (x^n - 1)^μ(d/n)
#     for n in Nemo.divisors(d)
#         m = Int(Nemo.moebius_mu(d ÷ n))
#         m == 0 && continue

#         # chord length = 2 * sin(n*π/h)
#         # We use BigFloat for the division to maintain precision
#         val_sin = 2 * sinpi(big(n) / h)
        
#         total_log_mag += m * log(abs(val_sin))
        
#         # phase = (n*π/h + π)/2
#         p_n = (big(n) * big(π) / h + big(π)) / 2
#         val_sin < 0 && (p_n += big(π))
        
#         total_phase += m * p_n
#     end
    
#     res = (total_log_mag, total_phase)
#     PHI_EVAL_CACHE[cache_key] = res
#     return res
# end


# """
#     mobius(n::Int)

# Computes the classical Möbius function μ(n) for integer factorization.
# Returns 0 if n has a squared prime factor, 1 if square-free with an even number of prime factors, 
# and -1 if square-free with an odd number of prime factors.
# """
# #compute the 
# function mobius(n::Int)
#     n == 1 && return 1
#     p = 0
#     for i in 2:isqrt(n)
#         if n % i == 0
#             (n % (i * i) == 0) && return 0
#             p += 1
#             n ÷= i
#             while n % i == 0
#                 n ÷= i
#             end
#         end
#     end
#     n > 1 && (p += 1)
#     return isodd(p) ? -1 : 1
# end

# """
#     evaluate_cyclotomic(d::Int, q::Number)

# Evaluates the d-th cyclotomic polynomial Φ_d(q) analytically at a complex or real parameter `q`.
# Utilizes Möbius inversion for strictly O(τ(d)) performance, bypassing polynomial expansion.
# """
# function evaluate_cyclotomic(d::Int, q::Number)
#     d == 1 && return q - one(q)
    
#     val = one(q)
#     for k in 1:d
#         if d % k == 0
#             mu = mobius(k)
#             if mu == 1
#                 val *= (q^(d ÷ k) - one(q))
#             elseif mu == -1
#                 val /= (q^(d ÷ k) - one(q))
#             end
#         end
#     end
#     return val
# end


# # ----- Internal Dispatch: Physical Level k (Discrete Geometry) --- #


# function evaluate_level(m::CycloMonomial, k::Int, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
#     m.sign == 0 && return zero(T)
    
#     # O(1) Algebraic Zero/Pole Check
#     if length(m.exps) >= k + 2
#         exponent = m.exps[k+2]
#         if exponent > 0
#             return zero(T) # Topological selection rule (exact zero)
#         elseif exponent < 0
#             throw(DomainError(k, "Topological pole: Level k=$k violates admissibility bounds."))
#         end
#     end
    
#     return setprecision(BigFloat, prec) do
#         θ = big(π) / (k + 2)
        
#         # Initialize with the overall sign and z^z_pow phase
#         log_mag = big(0.0)
#         phase = m.z_pow * θ
#         if m.sign == -1
#             phase += big(π)
#         end
        
#         # Add contributions from cyclotomic factors
#         for (d, e) in enumerate(m.exps)
#             (e == 0 || d < 1 || d == k + 2) && continue
            
#             # Fetch strictly by (d, k)
#             log_mag_phi, phase_phi = evaluate_phi_stable(d, k)
            
#             log_mag += e * log_mag_phi
#             phase += e * phase_phi
#         end
        
#         res = exp(log_mag) * cis(phase)
#         return T(res)
#     end
# end

# # function evaluate_level(res::GenericResult, k::Int, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
# #     res.pref_sq.sign == 0 && return zero(T)

# #     return setprecision(BigFloat, prec) do
# #         # 1. Squared prefactor -> take sqrt(abs()) since physical volume^2 is real
# #         val_pref_sq = evaluate_level(res.pref_sq, k, Complex{BigFloat}; prec=prec)
        
# #         # If the prefactor evaluates to an exact zero (via the d=k+2 check), abort early
# #         iszero(val_pref_sq) && return zero(T)
        
# #         val_pref = sqrt(abs(val_pref_sq))

# #         # 2. Sum the Racah series (Call internal function directly for performance)
# #         val_sum = zero(Complex{BigFloat})
# #         for term in res.series
# #             val_sum += evaluate_level(term, k, Complex{BigFloat}; prec=prec)
# #         end

# #         # 3. Combine
# #         final_val = val_pref * val_sum
        
# #         # 4. Safe casting
# #         if T <: Real
# #             return T(real(final_val))
# #         else
# #             return T(final_val)
# #         end
# #     end
# # end


# # function evaluate_level(res::GenericResult, k::Int, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
# #     res.pref_sq.sign == 0 && return zero(T)

# #     return setprecision(BigFloat, prec) do
# #         h = k + 2
# #         θ = big(π) / h

# #         # --- STEP 1: Batch Specialization (Pre-compute all Φ values) ---
# #         # Find max d across the whole result
# #         max_d = length(res.pref_sq.exps)
# #         for m in res.series
# #             max_d = max(max_d, length(m.exps))
# #         end

# #         # Pre-calculate log-mag and phase for every d needed
# #         # We use Vectors for cache-friendly O(1) access
# #         phi_lmag = zeros(BigFloat, max_d)
# #         phi_lphs = zeros(BigFloat, max_d)
        
# #         # Identify indices to avoid redundant trig calls
# #         needed = BitSet()
# #         _scan!(m) = for (d, e) in enumerate(m.exps); (e != 0 && d != h) && push!(needed, d); end
# #         _scan!(res.pref_sq)
# #         for m in res.series; _scan!(m); end

# #         for d in needed
# #             m_v, p_v = evaluate_phi_stable(d, k) # Uses warm cache internally
# #             phi_lmag[d] = m_v
# #             phi_lphs[d] = p_v
# #         end

# #         # --- STEP 2: Inner Evaluation Helper ---
# #         # We calculate (log_mag, phase) for a monomial to stay in real-arithmetic as long as possible
# #         function get_log_coords(m::CycloMonomial)
# #             # Check topological zero for SU(2)_k
# #             (length(m.exps) >= h && m.exps[h] > 0) && return (big(-Inf), big(0.0))
            
# #             lm = big(0.0)
# #             lp = m.z_pow * θ
# #             m.sign == -1 && (lp += big(π))
            
# #             for (d, e) in enumerate(m.exps)
# #                 (e == 0 || d == h) && continue
# #                 lm += e * phi_lmag[d]
# #                 lp += e * phi_lphs[d]
# #             end
# #             return (lm, lp)
# #         end

# #         # --- STEP 3: Summation ---
# #         val_sum = zero(Complex{BigFloat})
# #         for m in res.series
# #             lm, lp = get_log_coords(m)
# #             lm == -Inf && continue # Skip topological zeros
# #             val_sum += exp(lm) * cis(lp)
# #         end

# #         # --- STEP 4: Prefactor & Finalize ---
# #         pref_lm, pref_lp = get_log_coords(res.pref_sq)
# #         # Prefactor is squared in the struct, so we take sqrt(exp(log))
# #         val_pref = sqrt(exp(pref_lm)) 
        
# #         final = val_pref * val_sum
# #         return T <: Real ? T(real(final)) : T(final)
# #     end
# # end

# # struct SparseTerm
# #     sign::Int8
# #     z_pow::Int
# #     active::Vector{Pair{Int, Int}}
# # end

# # function evaluate_level(res::GenericResult, k::Int, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
# #     res.pref_sq.sign == 0 && return zero(T)

# #     return setprecision(BigFloat, prec) do
# #         h = k + 2
# #         θ = big(π) / h

# #         # --- PHASE 1: Build Local Dense Table (Fastest Access) ---
# #         max_d = length(res.pref_sq.exps)
# #         for m in res.series; max_d = max(max_d, length(m.exps)); end

# #         # Pre-compute Φ factors into Vectors
# #         lmag_table = zeros(BigFloat, max_d)
# #         lphs_table = zeros(BigFloat, max_d)
        
# #         needed = BitSet()
# #         _mark!(m) = for (d, e) in enumerate(m.exps); (e != 0 && d != h) && push!(needed, d); end
# #         _mark!(res.pref_sq)
# #         for m in res.series; _mark!(m); end

# #         for d in needed
# #             m_v, p_v = evaluate_phi_stable(d, k) 
# #             lmag_table[d] = m_v
# #             lphs_table[d] = p_v
# #         end

# #         # --- PHASE 2: Sparse Representation ---
# #         # Convert dense monomials into sparse index-exponent pairs 
# #         # to avoid iterating over thousands of zeros in the hot loop.
        

# #         function to_sparse(m::CycloMonomial)
# #             act = [d => m.exps[d] for d in 1:length(m.exps) if m.exps[d] != 0 && d != h]
# #             return SparseTerm(Int8(m.sign), m.z_pow, act)
# #         end

# #         sparse_series = [to_sparse(m) for m in res.series]
# #         sparse_pref = to_sparse(res.pref_sq)

# #         # --- PHASE 3: The Hot Summation ---
# #         # This loop is designed to be as "flat" as possible.
# #         val_sum = zero(Complex{BigFloat})
        
# #         # Temp variables for log-space calculation
# #         for st in sparse_series
# #             # Use local real variables (Fast)
# #             lm = big(0.0)
# #             lp = st.z_pow * θ
# #             st.sign == -1 && (lp += big(π))
            
# #             for (d, e) in st.active
# #                 lm += e * lmag_table[d]
# #                 lp += e * lphs_table[d]
# #             end
            
# #             # This is the ONLY place where we enter the heap-allocated value space
# #             val_sum += exp(lm) * cis(lp)
# #         end

# #         # --- PHASE 4: Finalize ---
# #         p_lm = big(0.0)
# #         p_lp = sparse_pref.z_pow * θ
# #         sparse_pref.sign == -1 && (p_lp += big(π))
# #         for (d, e) in sparse_pref.active
# #             p_lm += e * lmag_table[d]
# #         end
        
# #         val_pref = sqrt(exp(p_lm))
# #         return T(val_pref * val_sum)
# #     end
# # end





# # function evaluate_level(res::GenericResult, k::Int, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
# #     res.pref_sq.sign == 0 && return zero(T)

# #     return setprecision(BigFloat, prec) do
# #         h = k + 2
# #         θ = big(π) / h

# #         # --- PHASE 1: Fast Table Generation via Dirichlet Sieve ---
# #         # Find the maximum cyclotomic index needed
# #         max_d = length(res.pref_sq.exps)
# #         for m in res.series
# #             max_d = max(max_d, length(m.exps))
# #         end

# #         # Generate the entire table in one shot (The 2.2ms step!)
# #         lmag_table, lphs_table = build_phi_table(max_d, k)

# #         # --- PHASE 2: Sparse Representation ---
# #         # Convert dense monomials into sparse (index => exponent) arrays.
# #         # This prevents the inner loop from checking thousands of zeros.
# #         function to_sparse(m::CycloMonomial)
# #             act = Pair{Int, Int}[]
# #             for (d, e) in enumerate(m.exps)
# #                 if e != 0 && d != h
# #                     push!(act, d => e)
# #                 end
# #             end
# #             # Using a NamedTuple for clean, allocation-free grouping
# #             return (sign=Int8(m.sign), z_pow=m.z_pow, active=act)
# #         end

# #         sparse_series = [to_sparse(m) for m in res.series]
# #         sparse_pref = to_sparse(res.pref_sq)

# #         # --- PHASE 3: The Hot Summation ---
# #         val_sum = zero(Complex{BigFloat})
        
# #         for st in sparse_series
# #             # Local real variables for fast log-space accumulation
# #             lm = big(0.0)
# #             lp = st.z_pow * θ
# #             st.sign == -1 && (lp += big(π))
            
# #             # The innermost loop: Pure Real Addition
# #             for (d, e) in st.active
# #                 lm += e * lmag_table[d]
# #                 lp += e * lphs_table[d]
# #             end
            
# #             # We only convert back to Value-Space ONCE per term
# #             val_sum += exp(lm) * cis(lp)
# #         end

# #         # --- PHASE 4: Finalize ---
# #         p_lm = big(0.0)
# #         # Prefactor phase vanishes under the physical abs(), so we only need magnitude
# #         for (d, e) in sparse_pref.active
# #             p_lm += e * lmag_table[d]
# #         end
        
# #         # val_pref = sqrt(abs(val_pref_sq)) = sqrt(exp(p_lm)) = exp(p_lm / 2)
# #         val_pref = exp(p_lm / 2)
        
# #         final = val_pref * val_sum
# #         return T <: Real ? T(real(final)) : T(final)
# #     end
# # end



# function evaluate_level(res::GenericResult, k::Int, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
#     res.pref_sq.sign == 0 && return zero(T)

#     return setprecision(BigFloat, prec) do
#         h = k + 2
#         θ = big(π) / h

#         # --- PHASE 1: O(1) Table Lookup ---
#         # 0.0 ms if cached, 2.2 ms if cold!
#         lmag_table, lphs_table = get_phi_table(k)

#         # --- PHASE 2: Differential Summation (The Magic) ---
#         val_sum = zero(Complex{BigFloat})
#         isempty(res.series) && return zero(T)

#         # Evaluate the FIRST term completely
#         prev_m = res.series[1]
#         lm = big(0.0)
#         lp = prev_m.z_pow * θ
#         prev_m.sign == -1 && (lp += big(π))
        
#         for (d, e) in enumerate(prev_m.exps)
#             (e == 0 || d == h) && continue
#             lm += e * lmag_table[d]
#             lp += e * lphs_table[d]
#         end
        
#         val_sum += exp(lm) * cis(lp)

#         # Evaluate ALL SUBSEQUENT terms differentially
#         # We only add the differences in exponents, skipping 95% of the math.
#         @inbounds for i in 2:length(res.series)
#             curr_m = res.series[i]
            
#             # 1. Update Phase by Sign Difference
#             if curr_m.sign != prev_m.sign
#                 lp += big(π) # Adding pi flips the sign of cis()
#             end
            
#             # 2. Update Phase by z_pow Difference
#             dz = curr_m.z_pow - prev_m.z_pow
#             if dz != 0
#                 lp += dz * θ
#             end
            
#             # 3. Update by Cyclotomic Factor Differences
#             len_prev = length(prev_m.exps)
#             len_curr = length(curr_m.exps)
#             max_len = max(len_prev, len_curr)
            
#             for d in 1:max_len
#                 d == h && continue
#                 e_prev = d <= len_prev ? prev_m.exps[d] : 0
#                 e_curr = d <= len_curr ? curr_m.exps[d] : 0
#                 de = e_curr - e_prev
                
#                 # 'de' is only non-zero for ~8 elements out of hundreds!
#                 if de != 0
#                     lm += de * lmag_table[d]
#                     lp += de * lphs_table[d]
#                 end
#             end
            
#             # The only heavy BigFloat allocation left
#             val_sum += exp(lm) * cis(lp)
#             prev_m = curr_m
#         end

#         # --- PHASE 3: Finalize ---
#         p_lm = big(0.0)
#         for (d, e) in enumerate(res.pref_sq.exps)
#             (e == 0 || d == h) && continue
#             p_lm += e * lmag_table[d]
#         end
        
#         val_pref = exp(p_lm / 2)
#         final = val_pref * val_sum
        
#         return T <: Real ? T(real(final)) : T(final)
#     end
# end



















# # ----- Internal Dispatch: Analytic Continuation (Continuous Parameter) --- #

# # function evaluate_analytic(m::CycloMonomial, q::Number, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
# #     m.sign == 0 && return zero(T)
    
# #     return setprecision(BigFloat, prec) do
# #         # Evaluate z = q^(1/2) for the phase prefactor
# #         z = sqrt(Complex(q))
# #         val = Complex{BigFloat}(m.sign) * (z ^ m.z_pow)
        
# #         for (d, e) in enumerate(m.exps)
# #             (e == 0 || d < 1) && continue
# #             poly_val = evaluate_cyclotomic(d, q) # Evaluates Φ_d(q)
# #             val *= poly_val ^ e
# #         end
        
# #         return T(val)
# #     end
# # end

# # function evaluate_analytic(res::GenericResult, q::Number, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
# #     res.pref_sq.sign == 0 && return zero(T)
    
# #     return setprecision(BigFloat, prec) do
# #         # 1. Squared prefactor
# #         val_pref_sq = evaluate_analytic(res.pref_sq, q, Complex{BigFloat}; prec=prec)

# #         # Abort early if the polynomial evaluates exactly to zero at this arbitrary q
# #         iszero(val_pref_sq) && return zero(T)

# #         # Do NOT take abs() here; preserve complex branch cuts for analytic continuation
# #         val_pref = sqrt(val_pref_sq)

# #         # 2. Sum the Racah series (Call internal function directly for performance)
# #         val_sum = zero(Complex{BigFloat})
# #         for term in res.series
# #             val_sum += evaluate_analytic(term, q, Complex{BigFloat}; prec=prec)
# #         end

# #         # 3. Combine
# #         final_val = val_pref * val_sum
        
# #         return T(final_val)
# #     end
# # end

# function evaluate_analytic(m::CycloMonomial, q::Number, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
#     m.sign == 0 && return zero(T)
    
#     return setprecision(BigFloat, prec) do
#         # Force q into BigFloat precision immediately
#         q_big = Complex{BigFloat}(q)
        
#         # Initialize Logarithmic Accumulators
#         log_mag = big(0.0)
#         phase = big(0.0)
        
#         # 1. Sign Phase
#         if m.sign == -1
#             phase += big(π)
#         end
        
#         # 2. Prefactor z^z_pow (where z = sqrt(q))
#         z = sqrt(q_big)
#         log_z = log(z)
#         log_mag += m.z_pow * real(log_z)
#         phase += m.z_pow * imag(log_z)
        
#         # 3. Cyclotomic Factors (Log-Sum-Exp over the complex plane)
#         for (d, e) in enumerate(m.exps)
#             (e == 0 || d < 1) && continue
            
#             poly_val = evaluate_cyclotomic(d, q_big)
#             log_poly = log(poly_val)
            
#             log_mag += e * real(log_poly)
#             phase += e * imag(log_poly)
#         end
        
#         # Recombine safely
#         return T(exp(log_mag) * cis(phase))
#     end
# end

# function evaluate_analytic(res::GenericResult, q::Number, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
#     res.pref_sq.sign == 0 && return zero(T)
    
#     return setprecision(BigFloat, prec) do
#         # 1. Squared prefactor
#         val_pref_sq = evaluate_analytic(res.pref_sq, q, Complex{BigFloat}; prec=prec)

#         # Abort early if exact zero
#         iszero(val_pref_sq) && return zero(T)

#         # Do NOT take abs() here; preserve complex branch cuts for analytic continuation!
#         val_pref = sqrt(val_pref_sq)

#         # 2. Sum the Racah series (Call internal function directly for performance)
#         val_sum = zero(Complex{BigFloat})
#         for term in res.series
#             val_sum += evaluate_analytic(term, q, Complex{BigFloat}; prec=prec)
#         end

#         # 3. Combine
#         final_val = val_pref * val_sum
        
#         return T(final_val)
#     end
# end

# # Public API: Unified Evaluation Interface

# """
#     evaluate_generic(m::CycloMonomial, ::Type{T}=Complex{BigFloat}; k=nothing, q=nothing, prec=512)

# Evaluates a symbolic CycloMonomial. You must specify exactly one evaluation target:
# - `k=val`: Evaluates at the physical root of unity q = exp(iπ/(k+2)), enforcing SU(2)_k topological checks.
# - `q=val`: Evaluates analytically at a generic complex or real parameter.
# """
# function evaluate_generic(m::CycloMonomial, ::Type{T}=Complex{BigFloat}; k=nothing, q=nothing, prec=512) where {T}
#     if (isnothing(k) && isnothing(q)) || (!isnothing(k) && !isnothing(q))
#         throw(ArgumentError("Ambiguous evaluation. You must specify exactly one target: either `k=val` or `q=val`."))
#     end

#     if !isnothing(k)
#         return evaluate_level(m, k, T; prec=prec)
#     else
#         return evaluate_analytic(m, q, T; prec=prec)
#     end
# end

# """
#     evaluate_generic(res::GenericResult, ::Type{T}=Complex{BigFloat}; k=nothing, q=nothing, prec=512)

# Evaluates a full 3j or 6j GenericResult. You must specify exactly one evaluation target:
# - `k=val`: Evaluates at the physical root of unity, ensuring exact algebraic cancellation.
# - `q=val`: Analytically continues the quantum symbol to an arbitrary complex parameter.
# """
# function evaluate_generic(res::GenericResult, ::Type{T}=Complex{BigFloat}; k=nothing, q=nothing, prec=512) where {T}
#     if (isnothing(k) && isnothing(q)) || (!isnothing(k) && !isnothing(q))
#         throw(ArgumentError("Ambiguous evaluation. You must specify exactly one target: either `k=val` or `q=val`."))
#     end

#     if !isnothing(k)
#         return evaluate_level(res, k, T; prec=prec)
#     else
#         return evaluate_analytic(res, q, T; prec=prec)
#     end
# end