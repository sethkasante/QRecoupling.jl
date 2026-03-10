
#Symbolics.jl

# ============================================================
# Core Allocation-Free In-Place Operations
# ============================================================

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

# ============================================================
# Symbolic Engine (Single-Allocation Setup)
# ============================================================

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

# ============================================================
# Racah Series (Hypergeometric Ratio Method)
# ============================================================

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

# ============================================================
# Quantum 3j symbols (Optimized)
# ============================================================

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

"""
    evaluate_phi_stable(d::Int, k::Int)

Computes the Log-Magnitude and Phase of Φ_d(q) at q = exp(i * 2π / (k+2)).
Uses the trigonometric product formula for extreme numerical stability and caches the result.
"""
# Efficiently evaluate Φ_d(x) using the product of (x^k - 1) for x root of unity
function evaluate_phi_stable(d::Int, k::Int)
    cache_key = (d, k)
    if haskey(PHI_EVAL_CACHE, cache_key)
        return PHI_EVAL_CACHE[cache_key]
    end

    # q = exp(iϕ) where ϕ = 2π / (k+2)
    ϕ = big(2) * big(π) / (k + 2)
    
    total_log_mag = big(0.0)
    total_phase = big(0.0)
    
    for n in 1:d
        if d % n == 0
            m = Nemo.moebius_mu(d ÷ n)
            m == 0 && continue

            # Φ_d(x) = Π_{n|d} (x^n - 1)^μ(d/n)
            # We evaluate |x^n - 1| and angle(x^n - 1)
            # |e^{inϕ} - 1| = |2 * sin(nϕ/2)|
            # angle(e^{inϕ} - 1) = (nϕ + π)/2
            val_sin = 2 * sin(n * ϕ / 2)
            mag_n = abs(val_sin)
            phase_n = (n * ϕ + big(π)) / 2
            
            # If sin() is negative, the abs() flipped the sign. 
            # We must add π to the phase to compensate.
            if val_sin < 0
                phase_n += big(π)
            end
            
            total_log_mag += m * log(mag_n)
            total_phase += m * phase_n
        end
    end
    
    res = (total_log_mag, total_phase)
    PHI_EVAL_CACHE[cache_key] = res
    return res
end


function evaluate_level(m::CycloMonomial, k::Int, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
    m.sign == 0 && return zero(T)
    # O(1) Algebraic Zero/Pole Check
    if length(m.exps) >= k + 2
        exponent = m.exps[k+2]
        if exponent > 0
            return zero(T) # Topological selection rule (exact zero)
        elseif exponent < 0
            throw(DomainError(k, "Topological pole: Level k=$k violates admissibility bounds."))
        end
    end
    
    return setprecision(BigFloat, prec) do
        θ = big(π) / (k + 2)
        
        # Initialize with the overall sign and z^z_pow phase
        log_mag = big(0.0)
        phase = m.z_pow * θ
        if m.sign == -1
            phase += big(π)
        end
        
        # Add contributions from cyclotomic factors
        for (d, e) in enumerate(m.exps)
            (e == 0 || d < 1) && continue
            
            # Fetch strictly by (d, k)
            log_mag_phi, phase_phi = evaluate_phi_stable(d, k)
            
            log_mag += e * log_mag_phi
            phase += e * phase_phi
        end
        
        res = exp(log_mag) * cis(phase)
        return T(res)
    end
end

function evaluate_analytic(m::CycloMonomial, q::Number, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
    m.sign == 0 && return zero(T)
    
    # Evaluate z = q^(1/2) for the phase prefactor
    z = sqrt(Complex(q))
    val = Complex{BigFloat}(m.sign) * (z ^ m.z_pow)
    
    for (d, e) in enumerate(m.exps)
        (e == 0 || d < 1) && continue
        poly_val = evaluate_cyclotomic(d, q) # Evaluates Φ_d(q)
        val *= poly_val ^ e
    end
    
    return T(val)
end

"""
    evaluate_generic(m::CycloMonomial, k::Int, ::Type{T}=Complex{BigFloat}; prec=512)

Evaluates a symbolic cyclotomic monomial using log-sum-exp for numerical stability.
"""
function evaluate_generic(m::CycloMonomial, ::Type{T}=Complex{BigFloat}; k=OptInt, q=OptInt, prec=512) where {T}
    # 1. Intent Validation
    if (isnothing(k) && isnothing(q)) || (!isnothing(k) && !isnothing(q))
        throw(ArgumentError("Ambiguous evaluation. You must specify exactly one target: either `k=val` or `q=val`."))
    end

    # 2. Physical Level Evaluation (Discrete Geometry)
    if !isnothing(k)
        evaluate_level(m, k, T; prec=prec)
    else
        evaluate_analytic(m, q, T; prec=prec)
    end
end


function evaluate_level(res::GenericResult, k::Int, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
    res.pref_sq.sign == 0 && return zero(T)

    return setprecision(BigFloat, prec) do
        # 1. Squared prefactor -> take sqrt(abs())
        val_pref_sq = evaluate_level(res.pref_sq, k, Complex{BigFloat}; prec=prec)
        val_pref = sqrt(abs(val_pref_sq))

        # 2. Sum the Racah series
        val_sum = zero(Complex{BigFloat})
        for term in res.series
            val_sum += evaluate_generic(term, k, Complex{BigFloat}; prec=prec)
        end

        # 3. Combine
        final_val = val_pref * val_sum
        
        # 4. Safe casting
        if T <: Real
            return T(real(final_val))
        else
            return T(final_val)
        end
    end
end


function evaluate_analytic(res::GenericResult, q::Number, ::Type{T}=Complex{BigFloat}; prec=512) where {T}
    res.pref_sq.sign == 0 && return zero(T)
   
    
    return setprecision(BigFloat, prec) do
        # 1. Squared prefactor -> take sqrt(abs())
        val_pref_sq = evaluate_analytic(res.pref_sq, q, Complex{BigFloat}; prec=prec)

        # If the prefactor evaluates to an exact zero (via the d=k+2 check), abort early
        iszero(val_pref_sq) && return zero(T)

        val_pref = sqrt(abs(val_pref_sq))

        # 2. Sum the Racah series
        val_sum = zero(Complex{BigFloat})
        for term in res.series
            val_sum += evaluate_generic(term, q, Complex{BigFloat}; prec=prec)
        end

        # 3. Combine
        final_val = val_pref * val_sum
        
        # 4. Safe casting
        if T <: Real
            return T(real(final_val))
        else
            return T(final_val)
        end
    end
end


"""
    evaluate_generic(res::GenericResult, k::Int, ::Type{T}=Complex{BigFloat}; prec=512)

Evaluates a full 3j or 6j GenericResult at level k.
"""
function evaluate_generic(res::GenericResult, ::Type{T}=Complex{BigFloat}; k=OptInt, q=OptInt, prec=512) where {T}
    # 1. Intent Validation
    if (isnothing(k) && isnothing(q)) || (!isnothing(k) && !isnothing(q))
        throw(ArgumentError("Ambiguous evaluation. You must specify exactly one target: either `k=val` or `q=val`."))
    end

    # 2. Physical Level Evaluation (Discrete Geometry)
    if !isnothing(k)
        evaluate_level(res, k, T; prec=prec)
    else
        evaluate_analytic(res, q, T; prec=prec)
    end
end
