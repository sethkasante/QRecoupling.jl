
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

#--------- Evaluations -------------- 

# Efficiently evaluate Φ_d(x) using the product of (x^k - 1) formula
function evaluate_phi(d::Int, x::Complex{BigFloat})
    # Φ_d(x) = Π_{n|d} (x^n - 1)^μ(d/n)
    # This is much faster than expanding the polynomial coefficients
    res = Complex{BigFloat}(1.0)
    for n in 1:d
        if d % n == 0
            term = (x^n - 1.0)
            m = moebius_mu(d ÷ n) # Use Nemo.moebius_mu
            if m == 1
                res *= term
            elseif m == -1
                res /= term
            end
        end
    end
    return res
end

"""
    evaluate_symbolic(m::CycloMonomial, k::Int; T=Complex{BigFloat})
Evaluates a symbolic cyclotomic monomial at the root of unity q = exp(iπ/(k+2)).
"""
function evaluate_symbolic(m::CycloMonomial, k::Int, ::Type{T}=Complex{BigFloat}) where {T}
    m.sign == 0 && return zero(T)
    
    # Root of unity for level k
    θ = BIG_PI / (k + 2)
    z = cis(θ) # z = q^(1/2)
    
    # 1. Start with the sign and z-power
    # Use cispi or exp for high precision
    res = Complex{BigFloat}(m.sign * z^m.z_pow)
    
    # 2. Multiply by each Φ_d(q)
    for (d, exp_val) in enumerate(m.exps)
        exp_val == 0 && continue
        
        # Evaluate the d-th cyclotomic polynomial at q = z^2
        # Φ_d(q) = Π (q - ζ_d^k) where gcd(k,d)=1
        val_phi = evaluate_phi(d, z^2)
        res *= val_phi^exp_val
    end
    
    return T(res)
end

"""
    evaluate_symbolic(res::GenericResult, k::Int, ::Type{T}=Complex{BigFloat})
Evaluates a symbolic 6j or 3j result at the level-k root of unity q = exp(iπ/(k+2)).
Returns a numerical value of type T.
"""
function evaluate_symbolic(res::GenericResult, k::Int, ::Type{T}=Complex{BigFloat}) where {T}
    # 1. Handle the case where the symbol is identically zero (admissibility)
    if res.pref_sq.sign == 0
        return zero(T)
    end

    # 2. Evaluate the squared prefactor
    # Mathematically, for SU(2)_k, the prefactor squared is always real and non-negative
    val_pref_sq = evaluate_symbolic(res.pref_sq, k, Complex{BigFloat})
    val_pref = sqrt(abs(val_pref_sq))

    # 3. Evaluate and sum the Racah series terms
    val_sum = zero(Complex{BigFloat})
    for term in res.series
        val_sum += evaluate_symbolic(term, k, Complex{BigFloat})
    end

    # 4. Final combination and cast to requested type T
    # Most 6j/3j are real, but we return Complex by default to handle R-matrices
    final_val = val_pref * val_sum
    
    if T <: Real
        return T(real(final_val))
    else
        return T(final_val)
    end
end
