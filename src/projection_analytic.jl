
# ---------------------------------------------------------------------------
#                   --- Analytic Continuation ----
# Project Deferred Cyclotomic Representations (DCR) analytically for q ∈ ℂ
# Convention: [n]_q =q^{n}-q^{-n} / (q-q^{-1}) = q^{1-n} * \prod \Phi_d(q^2)
# Optimized for complex phase tracking and branch-cut stability.
# Thread-Safe
# ---------------------------------------------------------------------------




# ---- On-The-Fly: cyclotomic table builder --- --

"""
    build_analytic_table(max_d::Int, q_sq::T) where T <: Number
Computes exact complex values of Φ_d(q²) up to max_d.
"""
function build_analytic_table(max_d::Int, q_sq::T) where T
    table = Vector{T}(undef, max_d)
    max_d == 0 && return table
    
    # Φ_1(q^2) = q^2 - 1
    table[1] = q_sq - one(T)
    
    q_curr = q_sq * q_sq
    @inbounds for n in 2:max_d
        val = q_curr - one(T) # (q^2)^n - 1
        
        # Factor out the lower divisors
        for d in 1:(n-1)
            if n % d == 0
                if iszero(table[d])
                     val = zero(T)
                     break
                end
                val /= table[d]
            end
        end
        table[n] = val
        q_curr *= q_sq
    end
    return table
end




# -- Internal computations ---


"""
    _eval_mono_analytic(m::CyclotomicMonomial, q::T, table::Vector{T}) where T
Evaluates the monomial directly using exact multiplication.
M =  q^{q_pow} * ∏ Φ_d(q²)^e
"""
@inline function _eval_mono_analytic(m::CyclotomicMonomial, q::T, table::Vector{T}) where T
    m.sign == 0 && return zero(T)
    
    val = one(T)
    pe = m.phi_exps
    
    @inbounds for i in 1:length(pe)
        p = pe[i]
        d, e = p.first, p.second
        
        # unroll exponent to bypass complex `^` operator
        if e == 1
            val *= table[d]
        elseif e == -1
            iszero(table[d]) && throw(DomainError(d, "Division by zero: Φ_$d(q²) = 0."))
            val /= table[d]
        elseif e == 2
            val *= (table[d] * table[d])
        elseif e == -2
            iszero(table[d]) && throw(DomainError(d, "Division by zero: Φ_$d(q²) = 0."))
            val /= (table[d] * table[d])
        elseif e > 0
            val *= table[d]^e
        else
            iszero(table[d]) && throw(DomainError(d, "Division by zero: Φ_$d(q²) = 0."))
            val /= table[d]^(-e)
        end
    end
    
    return m.sign * (q^m.q_pow) * val
end


# -----  Core Evaluator (Single-Pass) --- 

function _eval_analytic_dcr(res::DCR, q::T, q_sq::T) where T
    # for cyclotomic polynomials
    table = build_analytic_table(res.max_d, q_sq)

    # project prefactors
    val_root = _eval_mono_analytic(res.root, q, table)
    val_rad  = _eval_mono_analytic(res.radical, q, table)
    val_base = _eval_mono_analytic(res.base, q, table)
    
    # root * √(radical) * base
    pref_val = val_root * sqrt(val_rad) * val_base

    iszero(pref_val) && return zero(T)

    sum_val = pref_val
    curr_term = pref_val
    
    ratios = res.ratios
    @inbounds for i in 1:length(ratios)
        r_val = _eval_mono_analytic(ratios[i], q, table)
        curr_term *= r_val
        iszero(curr_term) && break 
        sum_val += curr_term
    end

    return sum_val
end


# ---  Projection of cyclotomic monomials  -----

"""
    project_analytic(m::CyclotomicMonomial, q::Number)
Evaluates a single quantum monomial
"""
function project_analytic(m::CyclotomicMonomial, q::Number)
    T = typeof(q * 1.0) 
    
    m.sign == 0 && return zero(T)
    abs(q - 1.0) < 1e-12 && throw(ArgumentError("For q=1, use mode=:classical instead."))
    
    q_T = T(q)
    q_sq = q_T * q_T
    
    table = build_analytic_table(m.max_d, q_sq)
    return _eval_mono_analytic(m, q_T, table)
end

# ---  Projection of DCR -----

"""
    project_analytic(dcr::DCR, q::Number)
Fast, thread-safe for evaluating TQFT symbols at any generic q.
"""
function project_analytic(dcr::DCR, q::Number)
    T = typeof(q * 1.0)
    
    dcr.base.sign == 0 && return zero(T)
    abs(q - 1.0) < 1e-12 && throw(ArgumentError("For q=1, use mode=:classical instead."))
    
    q_T = T(q)
    q_sq = q_T * q_T

    return _eval_analytic_dcr(dcr, q_T, q_sq)
end