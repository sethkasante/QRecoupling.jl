# -------------------------------------------------------------------------------
# Discrete Evaluator: SU(2)_k Roots of Unity
# Maps symbolic CycloResults into high-precision floats at exact roots of unity.
# -------------------------------------------------------------------------------

# Cache the magnitude table based on type T and level k. 
const GLOBAL_SIEVE_CACHE = LRU{Tuple{DataType, Int}, Vector}(maxsize=4096)

"""
    build_phi_table(D_max::Int, k::Int)

Computes the log-magnitude for cyclotomic polynomials Φ_d(q) using a 
Fast Möbius Transform natively in the precision of `T`.
"""
function build_phi_table(D_max::Int, k::Int, ::Type{T}) where {T}
    h = k + 2
    V_mag = zeros(T, D_max)
    V_valid = trues(D_max) 
    
    h_T = T(h)
    
    for n in 1:D_max
        if n % h == 0
            V_valid[n] = false
            continue
        end

        # Natively compute in T to ensure fair precision benchmarking
        val_sin = 2 * sinpi(n / h_T)
        V_mag[n] = log(abs(val_sin))
    end
    
    # Fast Möbius Transform for Cyclotomic factor log-magnitudes
    @inbounds for d in 1:D_max
        !V_valid[d] && continue
        for m in (2 * d):d:D_max
            if V_valid[m]
                V_mag[m] -= V_mag[d]
            end
        end
    end
    
    return V_mag
end

function get_phi_table(k::Int, D_max::Int, ::Type{T}) where {T}
    key = (T, k)
    if haskey(GLOBAL_SIEVE_CACHE, key)
        lmag = GLOBAL_SIEVE_CACHE[key]::Vector{T}
        if length(lmag) >= D_max
            return lmag
        end
    end
    
    lmag = build_phi_table(max(D_max, k + 2), k, T)
    GLOBAL_SIEVE_CACHE[key] = lmag
    return lmag
end

# -------------------------------------------
# High-Performance Discrete Projection
# -------------------------------------------

"""
    _project_to_real(m::CycloMonomial, lmag_table::Vector{T}, h::Int, ::Type{T})

The ultra-fast hot loop. Evaluates the cyclotomic log-magnitudes strictly in `T`.
"""
@inline function _project_to_real(m::CycloMonomial, lmag_table::Vector{T}, h::Int, ::Type{T}) where {T}
    m.sign == 0 && return zero(T)

    exps = m.exps

    @inbounds for i in eachindex(exps)
        d, e = exps[i]
        if d >= h
            e > 0 && return zero(T)
            e < 0 && throw(DomainError(h-2, "Topological pole: Level k=$(h-2)."))
        end
    end
    
    lm = zero(T)
    @inbounds @simd for i in eachindex(exps)
        d, e = exps[i]
        lm += e * lmag_table[d]
    end
    
    val = exp(lm)
    return m.sign == 1 ? val : -val
end

@inline function _project_to_log(m::CycloMonomial, lmag_table::Vector{T}, h::Int, ::Type{T}) where {T}
    # If the base monomial is zero, its log is -Infinity
    m.sign == 0 && return typemin(T) 
    
    lm = zero(T)
    exps = m.exps
    
    # Check for topological zeros/poles
    @inbounds for i in eachindex(exps)
        d, e = exps[i]
        if d >= h
            e > 0 && return typemin(T) # zero -> log_mag = -Inf
            e < 0 && throw(DomainError(h-2, "Topological pole: Level k=$(h-2)."))
        end
    end
    
    @inbounds @simd ivdep for i in eachindex(exps)
        d, e = exps[i]
        lm += e * lmag_table[d]
    end
    
    return lm
end

#core evaluation of CycloResult
function _eval_cresult(res::CycloResult, k::Int, ::Type{T}, buffer::Vector{T}) where {T}
    h = k + 2
    lmag_table = get_phi_table(k, res.max_d, T)

    base_lm = _project_to_log(res.base_term, lmag_table, h, T)
    base_lm == typemin(T) && return zero(T)
    
    root_lm = _project_to_log(res.root, lmag_table, h, T)
    root_lm == typemin(T) && return zero(T)
    
    rad_lm = _project_to_log(res.radical, lmag_table, h, T)
    rad_lm == typemin(T) && return zero(T) 
    
    pref_lm = root_lm + (rad_lm / 2) + base_lm
    pref_sign = res.root.sign * res.base_term.sign 

    num_ratios = length(res.ratios)
    
    max_lm = pref_lm
    curr_lm = pref_lm
    actual_len = num_ratios
    
    @inbounds for i in 1:num_ratios
        r_lm = _project_to_log(res.ratios[i], lmag_table, h, T)
        if r_lm == typemin(T)
            actual_len = i - 1 
            break 
        end
        buffer[i] = r_lm  # ZERO allocations here now!
        curr_lm += r_lm
        max_lm = max(max_lm, curr_lm)
    end

    curr_lm = pref_lm
    curr_sign = pref_sign
    sum_val = curr_sign * exp(curr_lm - max_lm)
    
    @inbounds for i in 1:actual_len
        curr_lm += buffer[i]
        curr_sign *= res.ratios[i].sign
        sum_val += curr_sign * exp(curr_lm - max_lm)
    end

    return exp(max_lm) * sum_val
end

"""
    evaluate_level(res::CycloResult, k::Int, ::Type{T}=Float64; prec=512)

Evaluates the quantum symbol at the exact discrete SU(2)_k root of unity.
"""
function evaluate_level(res::CycloResult, k::Int, ::Type{T}=Float64; prec=512) where {T}
    (res.radical.sign == 0 || res.base_term.sign == 0) && return zero(T)

    #extract type of real
    RealT = T <: Complex ? real(T) : T
    buff = Vector{T}(undef, length(res.ratios))

    val_real = if RealT == BigFloat
        setprecision(BigFloat, prec) do 
            _eval_cresult(res,k,RealT,buff) 
        end
    else
        _eval_cresult(res,k,RealT,buff)
    end
    return T <: Complex ? T(val_real, 0) : T(val_real)
end

#evaluation of CycloMonomial
function _eval_cmonomial(m::CycloMonomial, k::Int,::Type{T}) where {T}
    h = k + 2
    max_d = isempty(m.exps) ? 1 : maximum(keys(m.exps))
    lmag_table = get_phi_table(k, max_d, T)
    lmag = _project_to_log(m, lmag_table, h, T)
    return m.sign * exp(lmag)
end

"""
    evaluate_level(m::CycloMonomial, k::Int, ::Type{T}=Float64; prec=512)
"""
function evaluate_level(m::CycloMonomial, k::Int, ::Type{T}=Float64; prec=512) where {T}
    m.sign == 0 && return zero(T)

    RealT = T <: Complex ? real(T) : T

    val_real = if RealT == BigFloat
        setprecision(BigFloat, prec) do 
            _eval_cmonomial(m, k, RealT) 
        end
    else
        _eval_cmonomial(m, k, T)
    end
    return T <: Complex ? T(val_real, 0) : T(val_real)
end



"""
    cyclo_to_numeric(res::CycloResult, [T=Float64]; k=nothing, q=nothing, theta=nothing, prec=512)

Internal numeric projection engine. Maps the deferred graph into high-performance 
floating-point arithmetic (hardware or arbitrary precision).
"""
function cyclo_to_numeric(res::CycloResult, ::Type{T}=Float64; k=nothing, q=nothing, theta=nothing, prec=512) where {T}
    targets_defined = (!isnothing(k)) + (!isnothing(q)) + (!isnothing(theta))
    if targets_defined != 1
        throw(ArgumentError("Specify exactly one evaluation target: k (integer), theta (real), or q (complex)."))
    end

    if !isnothing(k)
        return evaluate_level(res, Int(k), T; prec=prec)
    elseif !isnothing(theta)
        return evaluate_unit_circle(res, theta, T; prec=prec)
    else
        return evaluate_analytic(res, q, T; prec=prec)
    end
end
