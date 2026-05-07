# ----------------------------------------------------------
#     --- Project DCR to discrete level k ----
#  Evaluate at SU(2)_k roots of unity (q = exp(iπ/(k+2)))
# ----------------------------------------------------------


# Cache for the magnitude table log|Φ_d(q²)|
const MAG_SIEVE_CACHE = Dict{Tuple{DataType, Int, Int}, Any}()
const MAG_SIEVE_LOCK  = ReentrantLock()
const MAG_SIEVE_MAXSIZE = 7_000

"""
    get_mag_table(k::Int, max_d::Int, ::Type{T})
Unified cache for the magnitude table log|Φ_d(q²)|.
"""
function get_mag_table(k::Int, max_d::Int, ::Type{T}) where T

    key_prec = T === BigFloat ? precision(BigFloat) : 53
    key = (T, k, key_prec)

    h   = k + 2

    if haskey(MAG_SIEVE_CACHE, key)
        table = MAG_SIEVE_CACHE[key]::Vector{T}
        # return if it's large enough!
        if length(table) >= max_d
            return table
        end
    end

    return @lock MAG_SIEVE_LOCK begin
        existing = get(MAG_SIEVE_CACHE, key, nothing)

        if existing !== nothing && length(existing) >= max_d
            existing::Vector{T}               
        else
            length(MAG_SIEVE_CACHE) >= MAG_SIEVE_MAXSIZE && empty!(MAG_SIEVE_CACHE)
            new_table = build_mag_table(max(max_d, h), k, T)
            MAG_SIEVE_CACHE[key] = new_table
            new_table
        end
    end::Vector{T}
end


"""
    build_mag_table(D_max::Int, k::Int, ::Type{T})
Computes log|Φ_d(q²)| at q = exp(iπ/(k+2)). 
Uses log-domain Möbius inversion.
"""
function build_mag_table(D_max::Int, k::Int, ::Type{T}) where T
    h = k + 2
    table = Vector{T}(undef, D_max)
    h_inv = inv(T(h))
    
    # Initialize with log|q^2n - 1| = log|2*sin(πn/h)|
    @inbounds for n in 1:D_max
        if n == h # primitive root of unity  
            table[n] = -T(Inf) # pole
        else
            table[n] = log(2 * abs(sinpi(n * h_inv)))
        end
    end
    
    # multiplicative sieve (log-domain Möbius inversion)
    @inbounds for d in 1:D_max
        val = table[d]
        isinf(val) && continue

        for m in (2d):d:D_max
            if !isinf(table[m])
                table[m] -= val
            end
        end
    end
    return table
end




"""
    _log_mag_internal(m::CyclotomicMonomial, table::Vector{T}) where T
Ultra-fast log-magnitude evaluator for the hot loop.
"""
@inline function _log_mag_internal(m::CyclotomicMonomial, table::Vector{T}) where T
    m.sign == 0 && return -T(Inf)
    
    lm = zero(T)
    @inbounds for (d, e) in m.phi_exps
        lm += e * table[d]
    end
    
    return lm
end


#  --- Projection to discrete level k ---- 

"""
    project_discrete(m::CyclotomicMonomial, k::Int, ::Type{T}=Float64)
Evaluates [n]_q or dimensions at SU(2)_k.
Returns the real value.
"""
function project_discrete(m::CyclotomicMonomial, k::Int, ::Type{T}=Float64) where T
    m.sign == 0 && return zero(T)
    h = k + 2
    
    # check zero or pole 
    e_val = _phi_exponent(m, h)
    if e_val > 0
        return zero(T)
    elseif e_val < 0
        throw(DomainError(k, "Topological pole at level k=$k."))
    end
    
    # magnitude
    table = get_mag_table(k, m.max_d, T)
    lm = _log_mag_internal(m, table)
    
    return m.sign * exp(lm)
end



function _evaluate_discrete_dcr(res::DCR, k::Int, ::Type{T}) where {T}
    table = get_mag_table(k, res.max_d, T)

    # project prefactors
    lm_root = _log_mag_internal(res.root, table)
    lm_rad  = _log_mag_internal(res.radical, table)
    lm_base = _log_mag_internal(res.base, table)
    
    pref_lm = lm_root + (0.5 * lm_rad) + lm_base
    pref_sign = res.root.sign * res.base.sign 

    # exit if zero
    if isnan(pref_lm) || (isinf(pref_lm) && pref_lm < 0) || pref_sign == 0
        return zero(T)
    end

    ratios = res.ratios
    num_ratios = length(ratios)
    
    max_lm = pref_lm
    curr_lm = pref_lm
    curr_sign = pref_sign
    res_scaled = T(curr_sign)
    
    # Single-Pass LSE summation
    @inbounds for i in 1:num_ratios
        r_lm = _log_mag_internal(ratios[i], table)
        curr_lm += r_lm
        curr_sign *= ratios[i].sign
        
        if curr_lm > max_lm
            res_scaled = res_scaled * exp(max_lm - curr_lm) + curr_sign
            max_lm = curr_lm
        else
            res_scaled += curr_sign * exp(curr_lm - max_lm)
        end
    end

    return exp(max_lm) * res_scaled
end


"""
    project_discrete(res::DCR, k::Int, ::Type{T}=Float64)
Thread-safe, Zero-Allocation evaluator for full DCR symbols (e.g., 6j, 3j).
"""
function project_discrete(res::DCR, k::Int, ::Type{T}=Float64) where {T}
    h = k + 2
    
    # check zeros and poles
    # prevents NaN results from +Inf and -Inf in the LSE loop
    e_rad  = _phi_exponent(res.radical, h)
    e_root = _phi_exponent(res.root, h)
    e_base = _phi_exponent(res.base, h)
    
    e_net = (e_rad / 2.0) + e_root + e_base
    
    if e_net > 0.0
        return zero(T)
    elseif e_net < 0.0
        throw(DomainError(k, "Topological pole at level k=$k."))
    end
    
    if res.radical.sign == 0 || res.base.sign == 0 || res.root.sign == 0
        return zero(T)
    end
    
    return _evaluate_discrete_dcr(res, k, T)
end