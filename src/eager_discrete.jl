

# ---------------------------------------------------------------------------
# :discrete eager builder - uses Log-Sum-Exp (LSE) for Hypergeometric Summation
# Direct numerical computation of the Racah-Wigner 3j and 6j symbols 
# ---------------------------------------------------------------------------

# cache for internal use.. make thread safe
const _NUMERIC_MODEL_CACHE = Dict{Tuple{Int, DataType, Int}, Any}()
const _NUMERIC_MODEL_LOCK  = ReentrantLock()


#----- Model Construction ------- 

"""
    NumericSU2kModel{T<:AbstractFloat}

A highly optimized memory cache holding the precomputed log-q-factorial table 
for a specific level `k`. Reusing this model across multiple symbol evaluations 
bypasses redundant O(k) trigonometric initializations.
"""
struct NumericSU2kModel{T<:AbstractFloat}
    k::Int
    logqnfact::Vector{T}
end


function NumericSU2kModel(k::Int, T::Type{<:AbstractFloat}=Float64)
    k ≥ 1 || throw(ArgumentError("level k must be ≥ 1, got $k"))
    
    key_prec = T === BigFloat ? precision(BigFloat) : 53 
    key = (k, T, key_prec)
 
    tab = @lock _NUMERIC_MODEL_LOCK begin
        # evict if about to add a new key
        if !haskey(_NUMERIC_MODEL_CACHE, key)
            length(_NUMERIC_MODEL_CACHE) >= 5_000 && empty!(_NUMERIC_MODEL_CACHE)
            _NUMERIC_MODEL_CACHE[key] = build_logqnfact_table(k, T)
        end
        _NUMERIC_MODEL_CACHE[key]
    end

    return NumericSU2kModel{T}(k, tab)
end


function build_logqnfact_table(k::Int, T::Type{<:AbstractFloat})::Vector{T}
    N = k + 2 
    θ = one(T) / T(N)
    logsinθ = log(sinpi(θ))
    
    # compute log of qinteger: log([n]_q)
    logqn = Vector{T}(undef, N)
    logqn[1] = zero(T) 

    # use symmetries of sin
    half = N ÷ 2 
    @inbounds for n in 1:half
        val = log(sinpi(T(n) * θ)) - logsinθ
        logqn[n+1] = val
        logqn[N + 1 - n] = val
    end
    
    # log([n]_q!)
    tab = Vector{T}(undef, N) 
    tab[1] = zero(T) # log([0]_q!) = 0
    
    @inbounds for n in 2:N
        tab[n] = tab[n-1] + logqn[n]
    end
    
    return tab
end


# --- quantum recouplings symbols --- 


# Computes log(Δ) for admissible triangular coefficients.. using doubled Spins (J = 2j)
@inline function log_qΔ(J1::Int, J2::Int, J3::Int, tab::Vector{T})::T where {T}
    a = (J1 + J2 - J3) ÷ 2
    b = (J1 - J2 + J3) ÷ 2
    c = (-J1 + J2 + J3) ÷ 2
    d = (J1 + J2 + J3) ÷ 2
    return (tab[a+1] + tab[b+1] + tab[c+1] - tab[d+2]) / 2 
end



#---- Compute quantum 6j Symbol ---- 

"""
    q6j_direct(J1, J2, J3, J4, J5, J6, k, T)

Evaluates the quantum 6j-symbol using a direct Log-Sum-Exp (LSE) alternating summation.
Inputs must be doubled spins (J = 2j ∈ ℤ).
"""
function q6j_direct(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int, k::Int, T::Type{<:AbstractFloat})
    !qδtet(J1, J2, J3, J4, J5, J6, k) && return zero(T)
    
    model = NumericSU2kModel(k, T)
    table = model.logqnfact
    
    # Prefactor
    logT = log_qΔ(J1, J2, J3, table) + log_qΔ(J1, J5, J6, table) + 
           log_qΔ(J2, J4, J6, table) + log_qΔ(J3, J4, J5, table)

    # Bounds 
    α1 = (J1 + J2 + J3) ÷ 2; α2 = (J1 + J5 + J6) ÷ 2; α3 = (J2 + J4 + J6) ÷ 2; α4 = (J3 + J4 + J5) ÷ 2
    β1 = (J1 + J2 + J4 + J5) ÷ 2; β2 = (J1 + J3 + J4 + J6) ÷ 2; β3 = (J2 + J3 + J5 + J6) ÷ 2

    z_min = max(α1, α2, α3, α4)
    z_max = min(β1, β2, β3, model.k) 
    
    z_min > z_max && return zero(T)

    logmax = typemin(T)
    res_scaled = zero(T)
    curr_sign = iseven(z_min) ? one(T) : -one(T)
    
    @inbounds for z in z_min:z_max
        log_sz = table[z+2] - (table[z-α1+1] + table[z-α2+1] + table[z-α3+1] + table[z-α4+1] + table[β1-z+1] + table[β2-z+1] + table[β3-z+1])
        
        if log_sz > logmax
            # A new maximum was found. Scale the running sum down.
            res_scaled = res_scaled * exp(logmax - log_sz) + curr_sign
            logmax = log_sz
        else
            res_scaled += curr_sign * exp(log_sz - logmax)
        end
        
        curr_sign = -curr_sign
    end

    return exp(logT + logmax) * res_scaled
end


#---- Compute quantum 3j symbol ---- 

"""
    q3j_direct(J1, J2, J3, M1, M2, M3, k, T)
Inputs must be doubled spins (J = 2j ∈ ℤ) and doubled magnetic projections (M = 2m ∈ ℤ).
"""
function q3j_direct(J1::Int, J2::Int, J3::Int, M1::Int, M2::Int, M3::Int, k::Int, T::Type{<:AbstractFloat})
    (!qδ(J1, J2, J3, k) || M1 + M2 + M3 != 0) && return zero(T)
    
    model = NumericSU2kModel(k, T)
    table = model.logqnfact
    
    # prefactor
    log_delta = log_qΔ(J1, J2, J3, table)
    log_facts = table[(J1 + M1) ÷ 2 + 1] + table[(J1 - M1) ÷ 2 + 1] + 
                table[(J2 + M2) ÷ 2 + 1] + table[(J2 - M2) ÷ 2 + 1] + 
                table[(J3 + M3) ÷ 2 + 1] + table[(J3 - M3) ÷ 2 + 1]
    
    logT = log_delta + 0.5 * log_facts

    # bounds
    α1 = (J3 - J2 + M1) ÷ 2; α2 = (J3 - J1 - M2) ÷ 2
    β1 = (J1 + J2 - J3) ÷ 2; β2 = (J1 - M1) ÷ 2; β3 = (J2 + M2) ÷ 2

    z_min = max(-α1, -α2, 0)
    z_max = min(β1, β2, β3, model.k) 
    
    z_min > z_max && return zero(T)

    logmax = typemin(T)
    res_scaled = zero(T)
    phase_offset = α1 - α2 
    curr_sign = isodd(z_min + phase_offset) ? -one(T) : one(T)

    @inbounds for z in z_min:z_max
        log_sz = -(table[z+1] + table[α1+z+1] + table[α2+z+1] + table[β1-z+1] + table[β2-z+1] + table[β3-z+1])
        
        if log_sz > logmax
            res_scaled = res_scaled * exp(logmax - log_sz) + curr_sign
            logmax = log_sz
        else
            res_scaled += curr_sign * exp(log_sz - logmax)
        end
        
        curr_sign = -curr_sign
    end

    return exp(logT + logmax) * res_scaled
end