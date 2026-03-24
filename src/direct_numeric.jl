
# ---------------------------------------------------------------------------
# :numeric builder - uses (Log-Sum-Exp for Hypergeometric Summation)
# `Brute force` numerical computation of the Racah-Wigner 3j and 6j symbols 
# ---------------------------------------------------------------------------


# LRU Cached tables for computations
const LOGQFACT_CACHE    = LRU{Tuple{Int, DataType, Int}, Any}(maxsize = 4096)
# const Q6J_NUMERIC_CACHE = LRU{Tuple{NTuple{6, Float64}, Int, Int}, Any}(maxsize=50_000)


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


# Constructs or retrieves a cached `NumericSU2kModel` for level `k`. 
function NumericSU2kModel(k::Int; T::Type{<:AbstractFloat}=Float64, prec::Int=256)
    key_prec = T == Float64 ? 53 : prec 
    tab = get!(LOGQFACT_CACHE, (k, T, key_prec)) do
        logqnfact_table(k, T, key_prec) 
    end
    return NumericSU2kModel{T}(k, tab)
end

"""
    logqnfact_table(k::Int, T::Type, prec::Int) -> Vector{T}

Generates the exact logarithmic quantum factorial table for n ∈ [0, k+1].
Leverages the reflection symmetry of the sine function to halve the required 
trigonometric computations.
"""
function logqnfact_table(k::Int, T::Type{<:AbstractFloat}, prec::Int)::Vector{T}
    # enforce precision when using BigFloats
    if T === BigFloat
        return setprecision(prec) do
            build_logqnfact_table(k, T)
        end
    else
        return build_logqnfact_table(k, T)
    end
end

# build a table for log q-factorials for n=0,...,k+1.
function build_logqnfact_table(k::Int, T::Type{<:AbstractFloat})::Vector{T}
    N = k + 2 
    θ = one(T) / T(N)
    logsinθ = log(sinpi(θ))
    
    # compute log of qinteger: log([n]_q)
    logqn = Vector{T}(undef, N)
    logqn[1] = zero(T) 

    #use symmetries of sin
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


# Computes log(Δ) for a valid triangular coefficients.
@inline function log_qΔ(j1::Spin, j2::Spin, j3::Spin, tab::Vector{T})::T where {T}
    a = Int(j1 + j2 - j3)
    b = Int(j1 - j2 + j3)
    c = Int(-j1 + j2 + j3)
    d = Int(j1 + j2 + j3)
    
    return (tab[a+1] + tab[b+1] + tab[c+1] - tab[d+2]) / 2 
end


#---- Compute quantum 3j ssymbol ---- 

"""
    _q3j_stable(model::NumericSU2kModel, j1, j2, j3, m1, m2)

Evaluates the quantum 3j-symbol using a NaN-safe Log-Sum-Exp alternating summation.
"""
function _q3j_stable(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin)::T where {T}
    table = model.logqnfact
    
    # Prefactor
    log_delta = log_qΔ(j1, j2, j3, table)
    log_facts = table[Int(j1+m1)+1] + table[Int(j1-m1)+1] + 
                table[Int(j2+m2)+1] + table[Int(j2-m2)+1] + 
                table[Int(j3-m1-m2)+1] + table[Int(j3+m1+m2)+1]
    
    logT = log_delta + 0.5 * log_facts

    # Bounds
    α1 = Int(j3 - j2 + m1) 
    α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 - j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    z_min = max(-α1, -α2, 0)
    z_max = min(β1, β2, β3, model.k) 
    
    z_min > z_max && return zero(T)

    # find maximum of log
    logmax = typemin(T)
    @inbounds for z in z_min:z_max
        log_sz = -(table[z+1] + table[α1+z+1] + table[α2+z+1] + table[β1-z+1] + table[β2-z+1] + table[β3-z+1])
        logmax = max(logmax, log_sz)
    end

    # perform summation
    res_scaled = zero(T)
    phase_offset = α1 - α2 
    curr_sign = isodd(z_min + phase_offset) ? -one(T) : one(T)
    
    @inbounds for z in z_min:z_max
        log_sz = -(table[z+1] + table[α1+z+1] + table[α2+z+1] + table[β1-z+1] + table[β2-z+1] + table[β3-z+1])
        res_scaled += curr_sign * exp(log_sz - logmax)
        curr_sign = -curr_sign
    end

    return exp(logT + logmax) * res_scaled
end



#---- Compute quantum 6j Symbol ---- 

"""
    _q6j_stable(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)

Evaluates the quantum 6j-symbol using a direct Log-Sum-Exp (LSE) alternating summation.

# Performance Note:
This function is optimized for single-shot floating-point evaluations. 
It bypasses all abstract cyclotomic algebra, jumping straight into a pre-allocated 
CPU register hot-loop. It is mathematically immune to `NaN` and `Inf` overflows 
at large levels k.
"""
function _q6j_stable(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::T where {T}
    table = model.logqnfact
    
    # Prefactor
    logT = log_qΔ(j1, j2, j3, table) + log_qΔ(j1, j5, j6, table) + 
           log_qΔ(j2, j4, j6, table) + log_qΔ(j3, j4, j5, table)

    # Bounds
    α1 = Int(j1 + j2 + j3) 
    α2 = Int(j1 + j5 + j6)
    α3 = Int(j2 + j4 + j6)
    α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5)
    β2 = Int(j1 + j3 + j4 + j6)
    β3 = Int(j2 + j3 + j5 + j6)

    z_min = max(α1, α2, α3, α4)
    z_max = min(β1, β2, β3, model.k) # note: [k+2]!=0
    
    z_min > z_max && return zero(T)
    
    # find maximum of log
    logmax = typemin(T)
    @inbounds for z in z_min:z_max
        log_sz = table[z+2] - (table[z-α1+1] + table[z-α2+1] + table[z-α3+1] + table[z-α4+1] + table[β1-z+1] + table[β2-z+1] + table[β3-z+1])
        logmax = max(logmax, log_sz)
    end

    # perform alternating sum
    res_scaled = zero(T)
    curr_sign = iseven(z_min) ? one(T) : -one(T)
    
    @inbounds for z in z_min:z_max
        log_sz = table[z+2] - (table[z-α1+1] + table[z-α2+1] + table[z-α3+1] + table[z-α4+1] + table[β1-z+1] + table[β2-z+1] + table[β3-z+1])
        res_scaled += curr_sign * exp(log_sz - logmax)
        curr_sign = -curr_sign
    end

    return exp(logT + logmax) * res_scaled
end


#----- Internal dispatchers ----

q3j_direct(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, k::Int; T::Type{<:AbstractFloat}=Float64, prec::Int=256) = 
    _q3j_stable(NumericSU2kModel(k; T=T, prec=prec), j1, j2, j3, m1, m2)


"""
    q6j_direct(j1, j2, j3, j4, j5, j6, k; T=Float64, prec=256)

The public-facing fast numeric API for the quantum 6j-symbol. 
Automatically manages the `NumericSU2kModel` caching behind the scenes.
"""   
q6j_direct(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int; T::Type{<:AbstractFloat}=Float64, prec::Int=256) = 
    _q6j_stable(NumericSU2kModel(k; T=T, prec=prec), j1, j2, j3, j4, j5, j6)

export q6j_direct, q3j_direct