
#Numerics.jl

# Cache for precomputed log-q-factorials: key (k, DataType) -> Vector{BigFloat}
const LOGQFACT_CACHE = LRU{Tuple{Int,Int}, Vector{BigFloat}}(maxsize = 10240)

"""
    NumericSU2kModel(k::Int) -> NumericSU2kModel

Create an NumericSU2kModel with level k, computing log-q-factorial table with BigFloat precision.
Uses LRU cache to avoid recomputation.
"""
function NumericSU2kModel(k::Int, prec=256)
    tab = get!(LOGQFACT_CACHE, (k, prec)) do
        logqnfact_table(k, prec) 
    end
    return NumericSU2kModel(k, tab)
end

"""
    logqn_table(k::Int, ::Type{T}=BigFloat) -> Vector{T}

Compute table of log q-integers: logqn[n+1] = log([n]_q) for n = 0..k+1,
where q = exp(2iπ/(k+2)).

Uses trigonometric symmetry to compute only half the values.
"""
function logqn_table(k::Int, prec::Int=256)::Vector{BigFloat}
    N = k + 2
    setprecision(prec) do
        θ = big(pi) / BigFloat(N)
        logsinθ = log(sin(θ))

        tab = Vector{BigFloat}(undef, N)
        half = N ÷ 2
        tab[1] = zero(BigFloat) 

        @inbounds for n in 1:half
            val = log(sin(n*θ)) - logsinθ
            tab[N + 1 - n] = tab[n+1] = val
        end
        return tab
    end
end

"""
    logqnfact_table(k::Int, ::Type{T}=BigFloat) -> Vector{T}

Compute table of log quantum factorials: logqnfact[n+1] = log([n]_q!) for n = 0..k+1.
"""
# function logqnfact_table(k::Int, ::Type{T}=BigFloat)::Vector{BigFloat} where {T<:AbstractFloat}
#     return cumsum(logqn_table(k, T))
# end
function logqnfact_table(k::Int, prec::Int=256)::Vector{BigFloat}
    logqn = logqn_table(k, prec)
    # setprecision(prec) do
    tab = Vector{BigFloat}(undef, k+2) 
    tab[1] = logqn[1]  
    @inbounds for n in 2:k+2
        tab[n] = tab[n-1] + logqn[n]
    end
    # tab[k+3] = -Inf # Include [k+2]! = 0
    return tab
    # end
end



# ==== core log building blocks ======#

"""
    log_qΔ(j1, j2, j3, tab) -> T

Compute log of quantum triangle coefficient qΔ(j1,j2,j3) using log-factorial table.
"""
@inline function log_qΔ(j1, j2, j3, tab)::BigFloat
    a = Int(j1 + j2 - j3)
    b = Int(j1 - j2 + j3)
    c = Int(-j1 + j2 + j3)
    d = Int(j1 + j2 + j3)
    (tab[a+1] + tab[b+1] + tab[c+1] - tab[d+2]) / 2
end

"""
    logqtri_coeffs(j1, j2, j3, j4, j5, j6, tab) -> T

Sum of log quantum triangle coefficients for all four triangles of a tetrahedron.
"""
@inline function logqtri_coeffs(j1, j2, j3, j4, j5, j6, tab)::BigFloat
    log_qΔ(j1, j2, j3, tab) + 
    log_qΔ(j1, j5, j6, tab) +
    log_qΔ(j2, j4, j6, tab) + 
    log_qΔ(j3, j4, j5, tab)
end

"""
    log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, tab) -> T

Compute log of summand for given z in the Racah sum.
"""
@inline function log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, tab)::BigFloat
    lognum = tab[z+2]
    logden = tab[z-α1+1] + tab[z-α2+1] + tab[z-α3+1] +
        tab[z-α4+1] + tab[β1-z+1] + tab[β2-z+1] + tab[β3-z+1]
    return lognum - logden
end

"""
    _qracah6j(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6) -> BigFloat

Compute quantum 6j symbol using Racah summation formula.
Returns 0 if admissibility conditions are not met.
"""
function _qracah6j(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return zero(BigFloat)
    end

    table = model.logqnfact

    logT = logqtri_coeffs(j1, j2, j3, j4, j5, j6, table)

    α1 = Int(j1 + j2 + j3) 
    α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6) 
    α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5)
    β2 = Int(j1 + j3 + j4 + j6) 
    β3 = Int(j2 + j3 + j5 + j6)

    zrange = max(α1, α2, α3, α4):min(min(β1, β2, β3),model.k) #take care of [k+2]! terms!! z has to be less than k  
    res = zero(BigFloat)
    @inbounds for z in zrange # loop up to [k+1] terms only 
        logTsz = logT + log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, table) 
        res += iseven(z) ? exp(logTsz) : -exp(logTsz) 
    end
    return res
end


"""
    qracah6j(j1, j2, j3, j4, j5, j6, k::Int) -> BigFloat

Compute quantum 6j symbol {j1 j2 j3; j4 j5 j6} at level k.

Creates a temporary NumericSU2kModel; use _qracah6j directly for better performance 
when computing multiple symbols with the same k.
"""
qracah6j(j1, j2, j3, j4, j5, j6, k::Int) =
    _qracah6j(NumericSU2kModel(k), j1, j2, j3, j4, j5, j6)



qracah6j_numeric(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6) =
    _qracah6j(model, j1, j2, j3, j4, j5, j6)


#TODO: Add more functions Rsymbols, Fsymbols,Gsymbols, qclebschgordan, q3jsymbols 















# src/Numerics.jl

# Parameterize the model to accept Float64 OR BigFloat
struct NumericSU2kModel{T <: AbstractFloat}
    k::Int
    logqnfact::Vector{T}
end

# Precomputing BigFloat Pi 
const BIG_PI = big(π)

# Cache for precomputed log-q-factorials: key (k, prec) -> Vector{BigFloat}
const LOGQFACT_CACHE = LRU{Tuple{Int,Int}, Vector{BigFloat}}(maxsize = 10240)

# Default to Float64 for blistering speed, opt-in to BigFloat for extreme precision
function NumericSU2kModel(k::Int; T::Type{<:AbstractFloat}=Float64, prec::Int=256)
    tab = get!(LOGQFACT_CACHE, (k, T, prec)) do
        logqnfact_table(k, T, prec) 
    end
    return NumericSU2kModel{T}(k, tab)
end

function logqn_table(k::Int, T::Type, prec::Int)::Vector{T}
    N = k + 2
    
    # If using Float64, bypass setprecision
    if T === Float64
        θ = pi / N
        logsinθ = log(sin(θ))
        tab = Vector{T}(undef, N)
        tab[1] = zero(T)
        half = N ÷ 2
        @inbounds for n in 1:half
            val = log(sin(n * θ)) - logsinθ
            tab[n+1] = val
            tab[N + 1 - n] = val
        end
        return tab
    else
        setprecision(BigFloat, prec) do
            θ = BIG_PI / BigFloat(N)
            # ... identical to current BigFloat logic ...
            logsinθ = log(sin(θ))
            tab = Vector{T}(undef, N)
            tab[1] = zero(T)
            half = N ÷ 2
            @inbounds for n in 1:half
                val = log(sin(n * θ)) - logsinθ
                tab[n+1] = val
                tab[N + 1 - n] = val
            end
            return tab
        end
    end
end

function logqnfact_table(k::Int, T::Type, prec::Int)::Vector{T}
    logqn = logqn_table(k, T, prec)
    tab = Vector{T}(undef, k + 2) 
    tab[1] = logqn[1]  
    @inbounds for n in 2:k+2
        tab[n] = tab[n-1] + logqn[n]
    end
    return tab
end



# =================================
# Core Log Building Blocks
# =================================

@inline function log_qΔ(j1, j2, j3, tab)::BigFloat
    a = Int(j1 + j2 - j3)
    b = Int(j1 - j2 + j3)
    c = Int(-j1 + j2 + j3)
    d = Int(j1 + j2 + j3)
    return 0.5 * (tab[a+1] + tab[b+1] + tab[c+1] - tab[d+2])
end

@inline function logqtri_coeffs(j1, j2, j3, j4, j5, j6, tab)::BigFloat
    return log_qΔ(j1, j2, j3, tab) + 
           log_qΔ(j1, j5, j6, tab) +
           log_qΔ(j2, j4, j6, tab) + 
           log_qΔ(j3, j4, j5, tab)
end

@inline function log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, tab)::BigFloat
    lognum = tab[z+2]
    logden = tab[z-α1+1] + tab[z-α2+1] + tab[z-α3+1] +
             tab[z-α4+1] + tab[β1-z+1] + tab[β2-z+1] + tab[β3-z+1]
    return lognum - logden
end

# ===========================================
# High-Performance Numeric 6j Evaluator
# ===========================================

"""
    _qracah6j_stable(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6) -> BigFloat

Computes the quantum 6j symbol using a log-sum-exp shift to guarantee 
floating-point stability even at massive values of k.
"""
function _qracah6j_stable(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return zero(BigFloat)
    end
    
    table = model.logqnfact

    logT = logqtri_coeffs(j1, j2, j3, j4, j5, j6, table)

    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)

    # Capping z at model.k gracefully bypasses the [k+2]_q! = 0 wall.
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3, model.k) 
    
    # 1. Find the maximum exponent for log-sum-exp stabilization
    logmax = -Inf
    @inbounds for z in zrange
        logTsz = logT + log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, table)
        logmax = max(logmax, logTsz)
    end

    # 2. Compute the stable alternating sum
    res_scaled = zero(BigFloat)
    @inbounds for z in zrange
        logTsz = logT + log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, table)
        # Shift by logmax prevents exp() from hitting Infinity or flushing to zero
        val = exp(logTsz - logmax)
        res_scaled += iseven(z) ? val : -val
    end

    return exp(logmax) * res_scaled
end

qracah6j_numeric(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6) =
    _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)

#TODO: define qracah numeric for input k. 

# ============================================================
# Extension Suite: TQFT Functions (Numeric)
# ============================================================

"""
    qdim_numeric(j, model::NumericSU2kModel)

Quantum dimension [2j+1]_q using the log-tables.
"""
function qdim_numeric(j, model::NumericSU2kModel)
    n = Int(2j + 1)
    # [n]_q = [n]_q! / [n-1]_q!
    # exp(log([n]!) - log([n-1]!))
    return exp(model.logqnfact[n+1] - model.logqnfact[n])
end

"""
    fsymbol_numeric(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)

Returns the normalized F-matrix element.
F = (-1)^{j1+j2+j3+j4} √([2j3+1][2j6+1]) {j1 j2 j3; j4 j5 j6}_q
"""
function fsymbol_numeric(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)
    val_6j = _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
    
    d3 = qdim_numeric(j3, model)
    d6 = qdim_numeric(j6, model)
    
    phase = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    return phase * sqrt(d3 * d6) * val_6j
end

"""
    rmatrix_numeric(j1, j2, j3, k::Int)

Returns the braiding R-symbol. Since this is purely a phase, it is evaluated 
directly without the log-tables.
"""
function rmatrix_numeric(j1, j2, j3, k::Int)
    phase_exp = Int(j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1
    return s * cispi(phase_exp / (k + 2))
end

"""
    gsymbol_numeric(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)

The G-symbol (Tetrahedral Weight) used in state sums.
G = {6j} * √(Π [2j_i + 1]_q)
"""
function gsymbol_numeric(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)
    val_6j = _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
    
    dims_prod = qdim_numeric(j1, model) * qdim_numeric(j2, model) * qdim_numeric(j3, model) * qdim_numeric(j4, model) * qdim_numeric(j5, model) * qdim_numeric(j6, model)
                
    return val_6j * sqrt(dims_prod)
end