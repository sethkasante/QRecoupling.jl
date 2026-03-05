# src/Numerics.jl

# ============================================================
# Caches & Constants
# ============================================================

# Precomputing BigFloat Pi 
const BIG_PI = big(π)

# Cache for precomputed log-q-factorials: key (k, prec) -> Vector{BigFloat}
const LOGQFACT_CACHE = LRU{Tuple{Int,Int}, Vector{BigFloat}}(maxsize = 10240)

"""
    NumericSU2kModel(k::Int, prec::Int=256) -> NumericSU2kModel

Creates a NumericSU2kModel with level k, computing the log-q-factorial table 
with BigFloat precision. Uses LRU cache to avoid recomputation.
"""
function NumericSU2kModel(k::Int, prec::Int=256)
    tab = get!(LOGQFACT_CACHE, (k, prec)) do
        logqnfact_table(k, prec) 
    end
    return NumericSU2kModel(k, tab)
end

# ====================
# Log-Space Tables
# ====================

"""
    logqn_table(k::Int, prec::Int=256) -> Vector{BigFloat}

Computes table of log q-integers: logqn[n+1] = log([n]_q) for n = 0..k+1.
"""
function logqn_table(k::Int, prec::Int=256)::Vector{BigFloat}
    N = k + 2
    setprecision(BigFloat, prec) do
        θ = BIG_PI / BigFloat(N)
        logsinθ = log(sin(θ))

        # Allocate once
        tab = Vector{BigFloat}(undef, N)
        half = N ÷ 2
        tab[1] = zero(BigFloat) 

        @inbounds for n in 1:half
            val = log(sin(n * θ)) - logsinθ
            tab[n+1] = val
            tab[N + 1 - n] = val # Mirror symmetry
        end
        return tab
    end
end

"""
    logqnfact_table(k::Int, prec::Int=256) -> Vector{BigFloat}

Computes table of log quantum factorials: logqnfact[n+1] = log([n]_q!).
"""
function logqnfact_table(k::Int, prec::Int=256)::Vector{BigFloat}
    logqn = logqn_table(k, prec)
    
    setprecision(BigFloat, prec) do
        # We only need up to k+2 because z is capped at k.
        tab = Vector{BigFloat}(undef, k + 2) 
        tab[1] = logqn[1]  
        
        @inbounds for n in 2:k+2
            tab[n] = tab[n-1] + logqn[n]
        end
        return tab
    end
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