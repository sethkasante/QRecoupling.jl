module QRacahSymbol

# using LinearAlgebra

export _qracah6j, qδtet, logqnfact_table, log_qtri

# ============================================================
# Model + caches
# ============================================================

@inline δ(j1, j2, j3) = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)

# (j3<= j1+j2) && (j1 <= j2 + j3) && (j2 <= j3 + j1) 

@inline qδ(j1, j2, j3, k) = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 

@inline qδtet(j1, j2, j3, j4, j5, j6,k) = qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)
"""
    logqn_table(k::Int, ::Type{T}=Float64) -> Vector{T}

Returns table logqn[n] = log([n]_q) for n = 0..k+1,
with q = exp(iπ/(k+2)).
Index 1 corresponds to n=0.
"""
function logqn_table(k::Int, ::Type{T}=BigFloat) where {T<:AbstractFloat}
    N = k + 2
    logqn = Vector{BigFloat}(undef, N)

    θ = T(pi) / N
    logden = log(sin(θ))

    # symmetry: n ↔ k+2-n
    half = N ÷ 2
    logqn[1] = zero(T)
    @inbounds for n in 1:half
        val = log(sin(n * θ)) - logden
        logqn[N + 1 - n] = logqn[n+1] = val
    end
    return logqn
end
"""
    logqnfact_table(k::Int, ::Type{T}=Float64) -> Vector{T}

Returns table logqnfact[n+1] = log([n]_q!) for n = 0..k+1.
"""
function logqnfact_table(k::Int, ::Type{T}=BigFloat) where {T<:AbstractFloat}
    return cumsum(logqn_table(k, T))
end


@inline function log_qΔ(j1,j2,j3,tab)
    a = j1 + j2 - j3
    b = j1 - j2 + j3
    c = -j1 + j2 + j3
    d = j1 + j2 + j3 
    (tab[a+1] + tab[b+1] + tab[c+1] - tab[d+2])/2
end

@inline function logqtri_coeffs(j1, j2, j3, j4, j5, j6,tab)
    log_qΔ(j1, j2, j3,tab) + 
    log_qΔ(j1, j5, j6,tab) +
    log_qΔ(j2, j4, j6,tab) + 
    log_qΔ(j3, j4, j5,tab)
end

@inline function log_racah_summand(z, α1,α2,α3,α4, β1,β2,β3,tab)
    lognum = tab[z+2]
    logden = tab[z-α1+1] + tab[z-α2+1] + tab[z-α3+1] +
        tab[z-α4+1] + tab[β1-z+1] + tab[β2-z+1] + tab[β3-z+1]
    return lognum - logden
end

function _qracah6j(j1, j2, j3, j4, j5, j6,k)
    if !qδtet(j1, j2, j3, j4, j5, j6,k) 
        return zero(T)
    end
    table = logqnfact_table(k)
    T = logqtri_coeffs(j1, j2, j3, j4, j5, j6,table)

    α1, α2, α3, α4 = j1 + j2 + j3, j1 + j5 + j6, j2 + j4 + j6, j3 + j4 + j5
    β1, β2, β3 = j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6

    zrange = max(α1,α2,α3,α4):min(β1,β2,β3)
    res = big(0)
    
    for z in zrange
        logTsz = T + log_racah_summand(z,α1,α2,α3,α4,β1,β2,β3,table) 
        res += iseven(z) ? exp(logTsz) :  -exp(logTsz) 
    end
    Float64(res)
end

end #module 