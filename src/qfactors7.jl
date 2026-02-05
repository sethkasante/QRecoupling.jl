module QuantumSU2kFast

# using LinearAlgebra

# export SU2kModel,
#        q6j,
#        admissible_triple,
#        admissible_tet

# # ============================================================
# # Model + caches
# # ============================================================

# δ(j1, j2, j3) = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)

# # (j3<= j1+j2) && (j1 <= j2 + j3) && (j2 <= j3 + j1) 

# qδ(j1, j2, j3, k) = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 

# qδtet(j1, j2, j3, j4, j5, j6,k) = qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)
# """
#     logqn_table(k::Int, ::Type{T}=Float64) -> Vector{T}

# Returns table logqn[n] = log([n]_q) for n = 0..k+1,
# with q = exp(iπ/(k+2)).
# Index 1 corresponds to n=0.
# """
# function logqn_table(k::Int, ::Type{T}=Float64) where {T<:AbstractFloat}
#     N = k + 2
#     logqn = Vector{T}(undef, N)

#     θ = one(T) / N
#     logden = log(sinpi(θ))

#     # symmetry: n ↔ k+2-n
#     half = N ÷ 2
#     logqn[1] = zero(T)

#     @inbounds for n in 1:half
#         val = log(sinpi(n * θ)) - logden
#         logqn[N + 1 - n] = logqn[n+1] = val
#         # if n <= N && k + 2 - n <= N
#     end
#     return logqn
# end
# """
#     logqnfact_table(k::Int, ::Type{T}=Float64) -> Vector{T}

# Returns table logqnfact[n+1] = log([n]_q!) for n = 0..k+1.
# """
# function logqnfact_table(k::Int, ::Type{T}=Float64) where {T<:AbstractFloat}
#     return cumsum(logqn_table(k, T))
# end

# function log_qΔ(j1,j2,j3,table)
#     a = j1 + j2 - j3
#     b = j1 - j2 + j3
#     c = -j1 + j2 + j3
#     d = j1 + j2 + j3 + 1
#     (table[a+1] + table[b+1] + table[c+1] - table[d+1])/2
# end

# function logqtri_coeffs(j1, j2, j3, j4, j5, j6,table)
#     log_qΔ(j1,j2,j3,table) + log_qΔ(j1, j5, j6,table) +
#         log_qΔ(j2, j4, j6,table) + log_qΔ(j3, j4, j5,table)
# end

# function log_racah_summand(z, α1,α2,α3,α4, β1,β2,β3,table)
#     log_num = table[z + 2]
#     log_den = table[z - α1 + 1] + table[z-α2+1] + table[z-α3+1] +
#         table[z-α4+1] + table[β1-z+1] + table[β2-z+1] + table[β3-z+1]
#     return log_num - log_den
# end

# function _qracah6j(j1, j2, j3, j4, j5, j6,k)
#     if !qδtet(j1, j2, j3, j4, j5, j6,k) 
#         return zero(T)
#     end
#     table = logqnfact_table(k)
#     T = logqtri_coeffs(j1, j2, j3, j4, j5, j6,table)

#     α1, α2, α3, α4 = j1 + j2 + j3, j1 + j5 + j6, j2 + j4 + j6, j3 + j4 + j5
#     β1, β2, β3 = j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6

#     zrange = max(α1,α2,α3,α4):min(β1,β2,β3)
#     res = big(0)
#     for z in zrange
#         logTsz = T + log_racah_summand(z,α1,α2,α3,α4,β1,β2,β3,table) 
#         res += iseven(z) ? exp(logTsz) :  -exp(logTsz) 
#     end
#     res
# end


export q6j, log_racah_term, logqnfact_table

# ------------------------------------------------------------
# Admissibility
# ------------------------------------------------------------

@inline δ(j1,j2,j3) =
    isinteger(j1 + j2 + j3) && abs(j1 - j2) ≤ j3 ≤ j1 + j2

@inline qδ(j1,j2,j3,k) =
    δ(j1,j2,j3) && (j1 + j2 + j3 ≤ k)

@inline qδtet(j1,j2,j3,j4,j5,j6,k) =
    qδ(j1,j2,j3,k) &&
    qδ(j1,j5,j6,k) &&
    qδ(j2,j4,j6,k) &&
    qδ(j3,j4,j5,k)

# ------------------------------------------------------------
# log([n]_q) and log([n]_q!) tables
# ------------------------------------------------------------

"""
    logqn_table(k, T)

Returns logqn[n+1] = log([n]_q), n = 0..k+1
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
        v = log(sin(n * θ)) - logden
        logqn[n+1] = v
        logqn[N+1-n] = v
    end
    return logqn
end

"""
    logqnfact_table(k, T)

Returns logqnfact[n+1] = log([n]_q!), n = 0..k+1
"""
function logqnfact_table(k::Int, ::Type{T}) where {T<:AbstractFloat}
    logqn = logqn_table(k, T)
    return cumsum(logqn)
end

# ------------------------------------------------------------
# Triangle prefactor
# ------------------------------------------------------------

@inline function log_qΔ(j1,j2,j3,tab)
    a = j1 + j2 - j3
    b = j1 - j2 + j3
    c = -j1 + j2 + j3
    d = j1 + j2 + j3 + 1
    return (tab[a+1] + tab[b+1] + tab[c+1] - tab[d+1]) / 2
end

@inline function log_triangle_prefactor(j1,j2,j3,j4,j5,j6,tab)
    log_qΔ(j1,j2,j3,tab) +
    log_qΔ(j1,j5,j6,tab) +
    log_qΔ(j2,j4,j6,tab) +
    log_qΔ(j3,j4,j5,tab)
end

# ------------------------------------------------------------
# Racah summand
# ------------------------------------------------------------

@inline function log_racah_term(z, α, β, tab)
    lognum = tab[z + 2]
    logden =
        tab[z-α[1]+1] + tab[z-α[2]+1] + tab[z-α[3]+1] + tab[z-α[4]+1] +
        tab[β[1]-z+1] + tab[β[2]-z+1] + tab[β[3]-z+1]
    return lognum - logden
end

# ------------------------------------------------------------
# Stable log-alternating-sum
# ------------------------------------------------------------

function alternating_logsumexp(logvals::Vector{T}, signs::Vector{Int}) where {T}
    m = maximum(logvals)
    s = zero(T)
    @inbounds for i in eachindex(logvals)
        s += signs[i] * exp(logvals[i] - m)
    end
    return exp(m) * s
end

# ------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------

function q6j(j1,j2,j3,j4,j5,j6,k; T=Float64)
    if !qδtet(j1,j2,j3,j4,j5,j6,k)
        return zero(T)
    end

    tab = logqnfact_table(k, T)
    logT = log_triangle_prefactor(j1,j2,j3,j4,j5,j6,tab)

    α = (
        j1+j2+j3,
        j1+j5+j6,
        j2+j4+j6,
        j3+j4+j5
    )
    β = (
        j1+j2+j4+j5,
        j1+j3+j4+j6,
        j2+j3+j5+j6
    )

    zmin = maximum(α)
    zmax = minimum(β)

    nz = zmax - zmin + 1
    logvals = Vector{T}(undef, nz)
    signs   = Vector{Int}(undef, nz)

    @inbounds for (i,z) in enumerate(zmin:zmax)
        logvals[i] = logT + log_racah_term(z, α, β, tab)
        signs[i] = iseven(z) ? 1 : -1
    end

    return alternating_logsumexp(logvals, signs)
end

end # module
