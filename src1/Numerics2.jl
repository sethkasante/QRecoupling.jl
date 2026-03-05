
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