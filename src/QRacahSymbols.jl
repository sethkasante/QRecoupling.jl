
module QRacahSymbol

# QRacahSymbol.jl
using LRUCache

export qracah6j, 
        qδ,
        qdim, 
        SU2kModel,
        rsymbol, 
        fsymbol,
        gsymbol

# ================
# Model + caches
# ================

"""
    SU2kModel{T<:AbstractFloat}

A model for quantum SU(2)_k angular momentum algebra.

# Fields
- `k::Int`: The level parameter
- `logqnfact::Vector{T}`: Precomputed log of quantum factorials
"""
struct SU2kModel{T<:AbstractFloat}
    k::Int
    logqnfact::Vector{T}
end

"""
    SU2kModel(k::Int) -> SU2kModel

Create an SU2kModel with level k, computing log-q-factorial table with BigFloat precision.
Uses LRU cache to avoid recomputation.
"""
function SU2kModel(k::Int)::SU2kModel
    BigT = BigFloat
    tab = get!(LOGQFACT_CACHE, (k, BigT)) do
        logqnfact_table(k, BigT)
    end
    return SU2kModel{BigT}(k, tab)
end

# Cache for precomputed log-q-factorials: key (k, DataType) -> Vector{BigFloat}
const LOGQFACT_CACHE = LRU{Tuple{Int,DataType}, Vector{BigFloat}}(maxsize = 1024)

# ===== admissibility conditions ===#

"""
    δ(j1, j2, j3) -> Bool
    qδ(j1, j2, j3, k) -> Bool

Classical angular momentum triangle inequality: j3 must satisfy 
|j1-j2| ≤ j3 ≤ j1+j2 and j1+j2+j3 must be an integer.
Quantum triangle inequality: satisfies classical δ and sum j1+j2+j3 ≤ k.
"""
@inline ishalfInt(j)::Bool = isinteger(2*j) && j ≥ 0

@inline δ(j1, j2, j3)::Bool = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)

@inline qδ(j1, j2, j3, k)::Bool = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 

"""
    qδtet(j1, j2, j3, j4, j5, j6, k) -> Bool

Admissibility condition for a tetrahedral (6j symbol) configuration.
"""
@inline qδtet(j1, j2, j3, j4, j5, j6, k)::Bool = 
    ishalfInt(j1) && 
    ishalfInt(j2) && 
    ishalfInt(j3) && 
    ishalfInt(j4) && 
    ishalfInt(j5) && 
    ishalfInt(j6) && 
    qδ(j1, j2, j3, k) && 
    qδ(j1, j5, j6, k) && 
    qδ(j2, j4, j6, k) && 
    qδ(j3, j4, j5, k)

# ===== log q-integers and q-factorial tables ===#

# quantum integers 

"""
    qinteger(n, k) -> BigFloat

Quantum integer [n]_q for SU(2)_k at root of unity
q = exp(iπ/(k+2)).

Defined as:
    [n]_q = sin(nπ/(k+2)) / sin(π/(k+2))
"""
@inline function qinteger(n::Int, k::Int)::BigFloat
    θ = big(pi) / (k + 2)
    return sin(n * θ) / sin(θ)
end

"""
    qdim(j, k) -> BigFloat

Quantum dimension of SU(2)_k irrep with spin j:
dim_q(j) = [2j+1]_q
"""
@inline qdim(j, k::Int)::BigFloat = qinteger(Int(2j + 1), k)

"""
    logqn_table(k::Int, ::Type{T}=BigFloat) -> Vector{T}

Compute table of log q-integers: logqn[n+1] = log([n]_q) for n = 0..k+1,
where q = exp(iπ/(k+2)).

Uses trigonometric symmetry to compute only half the values.
"""
function logqn_table(k::Int, ::Type{T}=BigFloat)::Vector{BigFloat} where {T<:AbstractFloat}
    N = k + 2
    logqn = Vector{T}(undef, N)

    θ = T(pi) / N
    logden = log(sin(θ))

    # Use symmetry: n ↔ k+2-n
    half = N ÷ 2
    logqn[1] = zero(T)
    @inbounds for n in 1:half
        val = log(sin(n * θ)) - logden
        logqn[N + 1 - n] = logqn[n+1] = val
    end
    return logqn
end

"""
    logqnfact_table(k::Int, ::Type{T}=BigFloat) -> Vector{T}

Compute table of log quantum factorials: logqnfact[n+1] = log([n]_q!) for n = 0..k+1.
"""
function logqnfact_table(k::Int, ::Type{T}=BigFloat)::Vector{BigFloat} where {T<:AbstractFloat}
    return cumsum(logqn_table(k, T))
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
    _qracah6j(model::SU2kModel, j1, j2, j3, j4, j5, j6) -> BigFloat

Compute quantum 6j symbol using Racah summation formula.
Returns 0 if admissibility conditions are not met.
"""
function _qracah6j(model::SU2kModel{T}, j1, j2, j3, j4, j5, j6)::T where T
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return zero(T)
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

    zrange = max(α1, α2, α3, α4):min(β1, β2, β3)
    res = zero(BigFloat)
    for z in zrange
        logTsz = logT + log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, table) 
        res += iseven(z) ? exp(logTsz) : -exp(logTsz) 
    end
    return res
end

"""
    qracah6j(j1, j2, j3, j4, j5, j6, k::Int) -> BigFloat

Compute quantum 6j symbol {j1 j2 j3; j4 j5 j6} at level k.

Creates a temporary SU2kModel; use _qracah6j directly for better performance 
when computing multiple symbols with the same k.
"""
qracah6j(j1, j2, j3, j4, j5, j6, k::Int) =
    _qracah6j(SU2kModel(k), j1, j2, j3, j4, j5, j6)

qracah6j(model::SU2kModel, j1, j2, j3, j4, j5, j6) =
    _qracah6j(model, j1, j2, j3, j4, j5, j6)



# ============================================================
# R-symbols (braiding)
# ============================================================

"""
    Rsymbol(model, j1, j2, j3) -> Complex{BigFloat}

Braiding eigenvalue R^{j1 j2}_{j3}.
"""
function rsymbol(model::SU2kModel, j1,j2,j3)::Complex{BigFloat}
    if !qδ(j1,j2,j3, model.k)
        return zero(Complex{BigFloat})
    end

    θ = big(pi) / (model.k + 2)
    phase = exp(im * θ *
        (j1*(j1+1) + j2*(j2+1) - j3*(j3+1)))

    return phase
end


# ============================================================
# F-symbols and G-symbols
# ============================================================

"""
    Fsymbol(model, j1,j2,j3, j4,j5,j6) -> BigFloat

Fusion (F) symbol for SU(2)_k, defined by:

    F = (-1)^{Σj} √(dim_q(j3) dim_q(j6)) {6j}_q
"""
function _fsymbol(model::SU2kModel,
                 j1,j2,j3,
                 j4,j5,j6)

    q6j = _qracah6j(model, j1,j2,j3, j4,j5,j6)
    if iszero(q6j)
        return zero(BigFloat)
    end

    pref = sqrt(qdim(j3, model.k) * qdim(j6, model.k))
    phase = (-1)^(j1 + j2 + j4 + j5)

    return phase * pref * q6j
end

fsymbol(j1,j2,j3, j4,j5,j6, k::Int) =
    _fsymbol(SU2kModel(k),  j1,j2,j3, j4,j5,j6)

fsymbol(model::SU2kModel, j1,j2,j3, j4,j5,j6) =
    _fsymbol(model,  j1,j2,j3, j4,j5,j6)

"""
    gsymbol(model::SU2kModel, j1, j2, j12, j3, j123, j23) -> BigFloat

Compute G-symbol (recoupling coefficient) for SU(2)_k.
Relates two different couplings of three angular momenta:
|j1,j2;j12⟩|j3;j123⟩ = Σ_{j23} fsymbol(...) |j1;j13⟩|j2,j3;j23⟩

Uses quantum 6j symbol: F = {j1 j2 j12; j3 j123 j23}
"""
function _gsymbol(model::SU2kModel, j1, j2, j3, j4, j5, j6)
    sumj = j1 + j2 + j3 + j4 + j5 + j6
    !ishalfInt(sumj) && return zero(Complex{BigFloat})
    q6j = _qracah6j(model, j1, j2, j3, j4, j5, j6)
    return isinteger(sumj) ? ( (-1)^sumj * q6j ) : (-1+0im)^sumj * q6j 
end

gsymbol(j1, j2, j3, j4, j5, j6, k::Int) =
    _gsymbol(SU2kModel(k), j1, j2, j3, j4, j5, j6)


"""
    fsymbol(model::SU2kModel, j1, j2, j12, j3, j123, j23) -> BigFloat

Compute F-symbol (recoupling coefficient) for SU(2)_k.
Relates two different couplings of three angular momenta:
|j1,j2;j12⟩|j3;j123⟩ = Σ_{j23} fsymbol(...) |j1;j13⟩|j2,j3;j23⟩

Uses quantum 6j symbol: F = {j1 j2 j12; j3 j123 j23}
"""

# ============================================================
# q-3j symbols and Clebsch–Gordan coefficients
# ============================================================

"""
    q3j(model, j1,j2,j3, m1,m2,m3) -> BigFloat

Quantum 3j symbol for SU(2)_k.

Selection rules:
- m1 + m2 + m3 = 0
- |mi| ≤ ji
- (j1,j2,j3) admissible
"""
function q3j(model::SU2kModel, j1,j2,j3, m1,m2,m3)
    if m1 + m2 + m3 != 0
        return zero(BigFloat)
    end
    if !qδ(j1,j2,j3, model.k)
        return zero(BigFloat)
    end
    if abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3
        return zero(BigFloat)
    end

    # Express via q-CG normalization
    pref = (-1)^(Int(j1 - j2 - m3)) / sqrt(qdim(j3, model.k))
    return pref * qclebschgordan(model, j1,m1, j2,m2, j3,-m3)
end

"""
    qclebschgordan(model, j1,m1, j2,m2, j3,m3) -> BigFloat

Quantum Clebsch–Gordan coefficient ⟨j1 m1, j2 m2 | j3 m3⟩_q.
"""
function qclebschgordan(model::SU2kModel,
                        j1,m1,
                        j2,m2,
                        j3,m3)

    if m1 + m2 != m3
        return zero(BigFloat)
    end
    if !qδ(j1,j2,j3, model.k)
        return zero(BigFloat)
    end

    # Normalization factor
    norm = sqrt(qdim(j3, model.k))

    # Phase convention (Kirillov–Reshetikhin)
    phase = (-1)^(Int(j1 - j2 + m3))

    return phase * norm * q3j(model, j1,j2,j3, m1,m2,-m3)
end

end