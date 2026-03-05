# src/QuantumRacah.jl (find a better name?)
module QuantumRacah

using Nemo
using LRUCache

# Export Models & Types
export CycloMonomial, ExactSU2kModel, NumericSU2kModel
export Generic6jResult, Exact6jResult

# Export Primary User API
export q6j, qracah6j_generic, qracah6j_exact, qracah6j_numeric, qracah6j,
        qracah3j_generic

# Include architecture files
include("Types2.jl")
include("Symbolics2.jl")
include("ExactAlgebra2.jl")
include("Numerics2.jl")
include("admissible.jl")




# Initialize the global Evaluation Caches
# Keys: (Canonical_Spin_Tuple, k_level)
const Q6J_NUMERIC_CACHE = LRU{Tuple{NTuple{6, Float64}, Int}, BigFloat}(maxsize=50000)
const Q6J_EXACT_CACHE   = LRU{Tuple{NTuple{6, Float64}, Int}, Exact6jResult}(maxsize=10000)

# ============================================================
# The Unified Researcher API
# ============================================================

"""
    q6j(j1, j2, j3, j4, j5, j6, k::Union{Int, Nothing}=nothing; mode=:numeric, prec=256)

The universal quantum 6j-symbol evaluator.

# Modes:
- `:numeric` (Default): Returns a fast BigFloat. Requires `k`.
- `:exact`: Returns an Exact6jResult in the Nemo cyclotomic field. Requires `k`.
- `:generic`: Returns a Generic6jResult (CycloMonomials). `k` is ignored.
"""
function q6j(j1, j2, j3, j4, j5, j6, k::Union{Int, Nothing}=nothing; mode=:numeric, prec=256)
    
    # 1. Generic / Pure Symbolic Mode (No k needed)
    if mode == :generic || k === nothing
        pref_sq = qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
        series = q6jseries_symb(j1, j2, j3, j4, j5, j6)
        return Generic6jResult(pref_sq, series)
    end
    
    # Admissibility Check for level k
    if !qδtet(j1, j2, j3, j4, j5, j6, k) 
        return mode == :exact ? nothing : zero(BigFloat)
    end

    # 2. Exact Algebraic Mode (Nemo)
    if mode == :exact
        model = ExactSU2kModel(k)
        pref2_nf = qtricoeff2_exact(j1, j2, j3, j4, j5, j6, model)
        sum_nf = q6jseries_exact(j1, j2, j3, j4, j5, j6, model)
        return Exact6jResult(k, pref2_nf, sum_nf)
    end
    
    # 3. Fast Numeric Mode (BigFloat)
    if mode == :numeric
        model = SU2kModel(k, prec)
        return _qracah6j_stable(model, j1, j2, j3, j4, j5, j6) # From Numerics.jl
    end
    
    error("Unknown mode: $mode. Use :numeric, :exact, or :generic.")
end



#TODO: Add q3j, rmatrix, fsymbols, gsymbols and other tqft functions 

end # module