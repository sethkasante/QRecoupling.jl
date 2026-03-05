# src/QRacahSymbols.jl
module QRacahSymbols

using Nemo
using LRUCache

# Export Data Structures
export CycloMonomial, ExactSU2kModel, NumericSU2kModel
export GenericResult, ExactResult

# Export Primary API Functions
# export q6j, q3j, fsymbol, rmatrix, gsymbol, qdim
# export evaluate_exact, evaluate_classical

export q6j, qracah6j_generic, qracah6j_exact, qracah6j_numeric, qracah6j,
        qracah3j_generic, evaluate_cyclofield, qtricoeff2_exact

# Include architecture files
include("admissible.jl")
include("Symmetry.jl")
include("Types.jl")
include("Symbolics.jl")
include("ExactAlgebra.jl")
include("Numerics.jl")
include("TQFT.jl")

# Initialize the global Evaluation Caches
# Keys: (Canonical_Spin_Tuple, k_level)
const Q6J_NUMERIC_CACHE = LRU{Tuple{NTuple{6, Float64}, Int}, BigFloat}(maxsize=50000)
const Q6J_EXACT_CACHE   = LRU{Tuple{NTuple{6, Float64}, Int}, ExactResult}(maxsize=10000)

# ============================================================
# The Unified Researcher API
# ============================================================

"""
    q6j(j1, j2, j3, j4, j5, j6, k::Union{Int, Nothing}=nothing; mode=:numeric, prec=256)

The universal quantum 6j-symbol evaluator.

# Modes:
- `:numeric` (Default): Fast BigFloat log-sum-exp. Requires `k`.
- `:exact`: ExactResult in the Nemo cyclotomic field. Requires `k`.
- `:generic`: GenericResult (CycloMonomials). `k` is ignored.
- `:classical`: Ponzano-Regge Float64 limit (q -> 1).
"""
function q6j(j1, j2, j3, j4, j5, j6, k::Union{Int, Nothing}=nothing; mode=:numeric, prec=256)
    
    # 1. Pure Symbolic Mode (No k needed)
    if mode == :generic || k === nothing
        return qracah6j_generic(j1, j2, j3, j4, j5, j6)
    end
    
    # 2. Classical Limit (q=1)
    if mode == :classical
        return qracah6j_classical(j1, j2, j3, j4, j5, j6)
    end

    # Admissibility Check for level k
    if !qδtet(j1, j2, j3, j4, j5, j6, k) 
        return mode == :exact ? ExactResult(k, 0, 0) : zero(BigFloat)
    end

    # 3. Canonicalize spins for Cache Lookups
    c_spins = canonical_spins(j1, j2, j3, j4, j5, j6)
    cache_key = (c_spins, k)

    # 4. Exact Algebraic Mode (Nemo)
    if mode == :exact
        return get!(Q6J_EXACT_CACHE, cache_key) do
            model = ExactSU2kModel(k)
            qracah6j_exact(model, c_spins...)
        end
    end
    
    # 5. Fast Numeric Mode (BigFloat)
    if mode == :numeric
        return get!(Q6J_NUMERIC_CACHE, cache_key) do
            model = NumericSU2kModel(k, prec)
            _qracah6j_stable(model, c_spins...)
        end
    end
    
    error("Unknown mode: $mode. Use :numeric, :exact, :generic, or :classical.")
end

end # module