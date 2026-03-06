# src/QRacahSymbols.jl
module QRacahSymbols

using Nemo
using LRUCache

# Export Data Structures
export CycloMonomial, ExactSU2kModel, NumericSU2kModel
export GenericResult, ExactResult
export Spin

# Export Primary API Functions
export q6j, q3j, fsymbol, rmatrix, gsymbol, qdim
export evaluate_exact, evaluate_classical
export qracah6j_generic, qracah6j_exact, qracah6j_numeric, qracah3j_generic

# Include architecture files
include("Types.jl")
include("admissible.jl")
include("symmetry.jl")
include("Symbolics.jl")
include("ExactAlgebra.jl")
include("Numerics.jl")
include("tqft.jl")

# Initialize the global Evaluation & Model Caches
const Q6J_NUMERIC_CACHE = LRU{Tuple{NTuple{6, Float64}, Int, DataType}, Any}(maxsize=50000)
const Q6J_EXACT_CACHE   = LRU{Tuple{NTuple{6, Float64}, Int}, ExactResult}(maxsize=10000)

# Model Caches (to prevent rebuilding fields/tables on repeated calls)
const EXACT_MODEL_CACHE = LRU{Int, ExactSU2kModel}(maxsize=200)

# ============================================================
# The Unified Researcher API
# ============================================================

"""
    q6j(j1, j2, j3, j4, j5, j6, k::Union{Int, Nothing}=nothing; mode=:numeric, T=Float64, prec=256)

The universal quantum 6j-symbol evaluator.

# Modes:
- `:numeric` (Default): Fast log-sum-exp via `T` (Float64 or BigFloat). Requires `k`.
- `:exact`: ExactResult in the Nemo cyclotomic field. Requires `k`.
- `:generic`: GenericResult (CycloMonomials). `k` is ignored.
- `:classical`: Ponzano-Regge Float64 limit (q -> 1). `k` is ignored.
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Union{Int, Nothing}=nothing; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    
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
        return mode == :exact ? ExactResult(k, 0, 0) : zero(T)
    end

    # 3. Canonicalize spins for S4 Tetrahedral Cache Lookups
    c_spins = canonical_spins(j1, j2, j3, j4, j5, j6)
    
    # 4. Exact Algebraic Mode (Nemo)
    if mode == :exact
        cache_key_exact = (c_spins, k)
        return get!(Q6J_EXACT_CACHE, cache_key_exact) do
            # Pull the Exact Model from cache so the Nemo field is only generated ONCE per k
            model = get!(EXACT_MODEL_CACHE, k) do
                ExactSU2kModel(k)
            end
            qracah6j_exact(model, c_spins...)
        end
    end
    
    # 5. Fast Numeric Mode (Float64 / BigFloat)
    if mode == :numeric
        cache_key_numeric = (c_spins, k, T)
        return get!(Q6J_NUMERIC_CACHE, cache_key_numeric) do
            # The Numeric model has its own internal global LRU cache for tables
            model = NumericSU2kModel(k; T=T, prec=prec)
            _qracah6j_stable(model, c_spins...)
        end
    end
    
    error("Unknown mode: $mode. Use :numeric, :exact, :generic, or :classical.")
end



# ============================================================
# Quantum 3j API
# ============================================================

"""
    q3j(j1, j2, j3, m1, m2, m3=nothing, k=nothing; mode=:numeric, T=Float64, prec=256)

The universal quantum 3j-symbol evaluator. 
If `m3` is omitted, it is enforced as `-m1 - m2`.
"""
function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Union{Spin, Nothing}=nothing, k::Union{Int, Nothing}=nothing; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    
    # Auto-assign m3 and check conservation
    m3_val = m3 === nothing ? -m1 - m2 : m3
    if !iszero(m1 + m2 + m3_val) || !δ(j1, j2, j3)
        # Handle mathematically vanishing states
        if mode == :generic return GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[]) end
        if mode == :exact return ExactResult(k === nothing ? 0 : k, 0, 0) end
        return zero(T)
    end

    if mode == :generic || k === nothing
        return qracah3j_generic(j1, j2, j3, m1, m2, m3_val)
    end
    
    # Check bounds against level k
    if !qδ(j1, j2, j3, k) 
        return mode == :exact ? ExactResult(k, 0, 0) : zero(T)
    end

    if mode == :exact
        model = get!(EXACT_MODEL_CACHE, k) do
            ExactSU2kModel(k)
        end
        return qracah3j_exact(model, j1, j2, j3, m1, m2)
    end
    
    if mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return _qracah3j_stable(model, j1, j2, j3, m1, m2)
    end
    
    error("Unknown mode: $mode")
end

end # module