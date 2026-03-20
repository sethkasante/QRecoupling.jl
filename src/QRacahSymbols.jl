module QRacahSymbols

# ---- Dependencies ----
using LRUCache
using Nemo


# --- Global Type Aliases ----
# Define spin as Real to allow for both Integer and Half-Integer (Float) representations.
const Spin = Real
const OptInt = Union{Nothing, Int}


# -----  Core Architecture & Data Structures ------
include("types.jl")
include("admissibility.jl")
include("symmetries.jl")

# Abstract Math Builders
include("cyclo_builder.jl")

# Computation Engines
include("evaluator_direct.jl")     # Fast Log-Sum-Exp Floats
include("evaluator_discrete.jl")   # Exact Roots of Unity Sieve
include("evaluator_analytic.jl")   # Continuous Complex sweeps
include("evaluator_classical.jl")  # q=1 Limits & Zero-Allocation GMP
include("evaluator_exact.jl")      # Rigorous Nemo.jl Cyclotomic Fields

# tqft symbols 
include("topological_symbols.jl")

# Public User Interface
# Contains the master dispatchers and cache-clearing API
include("main_interfaces.jl")



# public apis export functions 

export q3j, 
       q6j, 
       qdim, 
       rmatrix, 
       fsymbol, 
       gsymbol

export cyclo_to_numeric, 
       evaluate_level_exact

#clear caches 
export clear_caches!, 
       clear_numeric_caches!, 
       clear_sieve_caches!, 
       clear_exact_caches!

# Builders and classical specific exports (for internal use)
export q6j_cyclo, q3j_cyclo
export q6j_classical, q3j_classical, q6j_classical_exact, q3j_classical_exact
export q6j_direct, q3j_direct


#constructors for results
export ExactResult, CycloResult, ClassicalResult

end # module QRacahSymbols