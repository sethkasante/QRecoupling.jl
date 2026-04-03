module QRacahSymbols

# ---- dependencies ----
using LRUCache
using Nemo
using PrecompileTools


# --- global type aliases ----
# Define spin as Real to allow for both Integer and Half-Integer (Float) representations.
const Spin = Real
const OptInt = Union{Nothing, Int}


# -----  core architecture & data structures ------
include("types.jl")
include("admissibility.jl")
include("symmetries.jl")

# DCR builders 
include("cyclo_builder.jl")

# Computation Engines
include("direct_numeric_lse.jl")     # fast Log-Sum-Exp algo
include("projection_discrete.jl")   # level k
include("projection_analytic.jl")   #complex sweeps
include("projection_classical.jl")  # q→1 limit 
include("projection_exact.jl")      # project DCR to cyclotomic field using Nemo
include("eager_exact.jl")          # eager Nemo computation

# tqft categories
include("topological_symbols.jl")


# Public User Interface
# Contains the master dispatchers and caches
include("main_interfaces.jl")



# public apis export functions 

export q3j, 
       q6j, 
       qint,
       qdim, 
       rmatrix, 
       fsymbol, 
       gsymbol

export q6j_exact, q3j_exact

#constructors for results
export ExactResult, 
       CycloResult, 
       ClassicalResult, 
       CycloExactResult

#evaluation of cyclotomic representation
export evaluate_cyclo

#clear caches 
export clear_caches!, 
       clear_numeric_caches!, 
       clear_sieve_caches!, 
       clear_exact_caches!


# precompile these 
@setup_workload begin
    k = 3
    @compile_workload begin
        q6j(1, 1, 1, 1, 1, 1, k; mode=:numeric)
        q6j(1, 1, 1, 1, 1, 1, k; mode=:exact)
        q6j(1, 1, 1, 1, 1, 1; mode=:classical_exact)
        q6j(1, 1, 1, 1, 1, 1; mode=:cyclo)
    end
end

end # module QRacahSymbols