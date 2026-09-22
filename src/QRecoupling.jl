module QRecoupling

using Nemo, LRUCache

"Spin labels: multiples of 1/2 given as integers, rationals or floats."
const Spin = Union{Integer, Rational, AbstractFloat}

# --- Internal Includes ---
include("types_cyclotomic.jl")
include("cyclotomic_values.jl")
include("types_projection.jl")
include("admissibility.jl")
include("symmetries.jl")

include("generic_series_dcr.jl")
include("recoupling_symbols_dcr.jl") # recoupling symbols and tqft invariants
include("qphases.jl")

include("projection_classical.jl")
include("projection_discrete.jl")
include("projection_exact.jl")
include("projection_analytic.jl")
include("eager_discrete.jl")
include("eager_exact.jl")

include("modular_arithmetic.jl")
include("factorial_rule.jl")
include("level_tables.jl")
include("double_word.jl")
include("factorial_kernels.jl")
include("multiword.jl")
include("classical_sums.jl")
include("classical_exact.jl")
include("families.jl")
include("workspace.jl")
include("factorial_series.jl")

include("recoupling_api.jl")
include("batch.jl")
include("queries.jl")
include("targets.jl")


# Export physics tqft and recoupling symbols api
export q6j, q3j, fsymbol, gsymbol, rmatrix, qdim
export qint, qfact, qbinomial, qseries, qeval

#projection
export project_discrete, project_exact, project_analytic
export project_classical, project_classical_exact

# evaluation targets, batching and structural queries
export EvalTarget, Level, Exact, At, Classical, Symbolic
export all_6j, iszero_at, issingular_at, level_spectrum

# cache management
export empty_caches!

# Export api for generic series
export CyclotomicMonomial, DCR, QPhase
export AffineFactorial, FactorialSum, EvaluationWorkspace
export add_qint!, add_qfact!, build_dcr!, build_series

end
