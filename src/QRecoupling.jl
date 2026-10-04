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

include("modular_arithmetic.jl")
include("factorial_rule.jl")
include("level_tables.jl")
include("double_word.jl")
include("factorial_kernels.jl")
include("multiword.jl")
include("classical_sums.jl")
include("classical_exact.jl")
include("exact_cleared.jl")   # exact level values from the rule, without field arithmetic
include("families.jl")
include("fmatrix.jl")       # the associativity matrix, from shared column coefficients
include("workspace.jl")
include("analytic_rules.jl")
include("factorial_series.jl")
include("generic_x.jl")     # exact values in x = q + 1/q: generic q, canonical forms, identity proofs
include("phi_form.jl")      # the closed form a user sees: Φ factors in q
include("exact_x.jl")       # exact level values in x = 2cos(π/h), and the radical display ladder
include("exact_radicals.jl") # radical expressions and the Lagrange descent
include("exact_display.jl")  # rendering of exact values
include("exact_arithmetic.jl") # arithmetic, sums by square class, exact equality
include("symbolic_value.jl") # what Symbolic() returns: the rule and its DCR, shown in x

include("modular.jl")      # modular data: twists, S and T, the total dimension
include("qproducts.jl")   # [n], [n]!, [n choose m] through the same rule evaluator
include("recoupling_api.jl")
include("batch.jl")
include("symbols.jl")   # the symbol interface: one table, dispatched
include("queries.jl")
include("targets.jl")
include("qcg.jl")       # quantum Clebsch–Gordan coefficients and the quantum 3j symbol: weighted rules
include("qcg_columns.jl") # coupling columns by the Casimir recurrence: qcg_matrix and the near-edge tier


# Export physics tqft and recoupling symbols api
export q6j, q3j, qcg, qcg_matrix, qcg_row, fsymbol, gsymbol, rmatrix, qdim
export fmatrix
export smatrix, tmatrix, bmatrix, twist
export qint, qfact, qbinomial, qseries, qeval

# evaluation targets
export Level, Exact, At, Classical, Symbolic

# exact values in x = q + 1/q: generic q, and identities proved for every level at once
export prove_identity, x_form, XValue, phi_form
export ExactX, ExactXSum, RadExpr, NoRadical, radical, numeric_value

# structural queries and batches
export all_6j, iszero_at, issingular_at, level_spectrum

# factorial rules, workspaces and returned types
export FactorialSum, AffineFactorial, EvaluationWorkspace, SymbolicValue, QPhase
export empty_caches!

end
