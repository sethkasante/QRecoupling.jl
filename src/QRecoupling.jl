module QRecoupling

using LRUCache
using PrecompileTools: @setup_workload, @compile_workload

"Spin labels: multiples of 1/2 given as integers, rationals or floats."
const Spin = Union{Integer, Rational, AbstractFloat}

# --- Internal Includes ---
include("helpers.jl")
include("classical_result.jl")
include("admissibility.jl")
include("symmetries.jl")
include("qphases.jl")

include("modular_arithmetic.jl")
include("factorial_rule.jl")
include("level_tables.jl")
include("double_word.jl")
include("factorial_kernels.jl")
include("multiword.jl")
include("classical_sums.jl")
include("classical_exact.jl")
include("families.jl")
include("fmatrix.jl")       # the associativity matrix, from shared column coefficients
include("workspace.jl")
include("analytic_rules.jl")
include("factorial_series.jl")
include("algebra_interface.jl") # stable wrappers and optional-backend entry points

include("modular.jl")      # modular data: twists, S and T, the total dimension
include("qproducts.jl")   # [n], [n]!, [n choose m] through the same rule evaluator
include("recoupling_api.jl")
include("batch.jl")
include("symbols.jl")   # the symbol interface: one table, dispatched
include("queries.jl")
include("targets.jl")
include("qcg.jl")       # quantum Clebsch–Gordan coefficients and the quantum 3j symbol: weighted rules
include("qcg_columns.jl") # coupling columns by the Casimir recurrence: qcg_matrix and the near-edge tier
include("qform.jl")     # values written in quantum integers, generic and at a level
include("level_bounds.jl") # outward enclosures used only by exact level certification
include("level_field.jl") # exact level values in quantum integers, in plain Julia
include("level_arithmetic.jl") # products, inverses, sums and exact zero tests of those values
include("qform_value.jl") # a form as a value: equality and arithmetic


# Export physics tqft and recoupling symbols api
export q6j, q3j, qcg, qcg_matrix, qcg_row, fsymbol, gsymbol, rmatrix, qdim
export fmatrix, fmatrix_labels
export smatrix, tmatrix, bmatrix, twist
export qint, qfact, qbinomial, qseries, qeval

# evaluation targets
export Level, Exact, At, Classical, Symbolic

# exact values in x = q + 1/q: generic q, and identities proved for every level at once
export prove_identity, xform, XValue, LevelXValue
export RadExpr, NoRadical, radical, numeric_value
export QForm

# structural queries and batches
export all_6j, iszero_at, issingular_at, level_spectrum

# factorial rules, workspaces and returned types
export FactorialSum, AffineFactorial, EvaluationWorkspace, QPhase
export empty_caches!

# Run each target once while precompiling, since a first call otherwise compiles for seconds. The tables it builds are
# cleared again and do not reach the image.
@setup_workload begin
    labels = (1, 1, 1, 1, 1, 1)
    @compile_workload begin
        q6j(labels...); q6j(Level(10), labels...); q6j(At(0.8), labels...); q6j(At(0.8 + 0.1im), labels...)
        q6j(Level(10), [labels, (2, 2, 2, 2, 2, 2)]); q6j(labels...; exact = true)
        fsymbol(Level(10), labels...); gsymbol(Level(10), labels...); fmatrix(1, 1, 1, 1; k = 10)
        rmatrix(Level(10), 1, 1, 1); qdim(Level(10), 1); twist(Level(10), 1); smatrix(10); tmatrix(10)
        qcg(Level(10), 1, 0, 1, 0, 2, 0); qcg(1, 0, 1, 0, 2, 0); q3j(Level(10), 1, 1, 2, 0, 0, 0); q3j(1, 1, 2, 0, 0, 0)
        qint(Level(10), 3); qfact(Level(10), 3); qbinomial(Level(10), 4, 2)
        q6j(Exact(), labels...); qcg(Exact(), 1, 0, 1, 0, 2, 0); qint(Exact(), 3)
        for f in (q6j(Symbolic(), labels...), q6j(Exact(10), labels...), qcg(Symbolic(), 1, 0, 1, 0, 2, 0), qcg(Exact(10), 1, 0, 1, 0, 2, 0), qint(Symbolic(), 3))
            sprint(show, MIME("text/plain"), f); iszero(f); f == f; hash(f); f == 1
            f.k === nothing ? (f(1); f(0.8)) : float(f)
        end
        a = q6j(Exact(10), labels...); b = q6j(Exact(10), 2, 2, 2, 2, 2, 2)
        a == b; a == 1//2; a * b; a + b; inv(a)
        iszero_at(q6j, 10, labels...); level_spectrum(labels...; k = 2:6)
    end
    empty_caches!()
end

end
