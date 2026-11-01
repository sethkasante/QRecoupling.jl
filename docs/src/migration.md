# Migration

## To v0.6

### Exact values are forms

`Exact()` and `Exact(k)` return a [`QForm`](@ref): the value written in quantum integers, from its formula.

| Before | Now |
|---|---|
| `q6j(Exact(), js...)` was a rational or radical | A form at every q. `== 1//6`, `Float64`, `numeric_value` and arithmetic still read its value at q = 1, and `classical(v)` gives that value exactly |
| `q6j(Exact(k), js...)` was an `ExactX`, and needed Nemo | A form at the level, with no backend. `==`, `iszero` and arithmetic are exact |
| `q3j` and `qcg` had no `Exact(k)` | `Exact()` and `Exact(k)` write them with their powers of q |
| `v.x_value`, `v.rad` | `QRecoupling.xpolynomial(v)`, `QRecoupling.radicand(v)`, `radical(v)`, with Nemo |
| `ExactXSum(k)` to accumulate a sum | `sum`, or start from `zero(v)` |
| `isempty(residual)` | `iszero(residual)`, or `==` |

`q6j(js...; exact = true)` and `Classical(exact = true)` still return the classical number, and cost less than a form when that is all you need. `q6j(js...; k = k, exact = true)` is the same as `q6j(Exact(k), js...)`.

Between two `Exact()` values, `==` now means equal at every q: `qint(Exact(), 6) == qfact(Exact(), 3)` is false, although both are 6 at q = 1. Compare `classical(a) == classical(b)` for the old meaning.

### Optional backend

Nemo is no longer a required dependency. Install it with `import Pkg; Pkg.add("Nemo")` and load it with `using QRecoupling, Nemo` when using `radical`, `x_form`, `phi_form` or `prove_identity`. Both package load orders are supported. Missing-backend errors explain how to enable it.

`ExactX`, `ExactXSum`, and `XValue` now parameterize their polynomial payload. `QRecoupling.exact_x(rule, k)` still gives an `ExactX`, with its arithmetic, when Nemo is loaded.

### Removed after their deprecation

| Removed | Use instead |
|---|---|
| `Exact(k; form = :canonical)` and the `CompositeExactResult` values it returned | `Exact(k)` |
| `form = :x` | `Exact(k)`, which is the same thing |
| `eager = true` | Remove the keyword |
| `QRecoupling.xvalue(value)` with one argument | `x_form(value)` |
| `qeval(m; k, exact = true)` for a monomial or DCR | `qeval(Exact(k), rule)` with a `FactorialSum` |
| `QRecoupling.project_exact`, `QRecoupling.evaluate_exact` on the removed type | `Exact(k)` and `numeric_value` |

The three-argument `QRecoupling.xvalue(rad, num, den)`, numerical evaluation of monomials and DCRs, and exact classical evaluation of them are unchanged.

## To v0.5

### `q3j` is the quantum 3j symbol

`q3j` is now the 3j symbol of U_q(sl₂): the q-Clebsch–Gordan coefficient [`qcg`](@ref) with the 3j phase and normalisation. It includes the deformation-dependent weights missing in v0.4. Classical values are unchanged; deformed values change, and level values are complex.

`q3j` and `qcg` support numerical evaluation and exact classical values through `Exact()`. `Exact(k)` and `Symbolic()` were not implemented for these weighted coefficients in v0.5.

For existing code that requires the former formula, `QRecoupling.q3j_factorial` remains available with its previous targets.

### Negative real q is complex

At negative real q every symbol uses the square-root branch of `complex(q)`, so results are complex and agree with the same point supplied as complex. Positive real q, classical values and levels are unchanged.

```@example migration
using QRecoupling, Nemo
q6j(1, 1, 1, 1, 1, 1; q = -0.8) == q6j(1, 1, 1, 1, 1, 1; q = complex(-0.8))
```

### `T` is a precision floor in `fmatrix`

As in the scalar functions, `T` in `fmatrix` no longer narrows a complex or higher-precision result.

### A smaller export list

These names are no longer exported. They still exist: call them as `QRecoupling.name`, or bring them in with `using QRecoupling: name`.

| Group | Names |
|---|---|
| Cyclotomic/DCR layer (compatibility; planned for retirement) | `CyclotomicMonomial`, `DCR`, `build_dcr!`, `add_qint!`, `add_qfact!`, `build_series`, `qint_mono`, `qfact_mono`, `qbinomial_mono`, `symbolic_terms`, `project_discrete`, `project_exact`, `project_analytic`, `project_classical`, `project_classical_exact` |
| Symbol interface, for adding symbols | `QSymbol`, `SixJ`, `ThreeJ`, `FSymbol`, `GSymbol`, `Tetrahedron`, `ThetaValue`, `symbol_rule`, `level_admissible`, `symbol_family`, `symbol_of`, `nlabels` |
| Lower-level exact access | `exact_x` (use `Exact(k)`), `generic_sixj` (use `x_form(q6j(Symbolic(), …))`), `xvalue` (the `XValue` constructor; use `x_form` to expand), `radical_form`, `has_radical_form`, `radical_levels`, `xpolynomial`, `radicand` (use `v.x_value`), `splits_completely` |
| Modular data helpers | `level_labels`, `central_charge`, `total_qdim`, `gauss_sum`, `monodromy`, `verlinde` |
| Other | `EvalTarget` |

`smatrix`, `tmatrix`, `bmatrix`, `twist`, `rmatrix` and `fmatrix` remain exported.

### DCR compatibility

The DCR constructors, monomial helpers, callback form of `qseries`, and raw projection functions remain available for existing code and are planned for retirement. New code should use `FactorialSum` or factorial-triple `qseries`, `qeval`, and the evaluation targets. `value.dcr` still constructs a compatibility representation lazily; callbacks must return cyclotomic monomials. The public `phi_form(value)` operation remains available for cyclotomic factorization.

### New

[`qcg`](@ref), [`qcg_matrix`](@ref) and [`qcg_row`](@ref) give coupling coefficients, whole coupling matrices and the coupled decomposition of a product state. `qcg` accepts its labels flat, `qcg(j1, m1, j2, m2, j, m)`, or as pairs, `qcg((j1, m1), (j2, m2), (j, m))`.

## To v0.4

### Numerical defaults now include the product helpers

Symbol functions, `qeval`, `qint`, `qfact`, and `qbinomial` evaluate classically when no q or k is supplied. Code that used their return values as symbolic objects must request `Symbolic()` explicitly.

```@example migration
using QRecoupling, Nemo
q6j(1,1,1,1,1,1)                      # approximately 1/6
q6j(Exact(),1,1,1,1,1,1)          # exact classical coefficient
s = q6j(Symbolic(),1,1,1,1,1,1)       # rule-backed symbolic x-form
@assert qint(5) == 5 # hide
qint(Symbolic(),5)
```

| Earlier pattern | v0.4 replacement |
|:--|:--|
| `q6j(js...)` for a symbolic expression | `q6j(Symbolic(), js...)` |
| `qfact(n)` as a symbolic expression | `qfact(Symbolic(), n)` |
| Monomial callbacks for factorial sums | `FactorialSum` or factorial-triple `qseries` |
| Exact level cyclotomic wrapper by default | `Exact(k)` returns `ExactX` |
| `eager=true` | Remove the keyword |

`SymbolicValue` retains a factorial rule; it does not have the old monomial multiplication interface. Build product rules using factorial prefactors or expand with `x_form` when algebra is needed.

### Exact classical evaluation

In v0.4 and v0.5, `Exact()` was an alias for `Classical(exact=true)`, and it also worked with the other recoupling symbols and product helpers. From v0.6 it returns a form; see above.

`exact=true` and `Classical(exact=true)` remain supported. `Exact()` takes no level or form option: use `Exact(k)` for a level.

### Exact levels use x

In v0.4 and v0.5, `Exact(k)` and `k=10, exact=true` returned an `ExactX`, a polynomial in x. From v0.6 they return a form; see above. Read coefficients with `QRecoupling.xpolynomial` and `QRecoupling.radicand`, request nested square roots with `radical(v)`, and convert numerically with `Float64`, `BigFloat`, or `numeric_value`. `radical` returns `RadExpr` or a reasoned `NoRadical` result.

`Exact(k; form=:canonical)` retained the former cyclotomic representation until v0.6. `form=:deferred`, `ExactLevelWorkspace`, `canonicalize_exact`, and dedicated Biedenharn–Elliott verifier functions from development drafts are not current APIs. Write identities from the public exact arithmetic as in the [identity tutorial](tutorials/identities.md). Raw DCR/monomial projection remains a separate compatibility path; it is not automatically converted into an x-form by `qeval(Symbolic(), object)`.

### Symbolic output is bounded and lazy

Use `x_form(value)` to expand a symbolic rule. The earlier `QRecoupling.xvalue(value)` spelling was removed in v0.6; `value.x_value` on an exact level result instead exposes its stored polynomials.

Symbolic expressions display their finite factorial rule in x without carrying out the sum, for both small and large labels. Radicals use ψ factors. `x_form(value)` explicitly expands and caches the generic formula, returning an owned copy.

Dimensions and theta graphs use this interface too. R matrices and twists retain `QPhase` objects with an explicit branch over x. A classical exact R value is an integer sign; a symbolic or exact-level R value remains a phase.

### Targets and scratch storage

`Level(k)`, `Exact()`, `Exact(k)`, `At(q)`, and `Classical()` select evaluation. Conflicting targets and keywords are rejected. `exact=true` requires a classical or level target; use `Symbolic()` for generic formulas.

`EvaluationWorkspace()` can be reused sequentially by supported scalar calls. Batches manage worker workspaces automatically. `eager=true` is deprecated and delegates to the standard rule evaluator; it does not enable a different fast path.
