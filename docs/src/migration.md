# Migration

## To v0.5

### `q3j` is the quantum 3j symbol

`q3j` is now the 3j symbol of U_q(sl₂): the q-Clebsch–Gordan coefficient [`qcg`](@ref) with the 3j phase and normalisation. Up to v0.4 it substituted quantum factorials into the classical formula without the q-power weights, which is not the U_q(sl₂) coefficient and is not orthogonal at q ≠ 1. Classical values are unchanged; values at real or complex q and at levels change, and level values are complex. The old function remains available, unexported, as `QRecoupling.q3j_factorial`, with all its targets (including `Exact(k)` and `Symbolic()`).

```@example migration
using QRecoupling
q3j(1, 1, 1, 1, -1, 0) == QRecoupling.q3j_factorial(1, 1, 1, 1, -1, 0)   # classical: identical
(q3j(1, 1, 1, 1, -1, 0; q = 0.8), QRecoupling.q3j_factorial(1, 1, 1, 1, -1, 0; q = 0.8))
```

### Negative real q is complex

At negative real q every symbol uses the square-root branch of `complex(q)`, so results are complex and agree with the same point supplied as complex. Positive real q, classical values and levels are unchanged.

```@example migration
q6j(1, 1, 1, 1, 1, 1; q = -0.8) == q6j(1, 1, 1, 1, 1, 1; q = complex(-0.8))
```

### `T` is a precision floor in `fmatrix`

As in the scalar functions, `T` in `fmatrix` no longer narrows a complex or higher-precision result.

### A smaller export list

These names are no longer exported. They still exist: call them as `QRecoupling.name`, or bring them in with `using QRecoupling: name`.

| Group | Names |
|---|---|
| Cyclotomic/DCR layer (deprecated, removed in v0.6) | `CyclotomicMonomial`, `DCR`, `build_dcr!`, `add_qint!`, `add_qfact!`, `build_series`, `qint_mono`, `qfact_mono`, `qbinomial_mono`, `symbolic_terms`, `project_discrete`, `project_exact`, `project_analytic`, `project_classical`, `project_classical_exact` |
| Symbol interface, for adding symbols | `QSymbol`, `SixJ`, `ThreeJ`, `FSymbol`, `GSymbol`, `Tetrahedron`, `ThetaValue`, `symbol_rule`, `level_admissible`, `symbol_family`, `symbol_of`, `nlabels` |
| Lower-level exact access | `exact_x` (use `Exact(k)`), `generic_sixj` (use `x_form(q6j(Symbolic(), …))`), `xvalue` (deprecated; use `x_form`), `radical_form`, `has_radical_form`, `radical_levels`, `xpolynomial`, `radicand` (use `v.x_value`), `splits_completely` |
| Modular data helpers | `level_labels`, `central_charge`, `total_qdim`, `gauss_sum`, `monodromy`, `verlinde` |
| Other | `EvalTarget`, `q3j_factorial` |

`smatrix`, `tmatrix`, `bmatrix`, `twist`, `rmatrix` and `fmatrix` remain exported.

### New

[`qcg`](@ref), [`qcg_matrix`](@ref) and [`qcg_row`](@ref) give coupling coefficients, whole coupling matrices and the coupled decomposition of a product state. `qcg` accepts its labels flat, `qcg(j1, m1, j2, m2, j, m)`, or as pairs, `qcg((j1, m1), (j2, m2), (j, m))`.

## To v0.4

### Numerical defaults now include the product helpers

Symbol functions, `qeval`, `qint`, `qfact`, and `qbinomial` evaluate classically when no q or k is supplied. Code that used their return values as symbolic objects must request `Symbolic()` explicitly.

```@example migration
using QRecoupling
q6j(1,1,1,1,1,1)                      # approximately 1/6
q6j(Exact(),1,1,1,1,1,1)          # exact classical coefficient
s = q6j(Symbolic(),1,1,1,1,1,1)       # rule-backed symbolic x-form
@assert qint(5) == 5 # hide
qint(Symbolic(),5)
```

| Earlier pattern | v0.4 replacement |
|:--|:--|
| `q6j(js...)` to construct a DCR | `q6j(Symbolic(), js...).dcr` |
| `qfact(n)` as a symbolic expression | `qfact(Symbolic(), n)` |
| `qfact(n)` inside a monomial callback | Prefer a factorial rule; otherwise `QRecoupling.qfact_mono(n)` |
| Exact level cyclotomic wrapper by default | `Exact(k)` returns `ExactX` |
| `eager=true` | Remove the keyword |

`qseries` remains a constructor: triples produce a `FactorialSum`, callbacks produce a DCR. Callback functions must return cyclotomic monomials; numeric product defaults cannot be substituted there. `SymbolicValue` is a rule-backed view, not a replacement monomial with the same multiplication interface. Build product rules using factorial prefactors or expand with `x_form` when algebra is needed.

### Exact classical evaluation

Prefer `q6j(Exact(), js...)` over `q6j(js...; exact=true)`, and `qeval(Exact(), rule)` for a classical factorial sum. `Exact()` is an alias for `Classical(exact=true)` and preserves the existing exact result types and fast evaluator. It also works with the other recoupling symbols and product helpers.

`exact=true` and `Classical(exact=true)` remain supported without deprecation warnings in v0.4. A future release may deprecate the keyword after a migration period; this release does not remove it. `Exact()` takes no level or form option: use `Exact(k)` for level arithmetic, or `Symbolic()` for a generic formula.

### Exact levels use x

```@example migration
v = q6j(Exact(10),1,1,1,1,1,1)
@assert v isa ExactX # hide
v
```

This also applies to `k=10, exact=true` for the recoupling and product functions. Read coefficients with `v.x_value`, `QRecoupling.xpolynomial`, and `QRecoupling.radicand`, request nested square roots with `radical(v)` (or `v.rad`), and convert numerically with `Float64`, `BigFloat`, or `numeric_value`. Printing shows the stored x-form; radicals and numerical approximations are explicit requests. `radical` returns `RadExpr` or a reasoned `NoRadical` result.

`Exact(k; form=:canonical)` retains the deprecated cyclotomic representation for compatibility. `form=:deferred`, `ExactLevelWorkspace`, `canonicalize_exact`, and dedicated Biedenharn–Elliott verifier functions from development drafts are not current APIs. Write identities from the public exact arithmetic as in the [identity tutorial](tutorials/identities.md). Raw DCR/monomial projection remains a separate compatibility path; it is not automatically converted into an x-form by `qeval(Symbolic(), object)`.

### Symbolic output is bounded and lazy

Use `x_form(value)` to expand a symbolic rule. The earlier `QRecoupling.xvalue(value)` spelling is deprecated; `value.x_value` on an exact level result instead exposes its stored polynomials.

Symbolic expressions display their finite factorial rule in x without carrying out the sum, for both small and large labels. Radicals use ψ factors. `x_form(value)` explicitly expands and caches the generic formula, returning an owned copy. `.dcr` constructs a compatibility DCR only when accessed; printing that DCR does no φ-form arithmetic.

Dimensions and theta graphs use this interface too. R matrices and twists retain `QPhase` objects with an explicit branch over x. A classical exact R value is an integer sign; a symbolic or exact-level R value remains a phase.

### Targets and scratch storage

`Level(k)`, `Exact()`, `Exact(k)`, `At(q)`, and `Classical()` select evaluation. Conflicting targets and keywords are rejected. `exact=true` requires a classical or level target; use `Symbolic()` for generic formulas.

`EvaluationWorkspace()` can be reused sequentially by supported scalar calls. Batches manage worker workspaces automatically. `eager=true` is deprecated and delegates to the standard rule evaluator; it does not enable a different fast path.
