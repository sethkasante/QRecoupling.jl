# Migration

## To v0.6

### Exact values are forms

`Exact()` is still the exact value at q = 1. `Exact(k)` and `Symbolic()` now return a [`QForm`](@ref): the value written in quantum integers, from its formula.

| Before | Now |
|---|---|
| `q6j(Symbolic(), js...)` was a `SymbolicValue` | A form, the formula at every q, with no value of its own. `f(q)` evaluates it, `f(1)` gives the exact value at q = 1, and `Float64(f)` and arithmetic raise an error. `xform` accepts it |
| `q6j(Exact(k), js...)` was an `ExactX`, and needed Nemo | A form at the level, with no backend. `==`, `iszero` and arithmetic are exact |
| `q3j` and `qcg` had no `Exact(k)` or `Symbolic()` | Both write them with their powers of q |
| `v.x_value`, `v.rad` | `xform(v)` to print it, `QRecoupling.xpolynomial(v)`, `QRecoupling.radicand(v)`, `radical(v)`, with Nemo |
| `ExactXSum(k)` to accumulate a sum | `sum`, or start from `zero(v)`. `ExactX` and `ExactXSum` are no longer exported; `QRecoupling.exact_x` still returns them, with Nemo |
| `isempty(residual)` | `iszero(residual)`, or `==` |

`q6j(js...; k = k, exact = true)` is the same as `q6j(Exact(k), js...)`, and `q6j(js...; exact = true)` the same as `q6j(Exact(), js...)`.

`==` on a `Symbolic()` value means equal at every q: `qint(Symbolic(), 6) == qfact(Symbolic(), 3)` and `qint(Symbolic(), 6) == 6` are false, although all are 6 at q = 1. `Exact()` compares values: `qint(Exact(), 6) == qfact(Exact(), 3)` is true.

### Optional backend

Nemo is no longer a required dependency. Install it with `import Pkg; Pkg.add("Nemo")` and load it with `using QRecoupling, Nemo` when using `radical`, `xform` or `prove_identity`. Both package load orders are supported. Missing-backend errors explain how to enable it.

`ExactX`, `ExactXSum`, and `XValue` now parameterize their polynomial payload. `QRecoupling.exact_x(rule, k)` still gives an `ExactX`, with its arithmetic, when Nemo is loaded.

### Removed after their deprecation

| Removed | Use instead |
|---|---|
| `Exact(k; form = :canonical)` and the `CompositeExactResult` values it returned | `Exact(k)` |
| `form = :x` | `Exact(k)`, which is the same thing |
| `eager = true` | Remove the keyword |
| `QRecoupling.xvalue(value)` with one argument | `xform(value)` |
| `x_form(f)` | `xform(f)`. `x_form` still works for now, with a deprecation warning |
| `phi_form(f)`, `splits_completely(f)` | Removed. `xform(f)` expands a formula in x, and `radical` gives nested roots |
| `qseries(z -> monomial, range)` | Removed. Use `qseries(factors, range)` or `FactorialSum`, with factorial triples |
| `CyclotomicMonomial`, `DCR` and their helpers, such as `build_dcr!`, `add_qfact!`, `build_series`, `qint_mono`, `qfact_mono`, `qbinomial_mono` | Removed. Use `FactorialSum` |
| `project_discrete`, `project_analytic`, `project_classical`, `project_classical_exact`, `value.dcr` | Removed. Use `qeval` with `Level(k)`, `At(q)`, `Classical()` or `Exact()` |
| `qeval(m; k, exact = true)` for a monomial or DCR | `qeval(Exact(k), rule)` with a `FactorialSum` |
| `QRecoupling.project_exact`, `QRecoupling.evaluate_exact` on the removed type | `Exact(k)` and `numeric_value` |

The three-argument `QRecoupling.xvalue(rad, num, den)` is unchanged.

Adding `Exact()` values with unlike roots gives an exact `ClassicalSum`, such as `1//2 + √15/6`, where it gave a `BigFloat`. Use `float` or `numeric_value` for a number.

Two numerical behaviours changed with the removal. `qint`, `qfact`, `qbinomial`, `qdim` and `theta_value` at q = ±i, which is the root of unity of level 0, now raise an error that names `Level(0)` unless the product is identically 0 or 1. They used to answer for some labels and raise for others. Values of these functions and of custom series at a level come from the same rule evaluator as the symbols, so a few last digits differ, and `T = Float32` is now rounded from `Float64` accuracy.

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
| Symbol interface, for adding symbols | `QSymbol`, `SixJ`, `ThreeJ`, `FSymbol`, `GSymbol`, `Tetrahedron`, `ThetaValue`, `symbol_rule`, `level_admissible`, `symbol_family`, `symbol_of`, `nlabels` |
| Lower-level exact access | `exact_x` (use `Exact(k)`), `generic_sixj` (use `xform(q6j(Symbolic(), …))`), `xvalue` (the `XValue` constructor; use `xform` to expand), `radical_form`, `has_radical_form`, `radical_levels`, `xpolynomial`, `radicand` (use `v.x_value`) |
| Modular data helpers | `level_labels`, `central_charge`, `total_qdim`, `gauss_sum`, `monodromy`, `verlinde` |
| Other | `EvalTarget` |

`smatrix`, `tmatrix`, `bmatrix`, `twist`, `rmatrix` and `fmatrix` remain exported.

### New

[`qcg`](@ref), [`qcg_matrix`](@ref) and [`qcg_row`](@ref) give coupling coefficients, whole coupling matrices and the coupled decomposition of a product state. `qcg` accepts its labels flat, `qcg(j1, m1, j2, m2, j, m)`, or as pairs, `qcg((j1, m1), (j2, m2), (j, m))`.

## To v0.4

### Numerical defaults now include the product helpers

Symbol functions, `qeval`, `qint`, `qfact`, and `qbinomial` evaluate classically when no q or k is supplied. Code that used their return values as symbolic objects must request `Symbolic()` explicitly.

```@example migration
using QRecoupling, Nemo
q6j(1,1,1,1,1,1)                      # approximately 1/6
q6j(Exact(),1,1,1,1,1,1)          # exact classical coefficient
s = q6j(Symbolic(),1,1,1,1,1,1)       # the formula at every q
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

`Symbolic()` writes a formula in quantum integers; it does not have the old monomial multiplication interface. Build product rules using factorial prefactors or expand with `xform`.

### Exact classical evaluation

`Exact()` is an alias for `Classical(exact=true)`, and it also works with the other recoupling symbols and product helpers.

`exact=true` and `Classical(exact=true)` remain supported. `Exact()` takes no form option: use `Exact(k)` for a level.

### Exact levels use x

In v0.4 and v0.5, `Exact(k)` and `k=10, exact=true` returned an `ExactX`, a polynomial in x. From v0.6 they return a form; see above. Read coefficients with `QRecoupling.xpolynomial` and `QRecoupling.radicand`, request nested square roots with `radical(v)`, and convert numerically with `Float64`, `BigFloat`, or `numeric_value`. `radical` returns `RadExpr` or a reasoned `NoRadical` result.

`Exact(k; form=:canonical)` retained the former cyclotomic representation until v0.6. `form=:deferred`, `ExactLevelWorkspace`, `canonicalize_exact`, and dedicated Biedenharn–Elliott verifier functions from development drafts are not current APIs. Write identities from the public exact arithmetic as in the [identity tutorial](tutorials/identities.md).
### Symbolic output is bounded and lazy

Use `xform(value)` to expand a symbolic rule. The earlier `QRecoupling.xvalue(value)` spelling was removed in v0.6; `value.x_value` on an exact level result instead exposes its stored polynomials.

`Symbolic()` displays the formula in quantum integers without carrying out the sum, for both small and large labels. `xform(value)` explicitly expands it in x, with radicals in ψ factors.

Dimensions and theta graphs use this interface too. R matrices and twists retain `QPhase` objects with an explicit branch over x. A classical exact R value or twist is an integer; a symbolic or exact-level one remains a `QPhase`.

### Targets and scratch storage

`Level(k)`, `Exact()`, `Exact(k)`, `At(q)`, and `Classical()` select evaluation. Conflicting targets and keywords are rejected. `exact=true` requires a classical or level target; use `Symbolic()` for generic formulas.

`EvaluationWorkspace()` can be reused sequentially by supported scalar calls. Batches manage worker workspaces automatically. `eager=true` is deprecated and delegates to the standard rule evaluator; it does not enable a different fast path.
