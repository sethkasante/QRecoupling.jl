# API reference

```@meta
CurrentModule = QRecoupling
```

## Evaluation targets

Omitting q and k requests classical numerical evaluation. `Exact()` requests exact classical rational/radical values, while `Symbolic()` retains the parameter. `Exact(k)` defaults to the real x basis; `Symbolic()` retains a factorial rule. See [Getting started](getting_started.md) for target selection and [Accuracy and performance](performance.md) for precision limits.

```@docs
QRecoupling.EvalTarget
Classical
Level
Exact
At
Symbolic
```

## Recoupling and product functions

Public spin labels are physical spins. `qint`, `qfact`, and `qbinomial` evaluate classically by default since v0.4; use their `Symbolic()` methods for a generic expression. The graph helpers below require the `QRecoupling.` qualifier.

```@docs
q6j
q3j
QRecoupling.q3j_factorial
qcg
qcg_matrix
qcg_row
fsymbol
gsymbol
qdim
rmatrix
QRecoupling.tetrahedron
QRecoupling.theta_value
qint
qfact
qbinomial
```

## Fusion and modular data

```@docs
fmatrix
fmatrix_labels
bmatrix
smatrix
tmatrix
QRecoupling.level_labels
twist
QRecoupling.monodromy
QRecoupling.central_charge
QRecoupling.total_qdim
QRecoupling.gauss_sum
QRecoupling.verlinde
```

## Rules, series, and repeated evaluation

```@docs
AffineFactorial
FactorialSum
qseries
qeval
EvaluationWorkspace
all_6j
empty_caches!
```

## Symbolic and exact x forms

The [exact-forms tutorial](tutorials/exact_forms.md) documents `SymbolicValue`, `x_form`, `XValue`, `ExactX`, `ExactXSum`, and `numeric_value` with examples. Coefficients and radical extraction are accessible independently of the bounded display.

```@docs
QRecoupling.exact_x
QRecoupling.xpolynomial
QRecoupling.radicand
QRecoupling.radical_form
QRecoupling.has_radical_form
QRecoupling.radical_levels
RadExpr
QRecoupling.generic_sixj
prove_identity
```

## Structural queries

`iszero_at` confirms unresolved cancellation candidates with exact arithmetic. The `:cancels` entries of `level_spectrum` are modular screening results unless you pass `prove = true`. See [Accuracy and performance](performance.md).

```@docs
iszero_at
issingular_at
level_spectrum
```

## Extending symbols

`SixJ`, `ThreeJ`, `FSymbol`, `GSymbol`, `Tetrahedron`, and `ThetaValue` identify supported symbol families. The rule interface uses physical labels; some lower-level algebraic constructors instead use doubled integers, as stated in their signatures.

```@docs
QRecoupling.QSymbol
QRecoupling.symbol_rule
QRecoupling.nlabels
QRecoupling.level_admissible
QRecoupling.symbol_family
QRecoupling.symbol_of
```

## Compatibility and low-level projections

DCR display is structural. These APIs support existing code and callback series and are planned for retirement; they are not required to obtain a numerical symbol or a generic x-form. `project_exact` is the legacy cyclotomic projector; prefer `Exact(k)` for current rule-backed exact level values. `Exact(k; form=:canonical)` remains deprecated compatibility behavior. The former deferred exact carrier and dedicated identity-verifier APIs are no longer part of the current package.

```@docs
QRecoupling.qint_mono
QRecoupling.qfact_mono
QRecoupling.qbinomial_mono
QRecoupling.CyclotomicMonomial
QRecoupling.DCR
QPhase
QRecoupling.symbolic_terms
phi_form
QRecoupling.splits_completely
QRecoupling.build_series
QRecoupling.build_dcr!
QRecoupling.add_qint!
QRecoupling.add_qfact!
QRecoupling.project_discrete
QRecoupling.project_exact
QRecoupling.project_classical
QRecoupling.project_classical_exact
QRecoupling.project_analytic
```
