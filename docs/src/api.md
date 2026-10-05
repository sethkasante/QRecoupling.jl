# API reference

```@meta
CurrentModule = QRecoupling
```

This reference covers the exported interface and two coefficient-access helpers. See [Migration](migration.md) for older interfaces and [Factorial rules and evaluation architecture](series.md) for implementation and extension details.

## Evaluation targets

Omitting q and k requests classical numerical evaluation. `Exact()` requests exact classical rational/radical values, `Exact(k)` uses the real x basis at a level, and `Symbolic()` retains a factorial rule. See [Getting started](getting_started.md) for target selection and [Accuracy and performance](performance.md) for precision limits.

```@docs
Classical
Level
Exact
At
Symbolic
```

## Recoupling and product functions

Public spin labels are physical spins. `qint`, `qfact`, and `qbinomial` evaluate classically by default; use their `Symbolic()` methods for a generic expression.

```@docs
q6j
q3j
qcg
qcg_matrix
qcg_row
fsymbol
gsymbol
qdim
rmatrix
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
twist
QPhase
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

The [exact-forms tutorial](tutorials/exact_forms.md#API-reference) documents `SymbolicValue`, `x_form`, `XValue`, `ExactX`, `ExactXSum`, `radical`, `NoRadical`, and `numeric_value` with examples. Use `phi_form` for cyclotomic factorization and `prove_identity` to compare generic expressions.

```@docs
RadExpr
phi_form
prove_identity
```

The following helpers return polynomial coefficients without expanding the display. They require the `QRecoupling.` qualifier.

```@docs
QRecoupling.xpolynomial
QRecoupling.radicand
```

## Structural queries

`iszero_at` confirms unresolved cancellation candidates with exact arithmetic. The `:cancels` entries of `level_spectrum` are modular screening results unless you pass `prove = true`. See [Accuracy and performance](performance.md).

These queries do not currently support the weighted `q3j` and `qcg` coefficients.

```@docs
iszero_at
issingular_at
level_spectrum
```
