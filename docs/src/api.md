# API reference

```@meta
CurrentModule = QRecoupling
```

## Evaluation targets

Omitting q and k requests classical numerical evaluation. `Exact()` requests exact
classical rational/radical values, while `Symbolic()` retains the parameter. `Exact(k)` defaults to the real x basis;
`Symbolic()` retains a factorial rule. See [Getting started](getting_started.md) for target
selection and [Accuracy and performance](performance.md) for precision limits.

```@docs
EvalTarget
Classical
Level
Exact
At
Symbolic
```

## Recoupling and product functions

Public spin labels are physical spins. `qint`, `qfact`, and `qbinomial` also evaluate
classically by default in v0.4; use their `Symbolic()` methods for a generic expression.
The graph helpers below require the `QRecoupling.` qualifier.

```@docs
q6j
q3j
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
level_labels
twist
monodromy
central_charge
total_qdim
gauss_sum
verlinde
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

The [exact-forms tutorial](tutorials/exact_forms.md) documents `SymbolicValue`, `x_form`,
`XValue`, `ExactX`, `ExactXSum`, and `numeric_value` with examples. Coefficients and radical
extraction are accessible independently of the bounded display.

```@docs
exact_x
xpolynomial
radicand
radical_form
has_radical_form
radical_levels
RadExpr
generic_sixj
prove_identity
```

## Structural queries

Cancellation queries use modular screening. They are useful for exploration but should not be presented as universal exact-zero certificates. Confirm cancellation candidates with
exact arithmetic. See [Accuracy and performance](performance.md).

```@docs
iszero_at
issingular_at
level_spectrum
```

## Extending symbols

`SixJ`, `ThreeJ`, `FSymbol`, `GSymbol`, `Tetrahedron`, and `ThetaValue` identify supported
symbol families. The rule interface uses physical labels; some lower-level algebraic
constructors instead use doubled integers, as stated in their signatures.

```@docs
QSymbol
symbol_rule
nlabels
level_admissible
symbol_family
symbol_of
```

## Compatibility and low-level projections

DCR display is structural. These APIs support existing code and callback series; they are
not required to obtain a numerical symbol or a generic x-form. `project_exact` is the legacy
cyclotomic projector; prefer `Exact(k)` for current rule-backed exact level values.
`Exact(k; form=:canonical)` remains deprecated compatibility behavior. The former deferred
exact carrier and dedicated identity-verifier APIs are no longer part of the current package.

```@docs
qint_mono
qfact_mono
qbinomial_mono
CyclotomicMonomial
DCR
QPhase
symbolic_terms
phi_form
splits_completely
build_series
build_dcr!
add_qint!
add_qfact!
project_discrete
project_exact
project_classical
project_classical_exact
project_analytic
```
