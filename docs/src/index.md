# QRecoupling.jl

**Quantum recoupling, from fast numerical coefficients to readable exact formulas.**

QRecoupling evaluates classical and quantum 3j and 6j symbols, fusion and braiding data,
and finite sums of symmetric q-factorials. Use the same labels at the classical limit,
at a root of unity, or at a real or complex deformation parameter. Request exact level
values in `x = q + q⁻¹`, or retain a factorial rule for symbolic exploration.

The package supplies local coefficients and change-of-basis matrices for angular-momentum
calculations, spin networks, Turaev–Viro models, and fusion-based tensor networks. It does
not assemble a triangulation or contract a complete tensor network for you.

## Start with one symbol

```@example home
using QRecoupling
classical = q6j(1, 1, 1, 1, 1, 1)
level = q6j(Level(10), 1, 1, 1, 1, 1, 1)
exact = q6j(Exact(10), 1, 1, 1, 1, 1, 1)
@assert isapprox(classical, 1/6) # hide
@assert isapprox(Float64(exact), level) # hide
exact
```

At level 10 this is `(2√3 − 3)/3`, also represented exactly by
`(2x² − 7)/3` with `x = 2cos(π/12)`.

```@example home
symbolic = q6j(Symbolic(), 1, 1, 1, 1, 1, 1)
x_form(symbolic)
```

Symbolic construction and display retain the factorial rule without carrying out its sum.
`x_form(symbolic)` explicitly expands it in x; large expanded formulas have a bounded display.

## Choose the calculation

| Task | Interface | Result |
|:--|:--|:--|
| Classical angular momentum | `q6j(js...)` | Numerical q = 1 value |
| Exact classical coefficients | `q6j(Exact(), js...)` | Rational/radical arithmetic |
| Root of unity | `q6j(Level(k), js...)` | Numerical q = exp(iπ/(k+2)) value |
| Exact level coefficient | `q6j(Exact(k), js...)` | Real x-basis value with radicals when needed |
| Real or complex deformation | `q6j(At(q), js...)` | Adaptive numerical evaluation |
| Generic formula | `q6j(Symbolic(), js...)` | Rule-backed symbolic x-form |
| Whole fusion transformation | `fmatrix(a,b,c,d; k=k)` | Matrix and its row/column labels |
| Finite factorial series | `FactorialSum(...)`, `qeval(...)` | The same evaluation targets |

## Why use it?

- **One mathematical rule, several evaluations.** Compact factorial rules feed direct
  numerical and exact calculations, avoiding large symbolic expansions on the numerical path.
- **Shared work for related symbols.** Batches reuse tables; F matrices share recurrence
  coefficients across columns instead of evaluating every entry as an independent sum.
- **Exact formulas you can inspect.** Generic expressions use x and ψ radical factors;
  level values reduce modulo the minimal polynomial of `2cos(π/(k+2))` using Nemo.
- **Control over expensive work.** Reuse workspaces, choose numerical precision, and
  request symbolic expansion or radical extraction explicitly.

See [Accuracy and performance](performance.md) for the numerical contracts and current
limitations, including the distinction between approximate checks and exact proofs.

## Install and learn

Julia 1.10 or later is supported. In Julia's package manager:

```text
pkg> add QRecoupling
```

This documentation describes the **v0.4 development/release-candidate API**. Until v0.4 is
registered, the registry may install v0.3.4; use the current checkout with `Pkg.develop(path=...)`
to run these examples. See [Migrating to v0.4](migration.md) for changed defaults.

- [Getting started](getting_started.md): labels, targets, and first calculations.
- [Recoupling symbols](tqft.md): normalization, dimensions, and phases.
- [Tensor networks and modular data](tutorials/tensor_networks.md): fusion bases and braiding.
- [Exact forms](tutorials/exact_forms.md) and [identity checks](tutorials/identities.md).
- [Finite series](tutorials/finite_series.md) and [factorial-rule architecture](series.md).
- [Applications](applications.md): research workflows built from these components.

## Citation

The original cyclotomic framework is described in Seth K. Asante,
[Deferred Cyclotomic Representation for Stable and Exact Evaluation of q-Hypergeometric Series](https://arxiv.org/abs/2604.13196).
The v0.4 factorial-rule kernels and x-form interface extend that framework.

```bibtex
@misc{Asante2026dcr,
  title={Deferred Cyclotomic Representation for Stable and Exact Evaluation of q-Hypergeometric Series},
  author={Seth K. Asante},
  year={2026},
  eprint={2604.13196},
  archivePrefix={arXiv},
  primaryClass={math-ph}
}
```
