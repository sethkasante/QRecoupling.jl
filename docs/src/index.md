# QRecoupling.jl

**Quantum recoupling, from fast numerical coefficients to readable exact formulas.**

QRecoupling evaluates classical and quantum Clebsch–Gordan coefficients, 3j and 6j symbols, fusion and braiding data, and finite sums of symmetric q-factorials. Use the same labels at the classical limit, at a root of unity, or at a real or complex deformation parameter. Request exact level values in `x = q + q⁻¹`, or retain a factorial rule for symbolic exploration.

The package supplies local coefficients and change-of-basis matrices for angular-momentum calculations, spin networks, Turaev–Viro models, and fusion-based tensor networks. It does not assemble a triangulation or contract a complete tensor network for you.

Nemo is optional. Numerical calculations, classical exact values, and unexpanded expressions
work without it. The exact level and polynomial examples below require `using Nemo`
after installing Nemo in your environment.

## Start with one symbol

```@example home
using QRecoupling, Nemo
classical_value = q6j(1, 1, 1, 1, 1, 1)
level = q6j(Level(10), 1, 1, 1, 1, 1, 1)
exact = q6j(Exact(10), 1, 1, 1, 1, 1, 1)
@assert isapprox(classical_value, 1/6) # hide
@assert isapprox(Float64(exact), level) # hide
exact
```

An exact value is written in quantum integers `[n]`, from the formula itself. `radical(exact)` gives `(2√3 − 3)/3`.

```@example home
symbolic = q6j(Symbolic(), 1, 1, 1, 1, 1, 1)
x_form(symbolic)
```

Symbolic construction and display retain the factorial rule without carrying out its sum. `x_form(symbolic)` explicitly expands it in x; large expanded formulas have a bounded display.

## Choose the calculation

| Task | Interface | Result |
|:--|:--|:--|
| Classical angular momentum | `q6j(js...)` | Numerical q = 1 value |
| Exact value at every q | `q6j(Exact(), js...)` | Quantum integers; `classical` gives the rational at q = 1 |
| Root of unity | `q6j(Level(k), js...)` | Numerical q = exp(iπ/(k+2)) value |
| Exact level value | `q6j(Exact(k), js...)` | Quantum integers, simplified, with zeros proved |
| Real or complex deformation | `q6j(At(q), js...)` | Adaptive numerical evaluation |
| Polynomials in x | `q6j(Symbolic(), js...)` | Rule-backed x-form, expanded by `x_form` with Nemo |
| Whole fusion transformation | `fmatrix(a,b,c,d; k=k)` | Matrix and its row/column labels |
| Finite factorial series | `FactorialSum(...)`, `qeval(...)` | The same evaluation targets |

Quantum CG and 3j coefficients support all of these except `Symbolic()`. Their exact values carry powers of q.

## Why use it?

- **One mathematical rule, several evaluations.** Compact factorial rules feed direct numerical and exact calculations, avoiding large symbolic expansions on the numerical path.
- **Shared work for related symbols.** Batches reuse tables; F matrices share recurrence coefficients across columns instead of evaluating every entry as an independent sum.
- **Exact formulas you can inspect.** Generic expressions use x and ψ radical factors; level values reduce modulo the minimal polynomial of `2cos(π/(k+2))` using Nemo.
- **Control over expensive work.** Reuse workspaces, choose numerical precision, and request symbolic expansion or radical extraction explicitly.

See [Accuracy and performance](performance.md) for the numerical contracts and current limitations, including the distinction between approximate checks and exact proofs.

## Install and learn

Julia 1.10 or later is supported. In Julia's package manager:

```text
pkg> add QRecoupling
```

To update an existing installation to the latest release, use `pkg> up QRecoupling`, and check the installed version with `pkg> status QRecoupling`. This documentation describes the latest release; see [Migration](migration.md) for changes between versions.

- [Getting started](getting_started.md): labels, targets, and first calculations.
- [Recoupling symbols](tqft.md): normalization, dimensions, and phases.
- [Tensor networks and modular data](tutorials/tensor_networks.md): fusion bases and braiding.
- [Exact forms](tutorials/exact_forms.md) and [identity checks](tutorials/identities.md).
- [Finite series](tutorials/finite_series.md) and [factorial-rule architecture](series.md).
- [Applications](applications.md): research workflows built from these components.

## Citation

To cite the package, use Seth K. Asante, [QRecoupling.jl: numerical and exact computation of classical and quantum recoupling coefficients](https://arxiv.org/abs/2610.24074), arXiv:2610.24074 (2026).

```bibtex
@article{Asanteqrecoupling,
  title={QRecoupling.jl: numerical and exact computation of classical and quantum recoupling coefficients},
  author={Asante, Seth K.},
  year={2026},
  doi={10.48550/arXiv.2610.24074},
  eprint={2610.24074},
  archivePrefix={arXiv},
  primaryClass={cs.MS}
}
```
