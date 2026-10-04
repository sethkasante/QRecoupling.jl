
[![DOI](https://zenodo.org/badge/1140106779.svg)](https://doi.org/10.5281/zenodo.19446095)
[![CI](https://github.com/sethkasante/QRecoupling.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sethkasante/QRecoupling.jl/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://sethkasante.github.io/QRecoupling.jl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# QRecoupling.jl

**Quantum recoupling, from fast numerical coefficients to readable exact formulas.**

Evaluate classical and quantum Clebsch–Gordan coefficients, 3j/6j symbols, fusion and braiding data, and finite q-factorial series in Julia. Use the same labels at `q = 1`, at a root-of-unity, or at a real or complex q. Build exact level values and generic formulas in `x = q + q⁻¹`.

`QRecoupling` provides local building blocks for angular momentum, spin networks, quantum topology, and fusion-based tensor networks. It does not contract full networks or assemble state sums automatically.

## Quick start

```julia
using QRecoupling

js = (1, 1, 1, 1, 1, 1)
q6j(js...)                      # classical value ≈ 1/6
q6j(Exact(), js...)             # exact classical rational = 1//6
q6j(Level(10), js...)           # value at q = exp(iπ/12)
v = q6j(Exact(10), js...)       # x-form: (2x² − 7)/3
radical(v)                     # explicitly request (2√3 − 3)/3
q6j(At(0.8 + 0.2im), js...)     # numerical complex deformation

s = q6j(Symbolic(), js...)      # symbolic form, unexpanded sum
x_form(s)                      # explicit x-form = (x² − 3) / (x⁴ − 3x² + 2)
```

The spin labels are physical spins: use `1//2` for a half-integer (also accepts `0.5`). The same target interface applies to Clebsch–Gordan coefficients, 3j, F/G symbols, dimensions, and q-integer/factorial/binomial products.

## What makes it useful?

- **Direct numerical kernels.** Factorial rules supply scaled ratios and adaptive precision without constructing or expanding a symbolic expression first.
- **Related coefficients share work.** Batches reuse tables; F matrices share recurrence coefficients across columns. Reusable workspaces reduce repeated scratch allocation.
- **Exact algebra in a readable basis.** Nemo-backed level values use the minimal polynomial of `2cos(π/(k+2))`; generic expressions use x with ψ factors under radicals.
- **Explicit control of symbolic cost.** Large expressions stay deferred. Expansion and radical extraction are requests, not hidden requirements for computing a number.

```julia
using LinearAlgebra
F, e, f = fmatrix(1, 1, 1, 1; k=6)
@assert transpose(F) * F ≈ I

labels = [(j,j,j,j,j,j) for j in 1:5]
values = q6j(Level(20), labels)
grid = q6j(Level(15:20), labels)  # labels × levels
```

At generic complex q, F matrices are complex orthogonal, not generally unitary: use `transpose(F)`, not the adjoint, in the algebraic identity.

`qcg(j1, m1, j2, m2, j, m)` is the quantum Clebsch–Gordan coefficient ⟨j₁m₁; j₂m₂|j m⟩ of U_q(sl₂), with `m = m1 + m2` by default; at q = 1 it is the ordinary coefficient:

```julia
qcg(1//2, 1//2, 1//2, -1//2, 1, 0)            # 1/√2
qcg((1//2, 1//2), (1//2, -1//2), (1, 0))      # the same, as (j, m) pairs
qcg(At(0.8), 1//2, 1//2, 1//2, -1//2, 1, 0)   # quantum value at q = 0.8
C, rows, cols = qcg_matrix(1, 1; q=0.8)       # whole coupling matrix, by recurrence
@assert C * transpose(C) ≈ I
```

## Custom finite series

```julia
n = 5
s = FactorialSum(0:n; prefactor=[n=>2],
                 factors=[(1,0,-2), (-1,n,-2)])
# [n]!² ∑ⱼ 1/([j]!² [n-j]!²)
@assert qeval(s) ≈ binomial(2n,n)
qeval(s; q=0.8)
qeval(Exact(20), s)
```

This interface covers finite factorial-product sums. Arbitrary parameterized q-Pochhammer factors and infinite-series convergence are not currently supported. A compatibility DCR callback interface remains available for cyclotomic monomials and explicit q powers.

## Install and upgrade version

Requires Julia 1.10 or later. In the package manager (press `]` at the Julia prompt):

```text
pkg> add QRecoupling
```

If you already have it, update to the latest release — older versions contain sign errors fixed in v0.3.4 — and check which version you have:

```text
pkg> up QRecoupling
pkg> status QRecoupling
```

Coming from v0.3: omitted q/k now means classical evaluation, including `qint`, `qfact`, and `qbinomial`. Use `Exact()` for exact classical values (rationals) and `Symbolic()` for generic output; `exact=true` is still supported. `Exact(k)` returns an `ExactX`, and `eager=true` is deprecated. DCRs remain available through `.dcr` and display only their structure. See the [migration guide](https://sethkasante.github.io/QRecoupling.jl/stable/migration/) and [changelog](https://github.com/sethkasante/QRecoupling.jl/blob/main/CHANGELOG.md).

## Documentation and scope

- [Getting started](https://sethkasante.github.io/QRecoupling.jl/stable/getting_started/)
- [Tensor networks and modular data](https://sethkasante.github.io/QRecoupling.jl/stable/tutorials/tensor_networks/)
- [Exact formulas](https://sethkasante.github.io/QRecoupling.jl/stable/tutorials/exact_forms/) and [identity checks](https://sethkasante.github.io/QRecoupling.jl/stable/tutorials/identities/)
- [Finite series](https://sethkasante.github.io/QRecoupling.jl/stable/tutorials/finite_series/)
- [Accuracy and performance](https://sethkasante.github.io/QRecoupling.jl/stable/performance/)
- [Research applications](https://sethkasante.github.io/QRecoupling.jl/stable/applications/)
- [Full documentation](https://sethkasante.github.io/QRecoupling.jl/stable/) (latest release; [`dev`](https://sethkasante.github.io/QRecoupling.jl/dev/) follows `main`)

**What the numbers mean.** Floating-point values, exact values, and zero tests are computed separately, each stated in the [accuracy guide](https://sethkasante.github.io/QRecoupling.jl/stable/performance/):

- *Floating point:* factorial-rule sums with compensated arithmetic; precision is raised when cancellation in the sum would otherwise cost digits loss.
- *Exact:* algebraic numbers in the x-basis; equality, including `ExactXSum`, is decided by algebraic arithmetic, never by a numerical tolerance.
- *Zero tests:* `iszero_at` confirms unresolved modular candidates with exact arithmetic. The `:cancels` entries of `level_spectrum` are screening results unless you pass `prove = true`.

The fast paths are Float64 and machine-size labels. High-degree exact sums and BigFloat F matrices are supported but cost more. Batches and level grids can run on several threads; do not clear caches while other threads are evaluating.

## Citation

For the original cyclotomic framework, cite Seth K. Asante,
[Deferred Cyclotomic Representation for Stable and Exact Evaluation of q-Hypergeometric Series](https://doi.org/10.21468/SciPostPhysCore.9.3.059), SciPost Phys. Core **9**, 059 (2026). The v0.4 factorial-rule numeric and exact implementations extend that framework.

```bibtex
@article{Asante2026dcr,
      title={Deferred cyclotomic representation for stable and exact evaluation of q-hypergeometric series},
      author={Asante, Seth K.},
      journal={SciPost Phys. Core},
      volume={9},
      pages={059},
      year={2026},
      doi={10.21468/SciPostPhysCore.9.3.059},
      eprint={2604.13196},
      archivePrefix={arXiv}
}
```
