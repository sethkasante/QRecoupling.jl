
[![DOI](https://zenodo.org/badge/1140106779.svg)](https://doi.org/10.5281/zenodo.19446095)
[![CI](https://github.com/sethkasante/QRecoupling.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sethkasante/QRecoupling.jl/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://sethkasante.github.io/QRecoupling.jl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# QRecoupling.jl

**Quantum recoupling, from fast numerical coefficients to readable exact formulas.**

Evaluate classical and quantum Clebsch–Gordan coefficients, 3j/6j symbols, fusion and braiding data, and finite q-factorial series in Julia. Use the same labels at `q = 1`, at a root-of-unity, or at a real or complex q. Retain exact level values and generic formulas in quantum integers, with optional polynomial and radical views.

`QRecoupling` provides local building blocks for angular momentum, spin networks, quantum topology, and fusion-based tensor networks. It does not contract full networks or assemble state sums automatically.

## Quick start

Numerical values and exact values need only QRecoupling. For nested radicals and expanded polynomials,
install the optional backend with `import Pkg; Pkg.add("Nemo")` and load `Nemo` alongside QRecoupling.

```julia
using QRecoupling

js = (1, 1, 1, 1, 1, 1)
q6j(js...)                      # classical value ≈ 1/6
q6j(Level(10), js...)           # value at q = exp(iπ/12)
q6j(At(0.8 + 0.2im), js...)     # numerical complex deformation

q6j(Exact(), js...)             # 1//6, exactly
q6j(Exact(3), js...)            # −1/[2]²  at level 3
iszero(q6j(Exact(10), 1, 3//2, 3//2, 2, 3//2, 3//2))   # true, proved

f = q6j(Symbolic(), js...)      # ([5] − 1)/[4]!  the formula at every q
f(0.8), f(1)                    # (0.1704…, 1//6)

using Nemo                      # the optional backend
radical(q6j(Exact(10), js...))  # (2√3 − 3)/3
xform(q6j(Exact(10), js...))   # (2x² − 7)/3  at x = 2cos(π/12)
xform(q6j(Symbolic(), js...))  # (x² − 3) / (x⁴ − 3x² + 2)
```

The spin labels are physical spins: use `1//2` for a half-integer (also accepts `0.5`). The same target interface applies to Clebsch–Gordan coefficients, 3j, F/G symbols, dimensions, and q-integer/factorial/binomial products.

`qcg` and `q3j` carry powers of q, which `Symbolic()` and `Exact(k)` print: `qcg(Symbolic(), 1//2, 1//2, 1//2, -1//2, 1)` is `q^(−1/2) √(1/[2])`. Their zero and equality tests are exact; arithmetic on these weighted level forms is not yet supported.

## What makes it useful?

- **Direct numerical kernels.** Factorial rules supply scaled ratios and adaptive precision without constructing or expanding a symbolic expression first.
- **Related coefficients share work.** Batches reuse tables; F matrices share recurrence coefficients across columns. Reusable workspaces reduce repeated scratch allocation.
- **Exact decisions without expansion.** Compact formulas support exact zero and equality tests without Nemo. Unweighted level forms also support exact arithmetic; polynomial expansion and nested radicals use the optional backend.
- **Explicit control of symbolic cost.** Large expressions stay deferred. Expansion and radical extraction are requests, not hidden requirements for computing a number.

```julia
using LinearAlgebra
F, e, f = fmatrix(Level(6), 1, 1, 1, 1)
@assert transpose(F) * F ≈ I

labels = [(j,j,j,j,j,j) for j in 1:5]
values = q6j(Level(20), labels)
grid = q6j(Level(15:20), labels)  # labels × levels
```

For parameter sweeps, each new level or q needs new tables. Large Float64 level and eligible positive-real-q recurrence tables use addition formulas with a small set of wide-precision seeds, reducing setup work. Reuse an `EvaluationWorkspace()` for scratch storage; see [Accuracy and performance](https://sethkasante.github.io/QRecoupling.jl/stable/performance/) for cold versus warm workloads.

At generic complex q, F matrices are complex orthogonal, not generally unitary: use `transpose(F)`, not the adjoint, in the algebraic identity.

`qcg(j1, m1, j2, m2, j, m)` is the quantum Clebsch–Gordan coefficient ⟨j₁m₁; j₂m₂|j m⟩ of U_q(sl₂), with `m = m1 + m2` by default; at q = 1 it is the ordinary coefficient:

```julia
qcg(1//2, 1//2, 1//2, -1//2, 1, 0)            # 1/√2
qcg((1//2, 1//2), (1//2, -1//2), (1, 0))      # the same, as (j, m) pairs
qcg(At(0.8), 1//2, 1//2, 1//2, -1//2, 1, 0)   # quantum value at q = 0.8
C, rows, cols = qcg_matrix(At(0.8), 1, 1)       # whole coupling matrix, by recurrence
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

This interface covers finite factorial-product sums. Arbitrary parameterized q-Pochhammer factors and infinite-series convergence are not currently supported.

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

See the [migration guide](https://sethkasante.github.io/QRecoupling.jl/stable/migration/) and [changelog](https://github.com/sethkasante/QRecoupling.jl/blob/main/CHANGELOG.md) for changes in old and current releases.

## Documentation and scope

- [Getting started](https://sethkasante.github.io/QRecoupling.jl/stable/getting_started/)
- [Tensor networks and modular data](https://sethkasante.github.io/QRecoupling.jl/stable/tutorials/tensor_networks/)
- [TensorKit sectors](https://sethkasante.github.io/QRecoupling.jl/stable/tutorials/tensorkit/) `SU2k(k)` and `SU2q(q)`, with TensorKitSectors loaded
- [Exact formulas](https://sethkasante.github.io/QRecoupling.jl/stable/tutorials/exact_forms/) and [identity checks](https://sethkasante.github.io/QRecoupling.jl/stable/tutorials/identities/)
- [Finite series](https://sethkasante.github.io/QRecoupling.jl/stable/tutorials/finite_series/)
- [Accuracy and performance](https://sethkasante.github.io/QRecoupling.jl/stable/performance/)
- [Research applications](https://sethkasante.github.io/QRecoupling.jl/stable/applications/)
- [Full documentation](https://sethkasante.github.io/QRecoupling.jl/stable/) (latest release; [`dev`](https://sethkasante.github.io/QRecoupling.jl/dev/) follows `main`)

**What the numbers mean:** Floating-point values, exact values, and zero tests are computed separately, each stated in the [accuracy guide](https://sethkasante.github.io/QRecoupling.jl/stable/performance/):

- *Floating point:* factorial-rule sums with compensated arithmetic; precision is raised when cancellation would otherwise cost significant digits. Final conversion can underflow a nonzero value to zero.
- *Exact:* a value is written in quantum integers from its formula and is never expanded to read it. Zeros and equality are decided from residues at word primes with a proved bound, never by a numerical tolerance. `float` of an exact form rounds from more digits than `Level(k)` keeps; `Level(k)` stays the direct route to a number.
- *Zero tests:* `iszero_at` confirms unresolved modular candidates with exact arithmetic. The `:cancels` entries of `level_spectrum` are screening results unless you pass `prove = true`.

The fast paths are Float64 and machine-size labels. High-degree exact sums and BigFloat F matrices are supported but cost more. Batches and level grids can run on several threads; do not clear caches while other threads are evaluating.

## Citation

To cite the package, use Seth K. Asante, [QRecoupling.jl: numerical and exact computation of classical and quantum recoupling coefficients](https://arxiv.org/abs/2610.24074), arXiv:2610.24074 (2026).

```bibtex
@article{Asante2026qrecoupling,
      title={QRecoupling.jl: numerical and exact computation of classical and quantum recoupling coefficients},
      author={Asante, Seth K.},
      year={2026},
      doi={10.48550/arXiv.2610.24074},
      eprint={2610.24074},
      archivePrefix={arXiv},
      primaryClass={cs.MS}
}
```
