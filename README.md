
[![DOI](https://zenodo.org/badge/1140106779.svg)](https://doi.org/10.5281/zenodo.19446095)
[![CI](https://github.com/sethkasante/QRecoupling.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sethkasante/QRecoupling.jl/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://sethkasante.github.io/QRecoupling.jl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# QRecoupling.jl

**QRecoupling.jl** is a high-performance Julia library for the **stable and scalable evaluation of quantum recoupling coefficients and q-hypergeometric series**, designed to overcome computational limitations of direct numerical and symbolic methods.

It is designed to overcome fundamental limitations of direct numerical and symbolic evaluation of $q$-deformed symbols, including catastrophic cancellation, expression swell, and redundant computation.


> **Main idea:** separate algebraic structure from numerical evaluation.

All quantities are first represented symbolically in a **Deferred Cyclotomic Representation (DCR)**, where exact cancellations happen automatically. Only after the expression is maximally reduced is it projected into a target field (numeric, exact algebraic, or classical asymptotic limit).

---

## Key Features

### • Deferred Cyclotomic Representation (DCR)
Instead of expanding quantum factorials into massive rational polynomials, the DCR encodes $q$-hypergeometric series using the sparse integer exponents of their cyclotomic factorization:
$$\mathcal{M} = \sigma q^P \prod_d \Phi_d(q^2)^{e_d} $$
- Multiplication and division are reduced to highly efficient integer vector addition/subtraction.
- Perfect square roots (like those in $\Delta$-triangle coefficients) are extracted exactly at the exponent level, bypassing the need for algebraic field extensions.

---

### • Universal Projection Framework

A single compiled DCR object can be evaluated across multiple regimes without recomputation:

| Regime | Description |
|------|-------------|
| **Root of unity ($k$)** | Fast numerical evaluation using Log-Sum-Exp |
| **Exact algebraic** | Evaluation in cyclotomic field $(\mathbb{Q}(\zeta_{2h}))$ via `Nemo.jl` |
| **Complex analytic** | Efficient evaluation for $q \in \mathbb{C}$ |
| **Classical limit** | Exact $q \to 1$ evaluation (recover Ponzano-Regge amplitudes) |

### • HPC-Ready & Zero-Allocation

The package is designed to be thread-safe and implements **zero-allocation** loops during large numerical evaluations.

---

### • Extensible TQFT Toolkit

The framework natively supports:

- $6j$ symbols (Racah–Wigner)
- $3j$ symbols  
- $F$-symbols (fusion)
- $R$-matrices (braiding)
- $G$-symbols (tetrahedral weights)

and is designed to extend to more $q$-deformed tensors.


## Installation

```julia
# Press ']' in the Julia REPL to enter the package manager
pkg> add QRecoupling
```
--- 
## Quick Start

Evaluate the core quantum $6j$ and $3j$-symbols. The evaluation mode is controlled via keyword arguments, dynamically routing the computation to the most optimal engine.

### DCR algebraic object construction
If no evaluation parameters are passed, the package builds the parameter-independent DCR object:
```julia
using QRecoupling

j = 1

# 1. Deferred Graph Construction (Returns a DCR Object)
julia> dcr6j = q6j(j, j, j, j, j, j)
DCR (Deferred Cyclotomic Representation)
 ├─ Range    : 3:4
 ├─ Max Index: d = 5
 ├─ Radical  : 1
 ├─ Root     : q¹² Φ₂⁻⁴ Φ₃⁻² Φ₄⁻²
 ├─ Base Term: -q⁻⁶ Φ₂² Φ₃ Φ₄
 └─ Sequence : 1 update ratios {R_z}
```
This representation is exact, minimal, and independent of evaluation field


### Projections 
Project the same abstract symbol into your required target field using the `project_dcr` function:
```julia
# 1. projection into discrete level `k` (Float64 by default)
julia> project_dcr(dcr6j,k=10)
0.1547005383792515

julia> j=1; 

# full evaluation (constructs dcr object internally and then project)
julia> q6j(j, j, j, j, j, j, k=10)
0.1547005383792515

# 2. exact algebraic projection in cyclotomic fields (ζ)
julia> project_dcr(dcr6j, k=10, exact=true)
Exact Algebraic Result in ℚ(ζ₂₄):
  Value: (-2//3*ζ^6 + 4//3*ζ^2 - 1)

julia> q6j(j, j, j, j, j, j, k=10, exact=true)
Exact Algebraic Result in ℚ(ζ₂₄):
-2//3*ζ^6 + 4//3*ζ^2 - 1

# 3. generic complex q projection
julia> q6j(j, j, j, j, j, j, q=exp(0.5im))
0.035851185150113485 + 1.969762350587362e-17im

# 4. Classical Ponzano-Regge Limit (q -> 1, WignerSymbols) 
julia> q6j(j, j, j, j, j, j, q=1, exact=true)
1//6
```
### Topological Tensors

`QRecoupling.jl` provides direct APIs for constructing the composite tensors necessary for 3D state sums, automatically handling internal phase shifts and quantum dimensions.
```julia
julia> k = 5;
#quantum dimensions
julia> qdim(1/2,k=k,exact=true)
-ζ^5 + ζ^4 - ζ^3 + ζ^2 + 1

# R-Matrix braiding
julia> rmatrix(1, 1, 1, k=5)
-0.9009688679024191 + 0.4338837391175581im

# F-Symbol (fusion)
julia> fsymbol(1, 1, 1, 1, 1, 1, k=5)
0.19806226419516196

# G-Symbol (tetrahedral weight for Turaev-Viro invariant)
julia> gsymbol(1, 1, 1, 1, 1, 1, k=5)
1.0000000000000007
```

## Generic $q$-Series
`QRecoupling.jl` can also be use to study generic basic $q$-hypergeometric series. Here's how to construct a DCR for a custom sequence, such as $\sum_{z=1}^{10} (-1)^z [z]_q!$:
```julia

# build the q-series 
julia> custom_series = qseries(1:10) do z
           return (-1)^z * qfact(z)
       end

julia> qeval(custom_series, k=10)
10527.615497522727

julia> qeval(custom_series, q=0.05+0.95im)
-0.8168346401544203 - 0.22668133997266324im
```

## More features
- **Memory Management:** `QRecoupling.jl` caches cyclotomic tables and numeric workspaces to speed up parameter sweeps. You can manually flush these by calling `empty_caches!()`.
- **Exact Algebra Computations:** The `exact=true` flag for quantum symbols returns a `CompositeExactResult`. You can multiply these by raw integers, floats, or other exact symbols.

## Documentation

For the complete API reference, interactive tutorials, and deep dives into the mathematical architecture, please see the [Official Documentation](https://sethkasante.github.io/QRecoupling.jl/).

## Citation

If you use `QRecoupling.jl` in your research, please cite the mathematical framework behind the evaluation algorithm:

**Deferred Cyclotomic Representation for Stable and Exact Evaluation of q-Hypergeometric Series**
Seth K. Asante (2026). *arXiv preprint arXiv:2604.13196*.

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

