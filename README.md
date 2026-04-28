
[![DOI](https://zenodo.org/badge/1140106779.svg)](https://doi.org/10.5281/zenodo.19446095)
[![CI](https://github.com/sethkasante/QRecoupling.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sethkasante/QRecoupling.jl/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://sethkasante.github.io/QRecoupling.jl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- > **Notice of Migration:** > `QRecoupling.jl` is the official, vastly expanded successor to the deprecated package `QRacahSymbols.jl`. The package has been completely rewritten and renamed to reflect its broader scope in quantum recoupling coefficients, $q$-hypergeometric series and their applications.

---  -->
# QRecoupling.jl

**QRecoupling.jl** is a high-performance Julia library for the **stable and scalable evaluation of quantum recoupling coefficients and q-hypergeometric series**, designed to overcome computational limitations of direct numerical and symbolic methods.

It is designed to overcome fundamental limitations of direct numerical and symbolic evaluation, including catastrophic cancellation, expression swell, and redundant computation.


> **Main idea:** separate algebraic structure from numerical evaluation.

All quantities are first represented symbolically in a **Deferred Cyclotomic Representation (DCR)**, and only later projected into a target field (numeric, exact, or asymptotic).

---

## Key Features

### • Deferred Cyclotomic Representation (DCR)
A sparse, combinatorial encoding of $q$-hypergeometric series with building blocks 
$$\mathcal{M} = \sigma q^P \prod_d \Phi_d(q^2)^{e_d} $$

- the representation is performed using the exponents of the cyclotomic basis {$q,\Phi_d(q^2)$} 
- multiplication and division translates to addition of exponents 

---

### • Universal Projection Framework

A single DCR object can be evaluated in multiple regimes:

| Regime | Description |
|------|-------------|
| **Root of unity ($k$)** | Fast numerical evaluation using Log-Sum-Exp |
| **Exact algebraic** | Evaluation in cyclotomic fields ($\zeta$) via `Nemo.jl` |
| **Complex analytic** | Efficient evaluation for $q \in \mathbb{C}$ |
| **Classical limit** | Exact $q \to 1$ reduction |

---

### • Extensible TQFT Toolkit

The framework naturally supports:

- $6j$ symbols (Racah–Wigner)
- $3j$ symbols  
- $F$-symbols (fusion)
- $R$-matrices (braiding)
- $G$-symbols (tetrahedral weights)

and is designed to extend to general tensor network observables.


## Installation

```julia
# Press ']' in the Julia REPL to enter the package manager
pkg> add QRecoupling
```
--- 
## Quick Start

Evaluate the core quantum $6j$ and $3j$ (Racah-Wigner) symbols. The **mode** keyword controls the computational engine. Output values may vary slightly based on the floating-point precision.

### DCR algebraic object construction
```julia
using QRecoupling

j = 1
k = 10 

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
For projection to target fields:
```julia
julia> j = 1;
# 1. projection into discrete level `k` (Float64 by default)
julia> q6j(j, j, j, j, j, j, k=10)
0.1547005383792515

# 2. exact algebraic projection in cyclotomic fields (ζ)
julia> q6j(j, j, j, j, j, j, k=10, exact=true)
Exact SU(2)₁₀ Symbol: 
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
`QRecoupling.jl` can also be use to study generic basic $q$-hypergeometric series. Here's how to construct a DCR for a custom sequence, such as $\sum_{z=1}^{10} [z]_q!$:
```julia

# build the q-series 
julia> custom_series = build_dcr(
           b -> nothing,                       
           (b, z) -> add_qfact!(b, z),      # base term at z_min   
           (b, z) -> add_qint!(b, z + 1),   # ratio: T_{z+1}/T_z = [z+1]_q   
           1, 10
          )
julia> project_dcr(custom_series,k=10)
29948.709646825515

julia> project_dcr(custom_series,q=0.5im)
1.0996624874742891e10 - 4.482327660636685e12im
```

## Cache Management
```julia
empty_caches!() # Clears all internal caches
```
Clears all internal caches (cyclotomic tables, numeric workspaces, etc.).

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

