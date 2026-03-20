# QRacahSymbols.jl

[![CI](https://github.com/sethkasante/QRacahSymbols.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sethkasante/QRacahSymbols.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://sethkasante.github.io/QRacahSymbols.jl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**QRacahSymbols.jl** is a high-performance, mathematically rigorous Julia library for the exact and numerically stable evaluation of quantum $\text{SU(2)}_k$ recoupling coefficients ($\{3j\}$- and $\{6j\}$-symbols and topological category) at arbitrary roots of unity. This package bridges the gap between **ultra-fast floating-point numerics** and **rigorous Computer Algebra System (CAS) exactness**.

By representing quantum prime factorizations as deferred sparse cyclotomic arrays, `QRacahSymbols.jl` completely bypasses floating-point overflow and traditional algebraic bottlenecks, allowing for blazing-fast discrete evaluations and exact topological proofs.

## Installation

```julia
# Press ']' in the Julia REPL to enter the package manager
pkg> add QRacahSymbols
```

## Core Features

`QRacahSymbols.jl` is built as a unified interface with highly specialized computational engines. 

* **Multi-Paradigm Engine**: Compute symbols dynamically as fast floats, rigorous cyclotomic algebraic numbers, or exact zero-allocation rationals.
* **Overflow Immunity**: Uses a specialized Log-Sum-Exp architecture and GMP C-calls to handle massive topological spins without NaN or Inf corruption.
* **Zero-Division CAS**: Extracts perfect algebraic squares symbolically before summation, providing division-free loops in exact cyclotomic fields.
* **TQFT Category Suite**: Support for quantum dimensions, R-matrices (braiding), F-symbols (fusion), and G-symbols (tetrahedral weights).


## Quick Start

Evaluate the core quantum $6j$ and $3j$ (Racah-Wigner) symbols. The **mode** keyword controls the computational engine. Output values may vary slightly based on the floating-point precision.
```julia
using QRacahSymbols

j = 1  # Spins 
k = 10 # Level k 

# 1. Fast numeric evaluation (default)
julia> q6j(j, j, j, j, j, j, k; mode=:numeric)
0.1547005383792515

# 2. Exact algebraic evaluation in rational cyclotomic fields (ζ)
julia> q6j(j, j, j, j, j, j, k; mode=:exact)
Exact SU(2)₁₀ Symbol:
  Value: -2//3*ζ^6 + 4//3*ζ^2 - 1

# 3. Classical Ponzano-Regge Limit (q -> 1 or k -> ∞)
julia> q6j(j, j, j, j, j, j; mode=:classical_exact)
1//6
```
### The Polymorphic Architecture
QRacahSymbols.jl utilizes different mathematically optimized engines invoked via the mode flag:
* `mode=:numeric`: The absolute speed limit. Uses pre-cached dense logarithmic arrays and a Log-Sum-Exp hot loop to rapidly compute single-shot floating point values for $\text{SU(2)}_k$.
* `mode=:exact`: The mathematician's tool. Maps the symbol into rigorous Nemo.jl cyclotomic fields. By extracting algebraic perfect squares symbolically, it provides $\mathcal{O}(1)$ radical extraction and impenetrable precision.
* `mode=:classical_exact`: Evaluates the un-deformed Ponzano-Regge limit. Bypasses Julia's Garbage Collector entirely by executing GMP integer math strictly in-place, dramatically outperforming standard libraries.
* `mode=:cyclo`: Builds the unevaluated algebraic graph of the symbol. Ideal for structural analysis of the hypergeometric ratio sequence.

### CycloMonomial Ouput

When operating in `:cyclo` mode, the engine dynamically collapses perfect squares to provide a human-readable, minimal representation of the symbol.
```julia
using QRacahSymbols

#Symbolic computation (very fast O(1) algebraic additions)
julia> res = q6j(1, 1, 1, 1, 1, 1; mode=:cyclo);
CycloResult (Hypergeometric Ratio Form)
  ├─ Max Φ_d(q) required : d = 5
  ├─ Overall Prefactor : -z⁶ Φ₂⁻² Φ₃⁻¹ Φ₄⁻¹
  └─ Ratios (R)     : 1 terms
       └─ R₁ : -z⁻⁴ Φ₅
```

### Topological Tensors

`QRacahSymbols.jl` provides direct APIs for constructing the composite tensors necessary for 3D state sums, automatically handling internal phase shifts and quantum dimensions.
```julia
julia> k = 5;
#quantum dimensions
julia> qdim(1/2, k; mode=:exact)
Exact SU(2)₅ Symbol:
  Value: -ζ^5 + ζ^4 - ζ^3 + ζ^2 + 1

# R-Matrix (braiding)
julia> rmatrix(1, 1, 1, k; mode=:numeric)
-0.6234898018587336 + 0.7818314824680298im

# F-Symbol (fusion)
julia> fsymbol(1, 1, 1, 1, 1, 1, k; mode=:numeric)
0.19806226419516174

# G-Symbol (tetrahedral weight for Turaev-Viro invariant)
julia> gsymbol(1, 1, 1, 1, 1, 1, k; mode=:numeric)
1.0000000000000002
```

## Memory Management

`QRacahSymbols.jl` caches prime sieves and cyclotomic phases under the hood when computing massive state sums over complex triangulations. If you change topological levels (k) drastically and need to free up RAM, you can clear the global caches dynamically.
```julia
clear_caches!()         # Clears all internal caches
clear_numeric_caches!() # Clears the Log-Sum-Exp tables
clear_exact_caches!()   # Clears dense Nemo.jl polynomials
```

## Documentation

For the complete API reference, interactive tutorials, and deep dives into the mathematical architecture, please see the [Official Documentation](https://sethkasante.github.io/QRacahSymbols.jl/).

## Citation

If you use `QRacahSymbols.jl` in your research, please cite the following paper:

> Seth K. Asante, *"QRacahSymbols.jl: Fast Exact Evaluation of Quantum Recoupling Coefficients using Cyclotomic Factorization"* (To appear soon, 2026). arXiv:XXXX.XXXXX.
