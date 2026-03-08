# QRacahSymbols.jl

[![CI](https://github.com/sethkasante/QRacahSymbols.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sethkasante/QRacahSymbols.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://sethkasante.github.io/QRacahSymbols.jl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**QRacahSymbols.jl** is a high-performance, mathematically rigorous Julia library for evaluating Quantum $6j$-symbols, $3j$-symbols, and topological category data for the SU(2)$_k$ quantum group.

Standard hypergeometric evaluations of quantum Racah symbols suffer from catastrophic cancellation and floating-point overflow at high spins. This package solves this by introducing a highly optimized **Three-Tier Architecture**, allowing researchers to seamlessly seamlessly jump between nanosecond floating-point tensor contractions and zero-precision-loss algebraic number field evaluations.

Standard hypergeometric evaluations of quantum Racah symbols suffer from catastrophic cancellation and floating-point overflow at high spins. **QRacahSymbols.jl** solves this by introducing a highly optimized **Three-Tier Architecture**, allowing researchers to seamlessly jump between nanosecond floating-point tensor contractions and zero-precision-loss algebraic number field evaluations.

## Core Features

* **The 3-Tier Engine**:
  * `Numeric`: Ultrafast table-based evaluation using log-sum-exp stabilization. Safe from factorial overflow up to massive spins ($j \approx 450$, $k > 20,000$).
  * `Exact`: Zero-precision-loss algebraic evaluation mapping strictly into `Nemo.jl` cyclotomic number fields ($\mathbb{Q}(\zeta_{2N}$).
  * `Generic`: Pure symbolic factorization into dense integer arrays of `CycloMonomial` representations.
* **TQFT Category Suite**: Natively evaluates Quantum Dimensions, R-matrices (braiding), F-symbols (fusion), and G-symbols (tetrahedral weights).
* **Geometric Canonicalization**: Internally exploits the 24-fold $S_4$ tetrahedral symmetry and massive LRU caches to bypass redundant computations in $O(1)$ time.

## Installation

```julia
# Press ']' in the Julia REPL to enter the package manager
pkg> add QRacahSymbols
```

## Quick Start

The master API intelligently routes your computation based on the provided mode. The level $k$ is passed as the final positional argument.
```julia
using QRacahSymbols

j = 2.5 # Spins are standard Reals (Float64, Int, or Rational)
k = 10  # Level of the SU(2)_k quantum group

# 1. High-Speed Numeric Evaluation (Default)
val_num = q6j(j, j, j, j, j, j, k; mode=:numeric) 

# 2. Exact Algebraic Evaluation in Cyclotomic Fields
val_exact = q6j(j, j, j, j, j, j, k; mode=:exact) 

# 3. Pure Symbolic Cyclotomic Polynomial Representation
val_symb = q6j(j, j, j, j, j, j; mode=:generic) # k-independent

# 4. Classical Ponzano-Regge Limit (q -> 1)
val_class = q6j(j, j, j, j, j, j; mode=:classical)

# F-Matrix (Fusion)
f_val = fsymbol(1, 1, 1, 1, 1, 1, k; mode=:numeric)

# R-Matrix (Braiding)
r_val = rmatrix(1, 1, 1, k; mode=:exact)

# Quantum Dimension
d_val = qdim(2.5, k; mode=:numeric, T=BigFloat) # Enforce BigFloat precision
```

## Evaluating TQFT Data
QRacahSymbols natively supports the local topological data required for state-sum invariants and tensor networks.
```julia
# F-Matrix (Fusion)
f_val = fsymbol(1, 1, 1, 1, 1, 1, k; mode=:numeric)

# R-Matrix (Braiding)
r_val = rmatrix(1, 1, 1, k; mode=:exact)

# Quantum Dimension
d_val = qdim(2.5, k; mode=:numeric, T=BigFloat) # Enforce BigFloat precision
```

## Architecture & Performance
Under the hood, QRacahSymbols.jl completely avoids heap allocations in its deep loops. The symbolic engine relies on in-place array mutation and the Hypergeometric Ratio Method to drop complexities from $O(N \log N)$ to $O(N)$, while the Exact engine natively builds quantum factorials directly in the algebraic field, eliminating generic array instantiation.

## Citation
If you use `QRacahSymbols.jl` in your research, please cite .... 