# QRacahSymbols.jl

[![CI](https://github.com/sethkasante/QRacahSymbols.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sethkasante/QRacahSymbols.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://sethkasante.github.io/QRacahSymbols.jl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**QRacahSymbols.jl** is a high-performance, mathematically rigorous Julia library for evaluating Quantum $6j$-symbols, $3j$-symbols, and topological category data for the $SU(2)_k$ quantum group.

Standard hypergeometric evaluations of quantum Racah symbols suffer from catastrophic cancellation and floating-point overflow at high spins. This package solves this by introducing a highly optimized **Three-Tier Architecture**, allowing researchers to seamlessly seamlessly jump between nanosecond floating-point tensor contractions and zero-precision-loss algebraic number field evaluations.



## Core Features

* **The 3-Tier Engine**:
  * `Numeric`: Ultrafast table-based evaluation using log-sum-exp stabilization. Safe from factorial overflow up to massive spins ($j \approx 400$, $k > 10,000$).
  * `Exact`: Zero-precision-loss algebraic evaluation mapping strictly into `Nemo.jl` cyclotomic number fields ($\mathbb{Q}(\zeta_{2N})$).
  * `Generic`: Pure symbolic factorization into dense integer arrays of `CycloMonomial` representations.
* **TQFT Category Suite**: Natively evaluates Quantum Dimensions, R-matrices (braiding), F-symbols (fusion), and G-symbols (tetrahedral weights).
* **Geometric Canonicalization**: Internally exploits the 24-fold $S_4$ tetrahedral symmetry and massive LRU caches to bypass redundant computations in $O(1)$ time.

## Installation

```julia
# Press ']' in the Julia REPL to enter the package manager
pkg> add QRacahSymbols
