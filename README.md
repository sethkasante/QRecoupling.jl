# QRacahSymbols.jl

[![CI](https://github.com/sethkasante/QRacahSymbols.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sethkasante/QRacahSymbols.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://sethkasante.github.io/QRacahSymbols.jl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**QRacahSymbols.jl** is a high-performance, mathematically rigorous Julia library for the exact and numerically stable evaluation of quantum $\text{SU(2)}_k$ recoupling coefficients ($\{3j\}$- and $\{6j\}$-symbols and topological category) at arbitrary roots of unity.

## The Catastrophic Cancellation Problem (and Solution)
Standard numerical implementations of the quantum Racah formula suffer from ***catastrophic cancellation*** at macroscopic spins due to the alternating summation of massive $q$-deformed factorials. This package (`QRacahSymbols.jl`) circumvents this crisis entirely using a novel **CycloMonomial architecture**, which rigorously factors the topological invariant into integer arrays of cyclotomic polynomials $(\Phi_d(q))$ prior to any floating-point evaluation.

This guarantees mathematically exact, zero-drift computations of tensor network amplitudes even in the deep semiclassical limit ($j \ge 500$), unlocking extreme-spin probes for research in mathematics and physics.

## Core Features

`QRacahSymbols.jl` is built as a unified interface with highly specialized computational engines:

* **The Polymorphic Engine**:
  * `:numeric`: Ultrafast evaluation using a *log-sum-exp* stabilization procedure for low spins. Safe from factorial overflow up to massive spins ($j \le 450$, for $k > 20,000$) using `BigFloat`. 
  * `:exact`: Algebraic evaluation mapping strictly into `Nemo.jl` cyclotomic number fields $\mathbb Q(\zeta_{2(k+2)})$.
  * `:generic`: Pure symbolic factorization into dense integer arrays of `CycloMonomial` representations.
  * `:classical`: Evaluates the un-deformed Ponzano-Regge limit ($q \to 1$ or $k \to \infty$). 
  
* **TQFT Category Suite**: Evaluates quantum dimensions, R-matrices (braiding), F-symbols (fusion), and G-symbols (tetrahedral weights).


## Installation

```julia
# Press ']' in the Julia REPL to enter the package manager
pkg> add QRacahSymbols
```

## Quick Start

The package exposes a unified interface through its polymorphic engines, allowing you to seamlessly switch between high-throughput floats, exact algebraic rings, and symbolic arrays.

### 1. Standard Evaluations
For low spins, the `:numeric` engine provides maximum speed using stack-allocated Float64 arithmetic, while the `:exact` engine gives mathematically rigorous results in cyclotomic fields:
```julia
using QRacahSymbols

j = 1 # Spins are standard Reals (Float64, Int, or Rational)
k = 10  # Level of the SU(2)_k quantum group

# 1. Fast numeric evaluation (default)
julia> q6j(j, j, j, j, j, j, k; mode=:numeric)
0.1547005383792515

# 2. Exact algebraic evaluation in rational cyclotomic Fields
julia> q6j(j, j, j, j, j, j, k; mode=:exact)
Exact SU(2)₁₀ Symbol:
  Prefactor(Δ²): 65//3*ζ^6 - 130//3*ζ^2 + 1351//36
  Racah Sum(Σ):  -14*ζ^6 + 28*ζ^2 + 24 

# 3. Symbolic cyclotomic factorization (default k-independent)
julia> q6j(j, j, j, j, j, j; mode=:generic)
√(z²⁴ Φ₂⁻⁸ Φ₃⁻⁴ Φ₄⁻⁴) × (-z⁻⁶ Φ₂² Φ₃ Φ₄  +  z⁻¹⁰ Φ₂² Φ₃ Φ₄ Φ₅)

# 4. Classical Ponzano-Regge Limit (q -> 1 or k -> ∞)
julia> q6j(j, j, j, j, j, j; mode=:classical)
0.16666666666666657
```
### 2. Exact Semiclassical Asymptotics (CycloMonomial Engine)

For macroscopic spins where standard floats fail, use the `:generic` engine to construct the level-independent algebraic invariant, and then defer specialization to arbitrary precision.
```julia
using QRacahSymbols

j = 550 # macroscopic spins

#1. Symbolic computation (very fast O(1) algebraic additions)
symb_gen = q6j(j, j, j, j, j, j; mode=:generic);

#2. Evaluate safely at physical level k = 50_000 using 512-bit precision
julia> evaluate_generic(symb_gen, Float64; k=50000, prec=512)
-3.627000456391523e-5

#3. Or analytically continue to a generic complex deformation parameter q
evaluate_generic(symb_gen, Complex{BigFloat}; q=0.5 + 0.1im)
```

## Documentation

For the complete API reference, interactive tutorials, and deep dives into the mathematical architecture, please see the [Official Documentation](https://sethkasante.github.io/QRacahSymbols.jl/).

## Citation

If you use `QRacahSymbols.jl` in your research, please cite our underlying methodological paper:

> Seth K. Asante, *"Exact Algebraic Evaluation of Quantum Recoupling Coefficients using Cyclotomic Representation"* (To appear soon, 2026). arXiv:XXXX.XXXXX.
