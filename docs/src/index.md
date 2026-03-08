# QRacahSymbols.jl Documentation
Welcome to the official documentation for **QRacahSymbols.jl**, a high-performance Julia library for evaluating quantum $6j$-symbols, $3j$-symbols, and topological category data for the $\text{SU(2)}_k$ quantum group.

## Overview
The calculation of quantum $6j$-symbols relies on the $q$-deformed Kirillov-Reshetikhin formula, an alternating hypergeometric sum over quantum factorials. For large spins and high levels, standard computational libraries inevitably fail due to catastrophic cancellation.

`QRacahSymbols.jl` resolves this by exposing a rigorous **Three-Tier Architecture**:

1. **`mode=:numeric`**: Employs an $O(1)$ table-lookup system fortified by a `log-sum-exp` shift, allowing for fast (nanoseconds) evaluation of tensor network invariants while remaining perfectly stable up to massive spin bounds.
2. **`mode=:exact`**: Leverages `Nemo.jl` to map the evaluation directly into the cyclotomic number field $\mathbb{Q}(\zeta_{2N})$, guaranteeing zero precision loss and returning exact topological invariants suitable for formal mathematical proofs.
3. **`mode=:generic`**: A highly optimized prime-factorization vector engine that returns the abstract polynomial decomposition of the symbol into cyclotomic roots ($\Phi_d(q)$).

## Installation
You can install the package directly from the Julia REPL. Press `]` to enter the Pkg prompt:
```julia
pkg> add QRacahSymbols
```

## Quick Start
The master `q6j` function dynamically dispatches to the most efficient computational engine based on whether you request a numerical level `k`, or specify a specific mode.


```julia
using QRacahSymbols
# 1. Fast Numerical Evaluation (Level k=20)
# Returns a Float64. Safe for spins up to j ≈ 450.
val_num = q6j(1, 1, 1, 1, 1, 1, 20)

# 2. Exact Algebraic Evaluation
# Returns a self-contained ExactValue holding a Nemo cyclotomic polynomial.
val_exact = q6j(1, 1, 1, 1, 1, 1, 20; mode=:exact)

# 3. Symbolic Cyclotomic Factorization
# Operates purely algebraically, bypassing all floating-point limitations.
val_symb = q6j(1, 1, 1, 1, 1, 1; mode=:generic)

# You can evaluate the symbolic result at a specific level k later:
evaluate_generic(val_symb, 20, Float64; prec=512)
```

## Next Steps

* Visit the [Mathematical Framework](@ref) page to understand how we bypass floating-point limits.
* Visit the [API Reference](@ref) page for details on function signatures, TQFT evaluators, and type stabilization strategies.