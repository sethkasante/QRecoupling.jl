```@meta
CurrentModule = QRecoupling
```
# QRecoupling.jl
Welcome to the official documentation for **QRecoupling.jl**, a high-performance Julia library for evaluating quantum $6j$-symbols, $3j$-symbols, and topological category data for the $\text{SU(2)}_k$ quantum group.

`QRecoupling.jl` is designed for researchers in Quantum Topology, Quantum Gravity, and fusion category theory. It provides highly optimized, numerically stable, and mathematically exact evaluations of quantum $\text{SU(2)}_k$ recoupling coefficients (such as $3j$- and $6j$-symbols) at arbitrary roots of unity.

By representing quantum prime factorizations as deferred sparse cyclotomic arrays, this package bypasses traditional algebraic bottlenecks and floating-point overflows, bridging the gap between very-fast numerical simulations and rigorous Computer Algebra System (CAS) proofs.

---

## Installation
You can install the package directly from the Julia REPL. Press `]` to enter the Pkg prompt, and run:
```julia
pkg> add QRecoupling
```

---

## The Multi-Paradigm Architecture
Evaluating quantum $6j$-symbols at high spins using standard floating-point evaluation suffers from catastrophic NaN/Inf overflows due to massive factorials, while traditional exact evaluation in cyclotomic fields scales terribly due to the $\mathcal{O}(N^3)$ complexity of algebraic division. 

`QRecoupling.jl` solves this by decoupling the **algebraic construction** from the **mathematical evaluation**.

Using the mode keyword, users can route the computation to specialized engines:
1. **Numeric Mode** (`mode=:numeric`) The absolute speed limit. Uses pre-cached dense logarithmic arrays and a highly specialized `Log-Sum-Exp` (LSE) hot loop to compute single-shot floating-point values. It is mathematically immune to overflow.
2. **Exact Mode** (`mode=:exact`) Maps the symbol into rigorous `Nemo.jl` cyclotomic fields. By extracting algebraic perfect squares symbolically before summation, it provides $\mathcal{O}(1)$ radical extraction and impenetrable, division-free loops.
3. **Classical Exact Mode** (`mode=:classical_exact`) Evaluates the un-deformed Ponzano-Regge limit ($q \to 1$). Bypasses Julia's standard Garbage Collector entirely by executing GMP integer math strictly in-place, yielding zero-allocation exact rational numbers.
4. **Cyclotomic Mode** (`mode=:cyclo`) Builds the unevaluated, sparse algebraic graph of the hypergeometric series. Ideal for structural analysis, zero-detection, and advanced parameter sweeps.

---

## Quick Start
The master `q6j` function dynamically dispatches to the most efficient computational engine based on whether you request a numerical level `k`, or specify a specific `mode`.
```julia
using QRecoupling
# Fast numerical evaluation (level k=20)
julia> q6j(1,1,1,1,1,1,20)
0.1640608932500723

# Exact algebraic evaluation: Returns an exact value holding a `Nemo` cyclotomic polynomial.
julia> val_exact = q6j(1,1,1,1,1,1,20; mode=:exact)
Exact SU(2)₂₀ Symbol:
  Value: 2*ζ^18 - 3*ζ^16 + ζ^14 + 2*ζ^12 - 2*ζ^10 - ζ^8 + 3*ζ^6 - 2*ζ^4 + 1

# Cyclotomic representation independent of k (in hypergeometric form)
julia> q6j(1,1,1,1,1,1; mode=:cyclo)
CycloResult (Hypergeometric Ratio Form)
  ├─ Max Φ_d(q) required : d = 5
  ├─ Overall Prefactor : -z⁶ Φ₂⁻² Φ₃⁻¹ Φ₄⁻¹
  └─ Ratios (R)     : 1 terms
       └─ R₁ : -z⁻⁴ Φ₅
```

---

## Topological symbols
Beyond basic $6j$-symbols, the package provides a first-class API for constructing the composite tensors useful for Turaev-Viro state sums and categorical proofs.
* **Quantum Dimensions**: `qdim(j, k)`
* **R-Matrices (Braiding)**: `rmatrix(j1, j2, j3, k)`
* **F-Symbols (Unitary Fusion)**: `fsymbol(j1, j2, j3, j4, j5, j6, k)`
* **G-Symbols (Symmetric Tetrahedral Weights)**: `gsymbol(j1, j2, j3, j4, j5, j6, k)`

---

## Next Steps

* Visit the [Theory & Architecture](@ref) page to understand how we bypass floating-point limits.
* Visit the [API Reference](@ref) page for details on function signatures, TQFT evaluators, and type stabilization strategies.