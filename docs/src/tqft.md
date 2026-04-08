# TQFT & Quantum Gravity

`QRecoupling.jl` is purpose-built to provide the exact, high-performance topological kernels required for 3D Quantum Gravity (spin foam models), string-Net condensates, and state-sum invariants like the Turaev-Viro model.

This page demonstrates how to compute physical observables, manipulate tensor networks, and seamlessly transition between quantum and classical spacetime limits.

## The Quantum Recoupling Symbols
The core topological vertex in $SU(2)_k$ recoupling theory is the quantum $6j$-symbol (or Racah-Wigner symbol). It represents the probability amplitude of a tetrahedral quantum geometry. Other useful tensors are quantum $3j$ symbols. 

To evaluate a $6j$-symbol, simply provide the six half-integer spins $j_i$ and the topological level $k$.

```julia
using QRecoupling

# Evaluate at a root of unity (Turaev-Viro regime)
val = q6j(1, 1, 1, 1, 1, 1, k=10)
println(val) # 0.1547005383792515

# Evaluate the EXACT cyclotomic field representation
exact_val = q6j(1, 1, 1, 1, 1, 1, k=10, exact=true)
println(exact_val) # -2//3*ζ^6 + 4//3*ζ^2 - 1
```

### Complex Deformations
If you are studying analytic continuation, or generic quantum groups, you can evaluate the geometry at any continuous complex deformation parameter $q$
```julia
complex_val = q6j(1, 1, 1, 1, 1, 1, q=exp(0.5im))
```
---

## The Classical Limit (Ponzano-Regge)
As the level $k \to \infty$, the quantum deformation parameter $q \to 1$. This limit reduces the Turaev-Viro spherical spacetime to a flat Ponzano-Regge geometry.

Because standard cyclotomic evaluation fails at $q=1$ (due to $0/0$ division singularities), QRecoupling.jl mathematically detects q=1 and routes the calculation through a safe, exact geometric reduction path
```julia
# Fast floating-point classical limit
classical_float = q6j(1, 1, 1, 1, 1, 1, q=1)

# Exact zero-loss rational classical limit
classical_rational = q6j(1, 1, 1, 1, 1, 1, q=1, exact=true)
1//6
```

## Topological Tensors
When constructing massive tensor networks, you often need composite symbols that include specific phase conventions or quantum dimension regularizations.

The API provides direct access to these networks
```julia
k = 5

# Wigner 3j symbol (coupling of angular momenta)
val_3j = q3j(1, 1, 1, 0, 0, 0, k=k)

# F-Symbol (fusion tree crossing matrix)
f_val = fsymbol(1, 1, 1, 1, 1, 1, k=k)

# G-Symbol (tetrahedrally symmetric invariant)
g_val = gsymbol(1, 1, 1, 1, 1, 1, k=k)
```

## High-Performance Tensor Contractions
To contract large networks (like an $S^3$ manifold or a knot complement), performance is critical. You can bypass the internal algebraic graph construction and use the `eager=true` flag to push floating-point data straight into the Log-Sum-Exp numerical solver.
```julia
# Bypasses DCR graph allocation for maximum speed inside tight loops
fast_val = q6j(10, 10, 10, 10, 10, 10, k=50, eager=true)
```
Note: `eager=true` is only available when an integer level k is provided.