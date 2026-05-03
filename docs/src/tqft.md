# TQFT & Quantum Gravity

`QRecoupling.jl` is purpose-built to provide the exact, high-performance topological kernels required for 3D Quantum Gravity (spin foam models), string-net condensates, and state-sum invariants like the Turaev-Viro model.

This page demonstrates how to compute physical observables, manipulate tensor networks, and seamlessly transition between quantum and classical spacetime limits.

---

## The Quantum Recoupling Symbols
The core topological vertex in $SU(2)_q$ recoupling theory is the quantum $6j$-symbol (or Racah-Wigner symbol). It represents the probability amplitude of a tetrahedral quantum geometry. 

To evaluate a $6j$-symbol, simply provide the six half-integer spins $j_i$ and the evaluation target (such as the topological level $k$).

```julia
using QRecoupling

# Evaluate at a discrete root of unity (Turaev-Viro regime)
val = q6j(1, 1, 1, 1, 1, 1, k=10)
println(val) 
# Output: 0.1547005383792515

# Evaluate exactly in the cyclotomic field ℚ(ζ)
exact_val = q6j(1, 1, 1, 1, 1, 1, k=10, exact=true)
println(exact_val) 
# Output: Exact Algebraic Result in ℚ(ζ₂₄): Value: -2//3*ζ^6 + 4//3*ζ^2 - 1
```

### Complex Deformations
If you are studying analytic continuation or generic quantum groups, you can evaluate the geometry at any continuous complex deformation parameter $q$.
```julia
cval = q6j(1, 1, 1, 1, 1, 1, q=exp(0.5im))
```

---

## The Classical Limit (Ponzano-Regge)
As the level $k \to \infty$, the quantum deformation parameter $q \to 1$. This limit reduces the Turaev-Viro spherical spacetime to a flat Ponzano-Regge geometry.

Setting `q = 1` dynamically routes the package to evaluate standard classical angular momentum coefficients.

```julia
# Fast floating-point classical limit
cl_float = q6j(1, 1, 1, 1, 1, 1, q=1)

# Exact, zero-loss rational classical limit
cl_rational = q6j(1, 1, 1, 1, 1, 1, q=1, exact=true)
println(cl_rational)
# Output: 1//6
```

---

## Topological Tensors
When constructing large tensor networks or modular tensor categories, you often need composite symbols that include specific geometric phases or quantum dimension regularizations.

`QRecoupling.jl` provides direct APIs for the full suite of TQFT observables:

```julia
k = 5

# 1. Quantum Dimensions [2j+1]_q
dim = qdim(1, k=k)

# 2. Wigner 3j Symbol (coupling of angular momenta)
val_3j = q3j(1, 1, 1, 0, 0, 0, k=k)

# 3. F-Symbol (fusion tree crossing matrix)
f_val = fsymbol(1, 1, 1, 1, 1, 1, k=k)

# 4. G-Symbol (tetrahedrally symmetric invariant for state sums)
g_val = gsymbol(1, 1, 1, 1, 1, 1, k=k)

# 5. R-Matrix (braiding eigenvalue)
r_val = rmatrix(1, 1, 1, k=k)
```
*Note: Requesting `exact=true` on phase factors like `rmatrix` automatically scales the cyclotomic field extension to safely encapsulate fractional powers of $q$.*

---

## High-Performance Tensor Contractions
To contract massive tensor networks (like an $S^3$ manifold or a dense knot complement), computation speed is critical. 

By default, `QRecoupling.jl` evaluates symbols by first compiling an exact symbolic graph (the DCR) to prevent overflow. However, if you are calling the same function millions of times in a tight loop, you can use the `eager=true` flag. This bypasses the symbolic graph allocation entirely and pushes floating-point data straight into the optimized Log-Sum-Exp numerical solver.

```julia
# Bypasses algebraic allocation for maximum speed inside tight loops
fast_val = q6j(10, 10, 10, 10, 10, 10, k=50, eager=true)
```
*(Note: `eager=true` is strictly available for discrete numerical evaluations where an integer level `k` is provided.)*