# Topological Symbols & Invariants

`QRecoupling.jl` provides the exact, high-performance building blocks required for computing spin networks, knot invariants, and 3D state-sum models.

This page demonstrates how to evaluate the core recoupling symbols, manipulate topological phases, and transition smoothly between discrete quantum topologies and continuous classical limits.

---

## The Wigner 6j-Symbol

The fundamental vertex in $\rm SU(2)_q$ recoupling theory is the quantum $6j$-symbol (or Racah-Wigner symbol). It represents the probability amplitude of a tetrahedral quantum geometry.

To evaluate a $6j$-symbol, simply provide the six half-integer spins $j_i$ and an evaluation target—such as the discrete Turaev-Viro level $k$.

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

### Generic Complex Deformations

If you are studying analytic continuation, hyperbolic geometries, or generic quantum groups, you can evaluate the geometry at any continuous complex deformation parameter $q$.

```julia
cval = q6j(1, 1, 1, 1, 1, 1, q=exp(0.5im))

```

---

## The Classical Limit (Ponzano-Regge)

As the level $k \to \infty$, the quantum deformation parameter $q \to 1$. This limit reduces the quantum spherical spacetime to a flat, classical Ponzano-Regge geometry.

Setting `q = 1` dynamically routes the package to evaluate standard classical angular momentum coefficients.

```julia
# Fast floating-point classical limit
cl_float = q6j(1, 1, 1, 1, 1, 1, q=1)

# Exact, zero-loss rational classical limit
cl_rational = q6j(1, 1, 1, 1, 1, 1, q=1, exact=true)
println(cl_rational)
# Output: Classical Result: 1//6

```

---

## Auxiliary Tensors & Symbols

When constructing large tensor networks, fusion categories, or knot invariants, you often need composite symbols that include specific geometric phases or dimension regularizations.

`QRecoupling.jl` provides direct APIs for the full suite of observables:

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

```

### Topological Phases (`QPhase`)

Tracking fractional framing anomalies and braiding eigenvalues is notorious for causing floating-point drift or algebraic type errors. `QRecoupling.jl` safely handles these via a deferred fractional architecture.

```julia
# R-Matrix (braiding eigenvalue for two crossings)
r_val = rmatrix(1, 1, 1, k=k)

# Exact Evaluation dynamically returns an isolated QPhase object
exact_phase = rmatrix(1/2, 1/2, 1, exact=true)
println(exact_phase)
# Output: q^(1//2)

```

---

## High-Performance quantum Wigner Symbols

To contract massive tensor networks (like an entire manifold triangulation), computation speed is critical.

By default, `QRecoupling.jl` evaluates symbols by first compiling a deferred cyclotomic graph (the DCR) to prevent overflow/NaN poisoning. However, if you are calling the same function millions of times in a tight loop and are confident in the scale of the boundaries, you can use the `eager=true` flag.

This bypasses the symbolic graph allocation entirely and pushes floating-point data straight into the optimized numerical solver.

```julia
# Bypasses algebraic graph allocation for maximum speed inside tight loops
fast_val = q6j(10, 10, 10, 10, 10, 10, k=5000, eager=true)
#Output: -0.002916327224229094

```

*(Note: `eager=true` is strictly available for discrete numerical evaluations where an integer level `k` is provided.)*

---
