
# Generic $q$-Series & Custom DCRs

While `QRecoupling.jl` is highly optimized for TQFT invariants and angular momentum algebra, its underlying engine is a completely generic, abstract mathematical compiler for basic hypergeometric series (${}_r\Phi_s$).

To study partition identities, quantum group observables built from $q$-hypergeometric polynomials or general finite $q$-series, you can use the `qseries` API to safely and exactly evaluate your own custom sequences. 

By compiling your mathematical sequence into a **Deferred Cyclotomic Representation (DCR)**, the package eliminates intermediate expression swell, catastrophic cancellation, and numerical overflow.

---

## The `qseries` Frontend

`QRecoupling.jl` provides a suite of symbolic primitives that are user friendly. These primitives natively support scalar multiplication, division, and alternating signs:
* `qint(n)` $\to [n]_q$
* `qfact(n)` $\to [n]_q!$
* `qbinomial(n, k)` $\to \begin{bmatrix} n \\ k \end{bmatrix}_q$

To compile a sequence, simply pass a `UnitRange` and a `do` block defining the $z$-th term to the `qseries` function. 

### Example 1: A Simple Sequence
Let's build a highly numerically unstable sequence—the sum of $q$-factorials:
$$S(q) = \sum_{z=1}^{10} [z]_q!$$

```julia
using QRecoupling

ex1_series = qseries(1:10) do z
    return qfact(z)
end

println(ex1_series)
```

### Example 2: Alternating Series and Algebra
Because `CyclotomicMonomial` objects natively hook into Julia's `Base` arithmetic, you can divide primitives and apply standard mathematical parity flags like `(-1)^z` directly.

$$Z(q) = \sum_{z=3}^{10} (-1)^z \frac{[z]_q!}{[z-3]_q!}$$

```julia
my_series = qseries(3:10) do z
    return (-1)^z * (qfact(z) / qfact(z - 3))
end
```

---

## The DCR Engine

When you call `qseries`, the package does not evaluate the sum. Instead, it algebraically compiles your sequence into a purely symbolic graph called a **DCR** (Deferred Cyclotomic Representation). 

A finite $q$-hypergeometric series takes the form $Z(q) = \text{Prefactor} \times \sum T(z)$. The DCR avoids generating massive polynomials by isolating three structural components:
1.  **The Root Prefactor:** A global multiplier pulled outside the sum.
2.  **The Base Term $T(z_{\text{min}})$:** The algebraic value of the very first term.
3.  **The Update Ratios $R_z$:** The multiplicative jump between adjacent terms, such that $T(z+1) = T(z) \times R_z$.

By strictly tracking the algebraic ratios $R_z$, `QRecoupling.jl` can evaluate sums with massive quantum factorials by utilizing the **Log-Sum-Exp** algorithm. This ensures numerical stability even when individual terms exceed standard `Float64` limits.

---

## Universal Evaluation (`qeval`)

Once your sequence is compiled into a DCR object, you can project it into any target field using the universal `qeval` API. The target parameters automatically route the DCR to the optimal numerical or exact solver without needing to recompile the algebra.

```julia
# Classical limit (q -> 1)
val_limit = qeval(my_series, q=1.0)

# Discrete root-of-unity (level k)
val_k = qeval(my_series, k=5)

# Exact algebraic (Nemo)
# Evaluate exactly in the cyclotomic field ℚ(ζ_2h)
val_exact = qeval(my_series, k=5, exact=true)

# Complex Analytic
# Evaluate across the complex plane
val_complex = qeval(my_series, q=0.5 + 0.2im)
```

---

## High-Performance Expert API: `build_dcr!`

**Note:** *For most use cases, `qseries` is the recommended API.*

In `qseries`, the engine automatically performs algebraic division to find the update ratios $R_z = T(z+1) / T(z)$. However, if you are computing millions of state sums, you may want to avoid the allocation overhead of these intermediate divisions.

For absolute, zero-allocation performance, developers can bypass the automatic algebra and use the low-level `build_dcr!` engine. This requires manually passing the mathematically simplified update ratios via closures over a pre-allocated `CycloBuffer`.

```julia
# Pre-allocate an integer buffer
buf = QRecoupling.CycloBuffer(20)

# Manually compile S(q) = \sum_{z=1}^{10} [z]_q!
ex2_series = build_dcr!(buf,
    b -> nothing,                       # Prefactor logic
    (b, z) -> add_qfact!(b, z),         # Base term at z_min
    (b, z) -> add_qint!(b, z + 1),      # Manually simplified ratio
    1, 10
)
```
*For examples of this zero-allocation pattern, see the internal implementation of `q6j` in the source code.*
