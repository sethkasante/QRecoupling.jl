# Generic $q$-Series & Custom DCRs

While `QRecoupling.jl` is highly optimized for TQFT physics, its underlying engine is a completely generic, abstract mathematical compiler for basic hypergeometric series ($_r\Phi_s$).

To study partition identities, Rogers-Ramanujan type sums, or quantum group observables built from basic hypergeometric series, you can use the **Deferred Cyclotomic Representation (DCR)** to safely and exactly evaluate your own custom $q$-series.

## What is a DCR?
A finite $q$-hypergeometric series takes the form:
$$Z(q) = \text{Prefactor} \times \sum_{z = z_{\text{min}}}^{z_{\text{max}}} T(z)$$

The DCR engine avoids generating massive, unstable polynomials by isolating three distinct components:
1. **The Root Prefactor:** A global multiplier outside the sum.
2. **The Base Term $T(z_{\text{min}})$:** The value of the very first term in the sum.
3. **The Update Ratios $R_z$:** The multiplicative jump between adjacent terms, such that $T(z+1) = T(z) \times R_z$.

By providing rules for these three components, the `build_dcr` compiler generates a purely algebraic graph that is completely immune to intermediate expression swell.

## Building a Custom Sequence

Let's build a highly unstable sequence: the sum of $q$-factorials:
$$S(q) = \sum_{z=1}^{10} [z]_q!$$

The update ratio between $T(z)$ and $T(z+1)$ is simply $[z+1]_q$. Here is how we compile this into a DCR object:

```julia
using QRecoupling

# Compile the sequence
custom_series = build_dcr(
    # 1. Prefactor (Nothing needed here)
    b -> nothing,                       
    
    # 2. Base Term (The first term when z = 1 is just [1]!)
    (b, z) -> add_qfact!(b, z),         
    
    # 3. Update Ratio (To get to the next term, multiply by [z+1])
    (b, z) -> add_qint!(b, z + 1),      
    
    # 4. Summation Bounds
    1, 10                               
)

println(custom_series)
```
---

## Universal Projection
Once the sequence is compiled into a DCR object, you can evaluate it in any field using the master `project_dcr` API. The target parameters automatically route the DCR to the correct numerical or exact solver.
```julia
# 1. Evaluate exactly at the classical limit (q -> 1)
val_limit = project_dcr(custom_series, q=1, exact=true)

# 2. Evaluate at a root of unity (k = 5)
val_k = project_dcr(custom_series, k=5)

# 3. Evaluate analytically at a complex deformation
val_complex = project_dcr(custom_series, q=0.5 + 0.2im)
```