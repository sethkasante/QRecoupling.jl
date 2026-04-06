# Theory and Architecture

Evaluating quantum $6j$-symbols at high spins is a notoriously hostile computational problem. The standard Racah-Wigner hypergeometric sum requires computing massive $q$-factorials, multiplying them, dividing them, and extracting algebraic square roots. 

Standard floating-point implementations inevitably suffer from catastrophic `NaN` or `Inf` overflows. Conversely, executing the sum in a Computer Algebra System (CAS) using dense cyclotomic polynomials triggers severe memory bloat and $\mathcal{O}(N^3)$ polynomial division overhead.

`QRecoupling.jl` bypasses both bottlenecks by decoupling the **topological graph construction** from the **mathematical evaluation**.

---

## The Core Innovation: Deferred Factorization

Instead of eagerly evaluating $q$-integers like $[n]_q$, the package intercepts the mathematics at the prime factorization level. Every component of the Racah sum is factored into irreducible cyclotomic polynomials $\Phi_d(q)$.

This data is stored in a highly optimized sparse array called a `CycloMonomial`, which tracks the integer exponents of each polynomial. 

The overall $6j$-symbol is deferred into a `CycloResult` struct representing the hypergeometric series:
$$\begin{Bmatrix} j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6 \end{Bmatrix} = \left( \Delta_{\text{root}} \sqrt{\Delta_{\text{rad}}} \right) M_0 \left[ 1 + R_1 + R_1 R_2 + \dots \right]$$
By keeping the series in this abstract ratio format, we completely eliminate algebraic division during the builder phase. 

---

## 1. The Exact Zero-Division Engine

When evaluated in a rigorous $\text{SU(2)}_k$ cyclotomic field using `Nemo.jl`, `QRecoupling.jl` employs a **Zero-Division Architecture**. Use as (`mode=:exact`). 

Because the hypergeometric structure is factored into prime cyclotomic bases *before* evaluation, the package executes two critical optimizations:
1. **Symbolic Square Roots:** Algebraic square roots are the heaviest operation in any CAS. By dividing the integer exponents of the `CycloMonomial` by 2 (via `divrem(e, 2)`), perfect squares are pulled out of the radical *symbolically* in $\mathcal{O}(1)$ time. The exact engine only ever evaluates the strictly square-free radical ($\Delta_{\text{rad}}$).
2. **Division-Free Sequences:** By precomputing the exact algebraic inverses of the cyclotomic bases at the requested root of unity, the hypergeometric sum $1 + R_1 + R_1 R_2 \dots$ is evaluated using strictly multiplication and addition. 

This enables exact topological proofs (like the Biedenharn-Elliott Pentagon Identity) at spins where traditional CAS approaches would exhaust system RAM.

---

## 2. The Log-Sum-Exp Numeric Engine

For high-speed floating-point simulations (e.g., Turaev-Viro invariant state sums over triangulations), the package relies on a specialized **Log-Sum-Exp (LSE) hot loop**. Use as (`mode=:numeric`).

Rather than using the `CycloResult` pipeline, the `:numeric` mode jumps straight to bare-metal CPU registers. 
1. The engine pre-caches a dense look-up table of logarithmic $q$-factorials: $\ln([n]_q!)$.
2. The boundaries of the Racah sum are analyzed to find the maximum logarithmic term ($M_{\text{max}}$).
3. The alternating sum is shifted by this maximum before exponentiation: $\sum \exp(\ln(M_z) - M_{\text{max}})$.

This guarantees absolute immunity against floating-point underflow/overflow, evaluating massive tetrahedral networks in nanoseconds.

---

## 3. The Classical GMP Engine

In the Ponzano-Regge classical limit ($q \to 1$), the $6j$-symbol reduces to purely rational arithmetic. However, standard Julia `BigInt` arithmetic creates massive Garbage Collection (GC) overhead when millions of fractions are allocated during a sum. Use as (`mode=:classical_exact`)

The `:classical_exact` engine utilizes direct `ccall` bindings to the underlying **GNU Multiple Precision (GMP)** C-library (`Base.GMP.MPZ`). 
* The hypergeometric ratios are tracked via dynamic integer prime sieves.
* The summation numerator is accumulated strictly in-place, bypassing Julia's GC entirely.

This makes the classical engine competitive with pure **C/C++** Wigner symbol libraries while maintaining mathematically exact `Rational{BigInt}` outputs.