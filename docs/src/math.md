# Theory & Architecture

The `QRacahSymbols.jl` package provides a high-performance framework for the evaluation of $q$-deformed recoupling coefficients in the $\text{SU(2)}_k$ modular tensor category.

## Quantum recoupling coefficients
The fundamental object of the library is the quantum $\{6j\}$-symbol: 
$$\begin{Bmatrix} j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6\end{Bmatrix}_q.$$

In the context of Topological Quantum Field Theory (TQFT), this symbol represents the associativity of fusion (the $F$-move) for three anyons of spins $j_1, j_2, j_3$ into a total spin $j_6$.

These symbols are defined by the Quantum Racah Formula, an alternating summation of massive products of quantum factorials. The building blocks are:

* **Quantum integers:** defined as $[n]_q = \frac{q^{n/2} - q^{-n/2}}{q^{1/2} - q^{-1/2}}$. At the level $k$, we typically specialize to the root of unity $q = e^{i \frac{2\pi}{k+2}}$. 
* **Quantum factorials:** given by $[n]_q! = \prod_{k=1}^n [k]_q$,  with $[0]_q!=1$. 

### The Numerical Precision Crisis
In standard floating-point arithmetic (`Float64`), the Racah sum encounters a **catastrophic cancellation** barrier. As the spins $j$ increase, the individual summands grow exponentially, while the final physical invariant oscillates and scales as $O(j^{-1/2})$.

To resolve a value of $10^{-10}$ from summands of magnitude $10^{50}$, a bit-depth is required that scales linearly with the spin. Without a specialized architecture, numerical "noise" completely obliterates the signal for spins $j > 50$.

## The CycloMonomial Architecture
To bypass these limits, `QRacahSymbols.jl` implements a `CycloMonomial` engine. Rather than projecting quantum factorials onto the complex plane prematurely, we "lift" the calculation into the cyclotomic ring $\mathbb{Z}[\zeta]$.

### Algebraic Factorization
We exploit the property that every quantum integer $[n]_q$ is a product of cyclotomic polynomials $\Phi_d(q)$. Specifically, for $z = q^{1/2}$:

$$[n]_q = \prod_{d|n, d>1} \Phi_d(z).$$

Consequently, any quantum factorial $[n]_q!$ can be uniquely represented as a symbolic monomial:
$$[n]_q! = s \cdot z^{p} \cdot \prod_d \Phi_d(z)^{e_d}$$

where $s \in \{-1, 1\}$ is a sign and $e_d$ are integer exponents.

### The Generic Engine (`:generic`)
In the `:generic` mode, the package represents the Racah sum as a collection of **Prime-Power Exponent Arrays**.
    1. **Symbolic Cancellation:** Massive factorial products in the numerator and denominator are cancelled by simply subtracting their exponent vectors. This is an $O(1)$ operation that precedes any numerical evaluation. 
    2. **Delayed Specialization:** The "subtraction" of large terms happens exactly in the symbolic space. Numerical evaluation only occurs at the very last step, using high-precision stabilization.
    3. **Structural Zeros:** This architecture natively recognizes the topological truncation of the fusion category. If a symbol is zero due to the level $k$ constraint, the engine encounters a symbolic factor of $\Phi_{k+2}(q)$, which is exactly zero at the root of unity.


## Topological Category Data

Beyond the $\{6j\}_q$-symbols, the package provides the complete set of structural constants required for Turaev-Viro state-sum models and anyonic braiding.

| Function | Physical Meaning | Output Type |
| --- | --- | --- |
| `q6j`, `q3j` | Quantum Wigner $6j$- and $3j$-symbols | `GenericResult` |
| `qdim` | Quantum Dimensions $[2j+1]_q$ | `CycloMonomial` |
| `fsymbol` | Unitary $F$-matrix (fusion basis transformation) | `Numeric/Generic` |
| `rmatrix` | $R$-matrix (braiding phases/statistics) | `Numeric/Generic` |
| `evaluate_generic` | Analytic continuation to arbitrary $q \in \mathbb{C}$ | `Complex/Real` |
