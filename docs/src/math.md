# Mathematical Framework

The **QRacahSymbols.jl** package implements the $q$-deformed Racah formula for the $SU(2)_k$ quantum group.

## The Quantum $6j$-Symbol
The core of the library is the evaluation of the quantum $6j$-symbol $\{j_1, j_2, j_3; j_4, j_5, j_6\}_q$, which represents the associativity of the fusion of three anyons in a topological quantum field theory.



### The Kirillov-Reshetikhin Formula
We evaluate the symbol using the following alternating sum:

$$\{j_1, j_2, j_3; j_4, j_5, j_6\}_q = \Delta(j_1, j_2, j_3) \Delta(j_1, j_5, j_6) \Delta(j_2, j_4, j_6) \Delta(j_3, j_4, j_5) \sum_z \frac{(-1)^z [z+1]_q!}{\prod_i [z - \alpha_i]_q! \prod_j [\beta_j - z]_q!}$$

Where:
* $[n]_q = \frac{q^{n/2} - q^{-n/2}}{q^{1/2} - q^{-1/2}}$ is the quantum integer.
* $\Delta(a, b, c)$ is the quantum triangle coefficient.
* The sum is stabilized in the `Numeric` engine using a **log-sum-exp** shift to prevent floating-point overflow.

## Cyclotomic Factorization
In the `Generic` mode, we exploit the fact that quantum integers $[n]_q$ are products of cyclotomic polynomials $\Phi_d(q)$. Any quantum $6j$-symbol can be represented as:

$$\text{Symbol} = \text{sign} \cdot z^{p} \cdot \prod_d \Phi_d(q)^{e_d}$$

where $z = q^{1/2}$. This allows for exact symbolic manipulation without ever choosing a specific value for $q$.