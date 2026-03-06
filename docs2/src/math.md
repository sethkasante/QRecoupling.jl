# Mathematical Framework



## Tetrahedral Symmetry
The quantum $6j$ symbol $\{j_1, j_2, j_3; j_4, j_5, j_6\}_k$ can be mapped to the edges of a tetrahedron. It exhibits a 24-fold $S_4$ geometric symmetry. `QRacahSymbols.jl` uses a zero-allocation canonicalization algorithm to map any input configuration to its lexicographically maximal symmetry partner, allowing aggressive LRU caching.

## The Racah Sum
The internal engines evaluate the Kirillov-Reshetikhin formula:
$$\sum_{z = \max(\alpha_i)}^{\min(\beta_j)} \frac{(-1)^z [z+1]_q!}{\prod_i [z - \alpha_i]_q! \prod_j [\beta_j - z]_q!}$$

To prevent floating-point overflow during numeric evaluations, the alternating sum is stabilized using a mathematically exact `log-sum-exp` shift.