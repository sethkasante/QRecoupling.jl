# What can you build with QRecoupling?

The package is most useful when a calculation needs many related recoupling coefficients, or needs to move between numerical experiments and exact formulas without rewriting the underlying symbol definitions.

## Fusion-based tensor networks

Use `fmatrix`, `rmatrix`, and `bmatrix` as local changes of basis and braiding operations. A practical first study is to measure how allowed channel counts and local tensor sizes grow with the level, then compare contractions with and without fusion-sector sparsity. The [tensor-network tutorial](tutorials/tensor_networks.md) provides tested local building blocks. A network optimizer or contraction engine is not included.

## Spin networks and state-sum models

Use q6j coefficients and quantum dimensions for tetrahedral and edge weights, with your model's phase and normalization conventions. Compare a small network at several levels with its classical limit, or check a local recoupling move before implementing a full state sum. Summing over triangulations, boundary conditions, and observables is application work.

## Exact identities and special levels

Explore when a coefficient vanishes, reduces to a rational number, or admits a short radical expression. `Symbolic()` and `x_form` expose generic formulas; `Exact(k)` specializes to the real algebraic basis of a level. Generic polynomial identities can establish a relation wherever its denominators and branch assumptions permit specialization. The [identity tutorial](tutorials/identities.md) separates such proofs from numerical checks.

## Deformations and asymptotics

Sweep k toward the classical limit or vary a real/complex q while keeping labels fixed. Batches reuse target tables; scalar BigFloat calculations provide higher-precision checks for difficult points. Near poles and roots of unity, distinguish sensitivity to q from error in the evaluation. Asymptotic fits and geometric interpretations remain the researcher's task.

## Finite series and new symbols

Encode factorial-product sums with `FactorialSum` and test the same expression classically, at a level, and at generic q. This supports terminating hypergeometric examples with integer factorial parameters and offers a starting point for new recoupling formulas. Arbitrary parameterized q-Pochhammer factors and infinite basic hypergeometric series are outside the current interface.

## Choosing a first project

| Question | Starting point | Useful validation |
|:--|:--|:--|
| How do fusion spaces grow with k? | `fmatrix_labels`, `QRecoupling.level_labels` | Check channel admissibility |
| Are my tensor phases consistent? | F/R matrices and scalar symbols | Orthogonality and local identities |
| Does a relation hold beyond sampled levels? | `x_form`, `prove_identity` | Exact polynomial residual and exceptional levels |
| Is a small numerical coefficient really zero? | Exact level or generic formula | Exact cancellation, not only a tolerance |
| Does a custom series recover a classical identity? | `FactorialSum`, `Exact()` | Independent factorial/binomial formula |

Performance depends on the labels, target, precision, and workload. Benchmark the intended application rather than extrapolating a scalar timing to a full tensor contraction.
