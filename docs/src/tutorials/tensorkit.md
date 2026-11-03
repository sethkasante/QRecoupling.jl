# TensorKit sectors

Loading [TensorKit](https://github.com/QuantumKitHub/TensorKit.jl) makes SU(2) at a level, or at real q, available as a tensor symmetry. `SU2k(k)` and `SU2q(q)` return sector types whose F- and R-symbols come from [`fsymbol`](@ref) and [`rmatrix`](@ref). The extension needs only TensorKitSectors, so Nemo is not involved.

```julia
using TensorKit, QRecoupling

I = SU2k(3)                                  # SU(2) at level 3, q = e^{iπ/5}
a = I(1)                                     # spin 1, also the τ of the Fibonacci category
collect(a ⊗ a)                               # [I(0), I(1)]: fusion is truncated by the level
dim(a)                                       # 1.618…, the quantum dimension
Fsymbol(a, a, a, a, a, a)                    # −0.618…
Rsymbol(a, a, I(0))                          # a phase

V = Vect[I](0 => 2, 1//2 => 2, 1 => 1)
t = randn(ComplexF64, V ⊗ V ← V ⊗ V)
norm(repartition(t, 1, 3)) ≈ norm(t)         # F-moves are unitary
braid(t, ((2, 1), (3, 4)), (1, 2, 3, 4))     # braiding, with levels for the strands
U, S, Vh = svd_compact(t)
```

The conventions are TensorKit's, F^{abc}_{d,ef} = (−1)^{a+b+c+d} √([2e+1][2f+1]) {a b e; c d f}, which is `fsymbol(a, b, e, c, d, f)`. The field `a.J` of a label is `2j`.

## What each type supports

| | `SU2k(k)` | `SU2q(q)`, q > 0 |
|---|---|---|
| Spins | 0, 1/2, …, k/2 | every multiple of 1/2 |
| Fusion | truncated by the level | as for SU(2) |
| F-symbols | real and orthogonal | real and orthogonal |
| R-symbols | complex phases | none: they are not unitary, which TensorKit assumes |
| Tensor operations | planar and braided | planar only |

The F-symbols are `Float64` in both, so planar tensors can stay real. `SU2q(1.0)` agrees with `SU2Irrep`. Fusion tensors are not provided, since the quantum dimension of a sector is not the size of a matrix block, so `convert(Array, t)` is unavailable, as for TensorKit's other anyons.

TensorKit warns that tensors with real data might be incompatible with `SU2k(k)`. That matters only for a tensor that is braided, since the R-symbols are complex. Nemo and TensorKit both export `dim`, so with both loaded write `TensorKit.dim`, or load Nemo with `import Nemo`.

## Large levels

An F-symbol takes 0.1 to 0.6 µs at levels up to 20000, with relative error below 5 × 10⁻¹⁴ on samples against high-precision values. Above level 32 a graded space stores its multiplicities in a dictionary instead of in its type, so large levels do not slow compilation.

## Checks

`test/tensorkit_environment.jl` runs TensorKit's own sector test suite on levels 1 to 40 and four values of q, every pentagon and hexagon equation at levels up to 6, sampled ones up to k = 5000, and tensor operations at a level and at real q. The integer spins of level 3 reproduce the chain spectra of TensorKit's `FibonacciAnyon`, and `SU2q(1.0)` the F-symbols of `SU2Irrep`. Chains of spins 1/2 at a level and at real q have the spectra of the U_q-invariant XXZ chain, and the ground-state energy densities of infinite chains found with MPSKit's VUMPS lie within 4 × 10⁻⁶ above the Bethe-ansatz values for levels 2 to 5.
