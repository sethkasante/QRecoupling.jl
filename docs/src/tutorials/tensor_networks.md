# Tensor networks and modular data

## Change a fusion basis

For four labels a, b, c, d, `fmatrix` changes between the intermediate channels of
`((a b)ₑ c)ᵈ` and `(a (b c)𝒇)ᵈ`. It returns the matrix and both channel lists.

```@example networks
using QRecoupling, LinearAlgebra
F, e, f = fmatrix(1,1,1,1; k=6)
@assert size(F) == (length(e), length(f)) # hide
@assert transpose(F)*F ≈ I # hide
@assert all(isapprox(F[i,j], fsymbol(1,1,e[i],1,1,f[j]; k=6)) for i in eachindex(e), j in eachindex(f)) # hide
(F, e, f)
```

Columns express f-basis states in the e basis, so `state_e = F * state_f`. Keep channel
labels with each tensor leg; array positions alone do not specify the fusion channel.

```@example networks
state_f = zeros(length(f)); state_f[1] = 1
state_e = F * state_f
@assert transpose(F)*state_e ≈ state_f # hide
state_e
```

For complex q, the algebraic relation is bilinear:

```@example networks
Fc, ec, fc = fmatrix(1,1,1,1; q=0.8+0.2im)
@assert transpose(Fc)*Fc ≈ I # hide
norm(transpose(Fc)*Fc - I)
```

Use `transpose`, not the adjoint, for this identity. Generic complex-q matrices need not
preserve a Hilbert-space norm. See [Accuracy and performance](../performance.md) before
requesting higher matrix element types.

## Braiding in the fusion basis

```@example networks
B, incoming, outgoing = bmatrix(1,1,1,1; k=6)
@assert B' * B ≈ I # hide
(size(B), incoming, outgoing)
```

`bmatrix` assembles the basis changes and channel R phases. The ordered channel lists
matter when the exchanged labels differ. The scalar `rmatrix` is a single-channel
braiding eigenvalue, not a full many-body operator.

## Modular data

```@example networks
S, labels = smatrix(3)
T, labels_T = tmatrix(3)
@assert labels == labels_T # hide
@assert S*S ≈ I # hide
@assert (S*T)^3 ≈ S*S # hide
(S, T, QRecoupling.total_qdim(3))
```

`tmatrix` includes the central-charge anomaly by default; `anomaly=false` retains only
the twists. `QRecoupling.level_labels`, `twist`, `QRecoupling.central_charge`, `QRecoupling.gauss_sum`, `QRecoupling.monodromy`, and
`QRecoupling.verlinde` expose the related data. `QRecoupling.verlinde` computes and rounds a numerical fusion
coefficient with a consistency check; it is not an exact symbolic summation routine.

## From local data to an application

A tensor-network application can enumerate allowed channels, cache the required F and R
data at a fixed level, attach labels to sparse blocks, and contract using its chosen tensor
library. QRecoupling supplies the local recoupling data. Network construction, contractions,
truncation, optimization, and observable definitions belong to the surrounding application.

Before scaling up, compare every matrix convention against scalar entries and check a
small identity or known network evaluation. These checks catch channel-order and phase
errors that high precision cannot fix.
