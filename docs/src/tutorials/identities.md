# Checking and proving identities

Numerical agreement tests conventions and implementations. An exact equality proves an identity for fixed labels, at one level or at every q. These are different kinds of evidence.

## Numerical orthogonality

```@example identities
using QRecoupling, Nemo, LinearAlgebra
F, e, f = fmatrix(1,1,1,1; k=6)
@assert transpose(F)*F ≈ I # hide
norm(transpose(F)*F - I)
```

This checks the whole recoupling matrix at one level to floating-point accuracy. It does not prove a symbolic identity for arbitrary labels or levels.

## Exact orthogonality at one level

For fixed outer spins equal to 1 and channel p = 1, orthogonality states `[3] * Σₓ [2x+1] {1 1 x; 1 1 1}² = 1`. At level 6 all channels x = 0, 1, 2 are allowed.

```@example identities
k = 6
lhs = sum(0:2) do x
    v = q6j(Exact(k), 1,1,x,1,1,1)
    qdim(Exact(k),1) * qdim(Exact(k),x) * v^2
end
@assert lhs == 1 # hide
lhs == 1
```

`==` is exact here. No floating comparison is involved, and no tolerance.

## The Biedenharn–Elliott relation

Here is a fixed-label pentagon check, assembled from the same arithmetic. With all nine external labels equal to 1, the phase of each summand is `(-1)^(9+x)`.

```@example identities
lhs_be = sum(0:2) do x
    v = q6j(Exact(k),1,1,x,1,1,1)
    (-1)^(9+x) * qdim(Exact(k),x) * v^3
end
rhs_be = q6j(Exact(k),1,1,1,1,1,1)^2
@assert lhs_be == rhs_be # hide
lhs_be == rhs_be
```

This verifies one admissible coefficient identity, not the full theorem for arbitrary labels. For general labels, derive the summation bounds, parity, and level truncation from the three left-hand symbols and check the right-hand symbols' admissibility. The test suite contains broader label samples.

## Two values at every q

`==` between two `Exact()` values decides equality at every q, without expanding either side. Here is a q-Vandermonde sum against one binomial.

```@example identities
vandermonde = qform(FactorialSum(0:3; prefactor = (3 => 1, 4 => 1), factors = ((1, 0, -2), (-1, 3, -1), (-1, 4, -1))); weight = (7, -24))
@assert vandermonde == qbinomial(Exact(), 7, 3) # hide
(vandermonde, qbinomial(Exact(), 7, 3), vandermonde == qbinomial(Exact(), 7, 3))
```

## Prove a generic identity in x

With Nemo loaded, the same diagonal orthogonality relation can be compared in expanded form before choosing a level:

```@example identities
lhs_x = let total = zero(XValue)
    for x in 0:2
        v = x_form(q6j(Symbolic(),1,1,x,1,1,1))
        weight = x_form(qdim(Symbolic(),1)) * x_form(qdim(Symbolic(),x))
        total = total + weight * v^2
    end
    total
end
certificate = prove_identity(lhs_x, one(XValue); kmax=30)
@assert certificate.verdict == :proved # hide
certificate
```

`prove_identity` compares normalized generic expressions by polynomial arithmetic. A proved identity is valid in the formal algebraic expression wherever the denominators are nonzero; complex numerical evaluation must retain consistent square-root branches. The returned `exceptional` list scans levels only through `kmax`. A canceled expression can have fewer visible denominator zeros than the original summands, so this list is not a substitute for checking their admissibility or singularities.

## When a zero test is a certificate

`iszero` and `==` give exact answers, independent of any numerical tolerance. A nonzero residue at a word prime proves a value nonzero at once. A zero is proved by vanishing at enough primes for a bound on the value to leave no other possibility. Equal values cost more than different ones, and the cost grows with the size of the formula. `iszero_at` performs exact confirmation of unresolved modular candidates itself. The `:cancels` entries of `level_spectrum` are screening results unless you pass `prove = true`.
