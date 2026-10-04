# Checking and proving identities

Numerical agreement tests conventions and implementations. An exact residual can prove an identity for fixed labels. A generic rational-function identity can prove it before choosing a level. These are different kinds of evidence.

## Numerical orthogonality

```@example identities
using QRecoupling, LinearAlgebra
F, e, f = fmatrix(1,1,1,1; k=6)
@assert transpose(F)*F ≈ I # hide
norm(transpose(F)*F - I)
```

This checks the whole recoupling matrix at one level to floating-point accuracy. It does not prove a symbolic identity for arbitrary labels or levels.

## Exact orthogonality at one level

For fixed outer spins equal to 1 and channel p = 1, orthogonality states `[3] * Σₓ [2x+1] {1 1 x; 1 1 1}² = 1`. At level 6 all channels x = 0, 1, 2 are allowed.

```@example identities
k = 6
lhs = let total = ExactXSum(k)
    for x in 0:2
        v = q6j(Exact(k), 1,1,x,1,1,1)
        total = total + qdim(Exact(k),1) * qdim(Exact(k),x) * v^2
    end
    total
end
one_at_level = qdim(Exact(k),0)
residual = lhs - one_at_level
@assert isempty(residual) # hide
isempty(residual)
```

An **empty exact residual** means every stored coefficient canceled exactly. No floating comparison is needed for this example.

## The Biedenharn–Elliott relation

Here is a fixed-label pentagon check, assembled from the same public arithmetic. With all nine external labels equal to 1, the phase of each summand is `(-1)^(9+x)`.

```@example identities
lhs_be = let total = ExactXSum(k)
    for x in 0:2
        v = q6j(Exact(k),1,1,x,1,1,1)
        total = total + (-1)^(9+x) * qdim(Exact(k),x) * v^3
    end
    total
end
rhs_be = q6j(Exact(k),1,1,1,1,1,1)^2
residual_be = lhs_be - rhs_be
@assert isempty(residual_be) # hide
isempty(residual_be)
```

This verifies one admissible coefficient identity, not the full theorem for arbitrary labels. For general labels, derive the summation bounds, parity, and level truncation from the three left-hand symbols and check the right-hand symbols' admissibility. The test suite contains broader label samples.

## Prove a generic identity in x

The same diagonal orthogonality relation can be compared before choosing a level:

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

`QRecoupling.generic_sixj` is a lower-level alternative that takes **doubled integer labels**. Prefer `x_form(q6j(Symbolic(), ...))` when working with physical spins consistently.

## When a zero test is a certificate

`ExactXSum` keeps terms in radical classes. At a particular level, different generic classes can represent dependent roots. `iszero`/`==` use structural and norm tests first, then exact algebraic numbers in the selected real embedding when those tests are inconclusive. A `true` result is an exact zero decision, independent of numerical tolerance or rational rescaling. The fallback can be more expensive for high-degree values.

An empty residual is the cheapest certificate. A nonempty residual can also vanish after specialization; use `iszero` to resolve it. Generic polynomial certificates instead prove an identity before specialization. `iszero_at` performs exact confirmation of unresolved modular candidates itself. The `:cancels` entries of `level_spectrum` remain screening results and can be confirmed with `iszero_at` or `Exact(k)`.
