# Factorial rules and evaluation architecture

## A common mathematical input

Most recoupling symbols are a factorial prefactor times a finite alternating sum. `FactorialSum` represents that structure directly. Each triple `(a,b,c)` means `[a*z+b]!^c`, and prefactor pairs `n=>c` mean `[n]!^c` outside the sum.

```@example rules
using QRecoupling
s = FactorialSum(1:3; factors=[(2,0,1)])
@assert isapprox(qeval(s), 746) # hide
(qeval(s), qeval(s; q=0.8), qeval(s; k=20))
```

Repeated factors are combined at construction. Factorial arguments must remain nonnegative integers over the finite range. General integer slopes, negative exponents, alternating signs, and an optional square-root prefactor are supported.

```@example rules
r = FactorialSum(0:4;
    factors=[(2,0,1), (1,0,-2), (-1,4,1)],
    prefactor=[4=>-1], alternating=true)
symbolic = qeval(Symbolic(), r)
```

This is a finite factorial-product language. Arbitrary parameterized `(a;Q)_z` factors, infinite-series convergence, and arbitrary summand functions are not part of this rule type. See [Finite series](tutorials/finite_series.md) for hypergeometric examples and conventions.

## Different targets, different arithmetic

| Target | Main route for factorial sums |
|:--|:--|
| Classical numerical | Adjacent-term ratios and adaptive summation |
| Classical exact | Integer Horner summation and prime-exponent cancellation |
| Level numerical | q-integer tables, valuation handling, compensated precision tiers |
| Real/complex q | Scaled ratios, compensated summation, adaptive BigFloat recomputation |
| Exact level | Reduced polynomial arithmetic in x, with ψ radical classes |
| Generic symbolic | Deferred rule; optional rational-function expansion in x |

Eligible families use three-term recurrences; selected single symbols can use a guarded near-edge recurrence. Difficult or unsuitable cases fall back to summation. Bare products such as dimensions have specialized paths, so sharing an interface does not imply identical internal arithmetic.

Numerical evaluation does not first expand the symbolic rational function. This matters: polynomial expansion can be much more expensive than computing one value, and evaluating large expanded polynomials near x = 2 can itself be ill-conditioned.

## The role of DCR

A Deferred Cyclotomic Representation stores a prefactor, a first summand, and adjacent-term ratios as factored cyclotomic monomials. It remains useful for legacy projections and callback series with explicit q powers. It is no longer the main user-facing symbolic form, and this compatibility layer is planned for retirement.

```@example rules
d = symbolic.dcr
@assert d isa QRecoupling.DCR # hide
d
```

Construction of `.dcr` is lazy and cached. Printing a DCR shows its structure; it neither sums nor factors its expression. `phi_form(symbolic)` remains an explicit, potentially expensive cyclotomic factorization request. Prefer `x_form(symbolic)` for reciprocal formulas.

For callback series, return **cyclotomic monomials**, not numeric q-factorials:

```@example rules
d = qseries(1:5) do z
    QRecoupling.qfact_mono(z)
end
@assert isapprox(qeval(d), sum(factorial(z) for z in 1:5)) # hide
qeval(d)
```

The `*_mono` constructors are qualified compatibility helpers. For new factorial sums, prefer the compact `FactorialSum` or factorial-triple `qseries` interface. A callback series with an interior zero term cannot in general supply the next ratio and is rejected.

## Extending the symbol interface

`QRecoupling.symbol_rule(QRecoupling.SixJ(), labels...)` exposes a symbol's rule using physical spins. A new `QRecoupling.QSymbol` subtype supplies `QRecoupling.symbol_rule`, `QRecoupling.nlabels`, and `QRecoupling.level_admissible`; `QRecoupling.symbol_family` and `QRecoupling.symbol_of` connect optional recurrence and function dispatch. Implementing a new formula this way lets it reuse the existing rule evaluators.

`QRecoupling.build_dcr!` and `QRecoupling.CycloBuffer` remain advanced compatibility tools. Reusing the buffer reduces scratch allocation, but constructing and storing the returned DCR still allocates.
