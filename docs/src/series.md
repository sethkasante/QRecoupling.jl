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

`qcg` and `q3j` reuse this factorial machinery with additional q-power weights stored separately. Those weights are not yet part of the public `FactorialSum` constructor or its symbolic representation.

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

## Symbolic expansion

`qeval(Symbolic(), rule)` retains the factorial rule without evaluating the sum. `x_form(value)` expands it in x; `phi_form(value)` requests cyclotomic factorization in q. Both expansions are explicit, and factorization can be substantially more expensive. See [Exact values in x](tutorials/exact_forms.md) for examples.

## Extending the symbol interface

`QRecoupling.symbol_rule(QRecoupling.SixJ(), labels...)` exposes a symbol's rule using physical spins. A new `QRecoupling.QSymbol` subtype supplies `QRecoupling.symbol_rule`, `QRecoupling.nlabels`, and `QRecoupling.level_admissible`; `QRecoupling.symbol_family` and `QRecoupling.symbol_of` connect optional recurrence and function dispatch. Implementing a new formula this way lets it reuse the existing rule evaluators.
