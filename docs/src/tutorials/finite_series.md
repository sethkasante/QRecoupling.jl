# Finite factorial and hypergeometric series

The compact rule interface can evaluate finite sums of symmetric q-factorials directly.
It does not require constructing a DCR for numerical evaluation.

## A terminating hypergeometric example

At q = 1 the classical Vandermonde identity gives

```math
{}_2F_1(-n,-n;1;1)=\sum_{j=0}^n\binom{n}{j}^2=\binom{2n}{n}.
```

Its symmetric q-factorial deformation is easy to construct:

```@example finite_series
using QRecoupling

n = 5
s = qseries([(1, 0, -2), (-1, n, -2)], 0:n;
            prefactor = [n => 2])

# [n]!^2 * sum(1 / ([j]!^2 * [n-j]!^2), j=0:n)
qeval(s)                         # 252.0
qeval(Exact(), s)              # exact classical result
qeval(s; q=0.8)                   # positive real deformation
qeval(s; q=0.8 + 0.2im)          # complex deformation
qeval(s; k=20)                    # level evaluation
qeval(Symbolic(), s)             # bounded symbolic x-form

# Equivalent callback construction, useful as a cross-check:
d = qseries(j -> QRecoupling.qbinomial_mono(n, j)^2, 0:n)
@assert isapprox(qeval(s; q=0.8), qeval(d; q=0.8))
@assert isapprox(qeval(s), binomial(2n,n)) # hide
```

The q-deformed sum above is not generally `qbinomial(2n,n)`. The classical identity
does not transfer to q without the appropriate weights.

## The finite q-binomial theorem

Standard basic hypergeometric notation usually uses the base `Q`, whereas this package uses
the symmetric integer `[n] = (q^n-q^(-n))/(q-q^(-1))`. With `Q=q^2`,

```math
[n]!=q^{-n(n-1)/2}\frac{(Q;Q)_n}{(1-Q)^n},\qquad
{n\brack j}_Q=q^{j(n-j)}\frac{[n]!}{[j]![n-j]!}.
```

The finite q-binomial theorem at series argument 1 is therefore

```math
\prod_{r=0}^{n-1}(1+q^{2r})
=\sum_{j=0}^n q^{j(n-1)}\frac{[n]!}{[j]![n-j]!}.
```

The extra q-power can be represented by a cyclotomic monomial:

```@example finite_series
n = 5
qpower(m) = CyclotomicMonomial(1, m, Pair{Int,Int}[], 0)
b = qseries(j -> qpower(j*(n-1)) * QRecoupling.qbinomial_mono(n,j), 0:n)

q = 0.8 + 0.2im
lhs = qeval(b; q=q)
rhs = prod(1 + q^(2r) for r in 0:n-1)
@assert isapprox(lhs, rhs)
(lhs, rhs)
```

The callback constructs a DCR here because the compact factorial rule currently has
no field for a summation-dependent q-power. See the [DLMF definition of basic
hypergeometric series](https://dlmf.nist.gov/17.4) and the
[q-binomial theorem](https://dlmf.nist.gov/17.2#iii) for the standard conventions.

## Scope and possible extensions

These are finite sums. Arbitrary parameter factors `(a;Q)_j`, an arbitrary series
argument `z^j`, and infinite-series convergence are not currently part of `FactorialSum`.
Terminating q-Chu–Vandermonde and q-Racah examples would be natural next tutorials
once parameterized q-Pochhammer factors and explicit q-power weights have a supported interface.

At a root of unity, distinguish a well-defined finite expression from individual terms
with poles. Cancellation of poles between summands is not generally supported by the
numerical factorial-rule evaluator. For the examples above, choosing `k=20` and `n=5`
keeps denominator factorials away from zeros.
