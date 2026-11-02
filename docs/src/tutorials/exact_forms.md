# Exact values

`Exact` always means one number, exactly. `Symbolic` means the formula, a function of q.

| Call | What it is | Number? | Arithmetic? |
|:--|:--|:--|:--|
| `q6j(Exact(), js...)` | the value at q = 1 | yes | yes |
| `q6j(Exact(k), js...)` | the value at level k | yes | yes |
| `q6j(Symbolic(), js...)` | the formula | through `f(q)` | no |

Zero and equal mean different things on the two sides. For `Exact(k)` they mean at this level, and for `Symbolic()` at every q. A 6j symbol can vanish at a level and not vanish as a formula: no admissible label with spins up to 3 is zero at every q, while 119 of them are zero at level 10.

## A number at q = 1

`Exact()` is a sign and a rational under a root, with its own arithmetic.

```@example exact_forms
using QRecoupling
r = q3j(Exact(), 1//2, 1//2, 0, 1//2, -1//2, 0)
@assert r^2 == 1//2 && q6j(Exact(), 1, 1, 1, 1, 1, 1) == 1//6 # hide
(q6j(Exact(), 1, 1, 1, 1, 1, 1), r, r^2)
```

## The formula at every q

`Symbolic()` writes the formula in quantum integers `[n]` and never expands it. It holds no number. Evaluate it with `c(q)` at any q that [`At`](@ref) takes. An exact 1 gives the exact value at q = 1, while `c(1.0)` gives a float.

```@example exact_forms
c = q6j(Symbolic(), 1, 1, 1, 1, 1, 1)
@assert c(1) == 1//6 # hide
@assert c(0.8) ≈ q6j(At(0.8), 1, 1, 1, 1, 1, 1) # hide
(c, c(1), c(1.0), c(0.8), c(0.8 + 0.2im))
```

A formula has no arithmetic and no `float` of its own, since these need a value of q. It knows its rule but not its labels, so it does not move to a level: use `q6j(Exact(k), js...)`, which applies the level's fusion rules.

`==` means equal at every q, for two formulas and for a formula and a number, which is a constant function. `[6]` and `[3]!` are both 6 at q = 1, and neither of them is 6.

```@example exact_forms
@assert qint(Symbolic(), 6) != qfact(Symbolic(), 3) && qint(Symbolic(), 6) != 6 && qint(Symbolic(), 6)(1) == 6 && qint(Symbolic(), 1) == 1 # hide
(qint(Symbolic(), 6), qfact(Symbolic(), 3), qint(Symbolic(), 6) == qfact(Symbolic(), 3), qint(Symbolic(), 6) == 6, qint(Symbolic(), 6)(1) == 6)
```

Nothing is expanded to decide this. Clearing denominators gives an integer polynomial in q, and a nonzero one cannot vanish at an integer larger than the sum of its coefficients. So the test reads that polynomial at one such integer, through its residues at word primes, and stops at the first nonzero residue. A vanishing sum is shown as 0. A value can vanish at q = 1 without vanishing elsewhere, and then its form is not 0 while `v == 0` is true.

## A level

`Exact(k)` is the same value at q = e^{iπ/(k+2)}. The labels fold to `n ≤ (k+2)/2`, the relations among the q-integers simplify the form up to k + 2 = 1100, and zeros are found and proved.

```@example exact_forms
v = q6j(Exact(10), 1, 1, 1, 1, 1, 1)
```

```@example exact_forms
@assert iszero(q6j(Exact(10), 1, 3//2, 3//2, 2, 3//2, 3//2)) # hide
(q6j(Exact(3), 1, 1, 1, 1, 1, 1), q6j(Exact(10), 1, 3//2, 3//2, 2, 3//2, 3//2))
```

Two equal values can print differently, so `==` decides equality and the text does not. Nothing is expanded to do it. `==` compares the squares of the two values at places, primes where q is a root of unity, reading each from the rule's factorials. One different residue settles it, in about 5 µs at spins 100, k = 401. Equal squares need more places, as many as the norm of the difference allows, about 6.5 ms at that size, and the certified number then gives the sign.

Arithmetic with `*`, `+`, `-`, `/` and `^` is exact at a level, needs no backend, and gives a value that prints in the basis `[1], …, [d]`.

```@example exact_forms
a = q6j(Exact(12), 1, 1, 1, 1, 1, 1)
@assert a * inv(a) == 1 # hide
(a^2, a * inv(a) == 1)
```

A form holds its rule, which is all that `==`, `iszero` and `float` read, so it costs about 0.15 µs at any size. What it prints is written the first time it is shown, about 490 µs at spins 100, k = 401.

Use `Level(k)` for numbers. It costs about 1 µs at that size, and is accurate to about 10⁻¹⁶ in the median. `float(v)` rounds from more digits than `Level(k)` keeps: it takes under a microsecond on a typical form, and about 260 µs at that size, where the sum cancels heavily. `numeric_value(v; bits)` gives more digits.

```@example exact_forms
(Float64(v), numeric_value(v; bits = 128))
```

## Powers of q

`qcg` and `q3j` carry powers of q, printed in front and on each term. At a level they reduce by q^h = −1, and the value is complex.

```@example exact_forms
qcg(Exact(9), 3//2, 1//2, 1, 0, 3//2, 1//2)
```

`==` and `iszero` decide these forms like any other. They have no arithmetic at a level.

## Radicals

With Nemo loaded, `radical` rewrites a level value in nested square roots.

```@example exact_forms
using Nemo  # optional backend for nested radicals
radical(v)
```

```@example exact_forms
no_surd = radical(q6j(Exact(5), 1, 1, 1, 1, 1, 1))
@assert no_surd isa NoRadical # hide
@assert no_surd.kind == :none # hide
no_surd
```

A radical exists only when the square of the value has degree a power of two over ℚ, which is rare beyond small levels. `radical` returns a `RadExpr` when it succeeds, otherwise a `NoRadical` with a reason:

| `kind` | Meaning |
|:--|:--|
| `:none` | No expression in real nested square roots exists |
| `:untried` | The degree budget prevented the search |
| `:long` | The expression exceeds an explicitly requested length budget |
| `:failed` | The descent could not finish its required sign checks |

Here “radical” means **nested square roots**, not arbitrary nth roots. The default `degree_limit` is 32; `maxlen=0` imposes no expression-length cap. Increasing the degree limit can be very expensive in both time and memory.

## Polynomials in x

With Nemo loaded, `QRecoupling.xpolynomial(v)` and `QRecoupling.radicand(v)` give a level value as `P(x)√R(x)` with `x = 2cos(π/(k+2))`, and `x_form` expands a generic expression in `x = q + q⁻¹`. Both carry out the sum, which grows quickly with the spins and is not needed to read a value or to decide one.

```@example exact_forms
s = q6j(Symbolic(), 1, 1, 1, 1, 1, 1)
g = x_form(s)       # (x² − 3)/(x⁴ − 3x² + 2)
```

`Symbolic()` writes the formula without carrying out the sum, and `x_form` expands it. Radicals are displayed as products of ψ factors.

Generic square roots specify an algebraic expression. To evaluate with the package's numerical branch convention, use `s(q)` or a level target on the original symbolic value.

## API reference

```@docs
QForm
x_form
XValue
QRecoupling.ExactX
QRecoupling.ExactXSum
radical
NoRadical
numeric_value
```
