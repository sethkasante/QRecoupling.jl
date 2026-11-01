# Exact values

A recoupling value already has a short exact description, its formula. `Exact()` and `Exact(k)` write that formula in quantum integers `[n]`, and never expand it.

```@example exact_forms
using QRecoupling, Nemo
c = q6j(Exact(), 1, 1, 1, 1, 1, 1)
```

## Every q, and q = 1

`Exact()` needs no value of q, and holds at all of them. Evaluate it with `c(q)` at a real positive q. Used as a number, it means its value at q = 1, and [`classical`](@ref) gives that value exactly.

```@example exact_forms
@assert c == 1//6 # hide
@assert c(0.8) ≈ q6j(At(0.8), 1, 1, 1, 1, 1, 1) # hide
(c == 1//6, classical(c), Float64(c), c(0.8))
```

`classical` returns a sign and a rational under a root, with its own arithmetic.

```@example exact_forms
r = q3j(Exact(), 1//2, 1//2, 0, 1//2, -1//2, 0)
@assert classical(r)^2 == 1//2 # hide
(r, classical(r), classical(r)^2)
```

Between two `Exact()` values, `==` means equal at every q. The integer 6 is `[6]` and `[3]!` at q = 1, but those are different functions of q.

```@example exact_forms
@assert qint(Exact(), 6) != qfact(Exact(), 3) # hide
(qint(Exact(), 6), qfact(Exact(), 3), qint(Exact(), 6) == qfact(Exact(), 3))
```

Nothing is expanded to decide this. Clearing denominators gives an integer polynomial in q, and a nonzero one cannot vanish at an integer larger than the sum of its coefficients. So the test reads that polynomial at one such integer, through its residues at word primes, and stops at the first nonzero residue. A vanishing sum is shown as 0. A value can vanish at q = 1 without vanishing elsewhere, and then its form is not 0 while `v == 0` is true.

## A level

`Exact(k)` is the same value at q = e^{iπ/(k+2)}. The labels fold to `n ≤ (k+2)/2`, the relations among the q-integers simplify the form, and zeros are found and proved.

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

Use `Level(k)` for numbers. It costs about 1 µs where the form costs 450 µs at spins 100, k = 401, and is accurate to about 10⁻¹⁶ in the median. `numeric_value(v; bits)` gives more digits.

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

`Symbolic()` retains a factorial rule, without evaluating the sum during construction or display. Radicals are displayed as products of ψ factors, where ψ_e(x) is the minimal polynomial of 2cos(2π/e). `phi_form(value)` requests cyclotomic factorization explicitly. Braiding and twist phases preserve their branch through `u² − x·u + 1 = 0`, with `u = q`: x alone cannot distinguish q from its reciprocal.

Generic square roots specify an algebraic expression. To evaluate with the package's numerical branch convention, use `qeval(s; q=...)` or a level target on the original symbolic value.

## API reference

```@docs
QForm
classical
SymbolicValue
x_form
XValue
ExactX
ExactXSum
radical
NoRadical
numeric_value
```
