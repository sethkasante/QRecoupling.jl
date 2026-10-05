# Exact values in x

For exact classical coefficients, use `Exact()`:

```@example exact_forms
using QRecoupling
c = q6j(Exact(), 1, 1, 1, 1, 1, 1)
@assert c == 1//6 # hide
r = q3j(Exact(), 1//2, 1//2, 0, 1//2, -1//2, 0)
@assert r^2 == 1//2 # hide
(c, r)
```

This evaluates at q = 1 with the direct rational/radical evaluator. It does not retain a symbolic parameter or construct a level field.

Use the same reciprocal variable for level values and generic formulas:

```@example exact_forms
using QRecoupling

v = q6j(Exact(10), 1, 1, 1, 1, 1, 1)
# (2√3 − 3)/3 = (2x² − 7)/3, where x = 2cos(π/12)

s = q6j(Symbolic(), 1, 1, 1, 1, 1, 1)
g = x_form(s)       # (x² − 3)/(x⁴ − 3x² + 2), where x = q + q⁻¹
```

`Exact(k)` performs arithmetic at a specified level. `Symbolic()` retains only a factorial rule, without evaluating the sum during construction or display. It displays a finite sum using `F(n) = ∏[r=1:n] U_{r-1}(x/2)`, the symmetric q-factorial written entirely in x. This avoids an expensive polynomial expansion just because a value was printed.

```@example exact_forms
large = q6j(Symbolic(), 8, 8, 8, 8, 8, 8)
# A finite factorial sum in x; no large expansion on display.

expanded = x_form(large)   # explicitly pay for the full generic sum
again = x_form(large)      # reuse the cached calculation; return an owned copy
nothing # hide
```

Dimensions and theta graphs use the same symbolic x-form interface. Radicals are displayed as products of ψ factors, where ψ_e(x) is the minimal polynomial of 2cos(2π/e). `phi_form(value)` requests cyclotomic factorization explicitly. Braiding and twist phases preserve their branch through `u² − x·u + 1 = 0`, with `u = q`: x alone cannot distinguish q from its reciprocal.

Generic square roots specify an algebraic expression. To evaluate with the package's numerical branch convention, use `qeval(s; q=...)` or a level target on the original symbolic value. A polynomial identity does not license arbitrary changes of square-root branches at complex q.

## Stored x-form and explicit radicals

`Exact(k)` displays the stored `P(x)√R(x)` form. Printing does not run radical extraction or compute a numerical approximation. Large polynomials are summarized by size; their exact coefficients remain available.

```@example exact_forms
v.x_value               # named pair (P, R) of stored polynomials
QRecoupling.xpolynomial(v)          # P coefficients, in ascending powers of x
QRecoupling.radicand(v)             # R coefficients, in ascending powers of x
radical(v)              # explicitly request a nested-square-root expression
```

`v.rad` is shorthand for `radical(v)`. It is a computation, not a stored display field. Treat the polynomials exposed by `v.x_value` as read-only; copy them before mutation. For sums, `s.x_value` gives the corresponding collection of polynomial pairs.

```@example exact_forms
no_surd = radical(q6j(Exact(5),1,1,1,1,1,1))
@assert no_surd isa NoRadical # hide
@assert no_surd.kind == :none # hide
no_surd
```

`radical` returns a `RadExpr` when it succeeds, otherwise a `NoRadical` with a reason:

| `kind` | Meaning |
|:--|:--|
| `:none` | No expression in real nested square roots exists |
| `:untried` | The degree budget prevented the search |
| `:long` | The expression exceeds an explicitly requested length budget |
| `:failed` | The descent could not finish its required sign checks |

Here “radical” means **nested square roots**, not arbitrary nth roots. A `:none` result is not a claim that no representation using more general radicals exists. The default `degree_limit` is 32; `maxlen=0` imposes no expression-length cap. Increasing the degree limit can be very expensive in both time and memory.

To request an approximate value alongside the stored x-form:

```@example exact_forms
show(IOContext(stdout, :approximate=>true), MIME"text/plain"(), v)
```

Approximation is now opt-in at every degree. The old display controls `:radicals` and `:radical_degree_limit` no longer select radical extraction; call `radical` explicitly.

## Numerical conversion

```@example exact_forms
Float64(v)
setprecision(BigFloat, 512) do
    BigFloat(v)
end
numeric_value(v; bits=128)
```

Single `ExactX` values use adaptive guard precision for their selected real embedding and round to the requested output precision. The polynomial cancellation estimate is not an interval certificate. Negative radicands are rejected rather than silently changed to zero. General factorial rules can produce complex values even at a level; a real x-value cannot silently replace the required phase.

The numerical evaluation of a sum of several exact radical classes has a separate cancellation problem; the single-value precision guarantee above should not be assumed for arbitrary `ExactXSum` expressions. Equality and `iszero` of a sum are still decided exactly, never by a numerical tolerance; see [Checking and proving identities](identities.md).

## API reference

```@docs
SymbolicValue
x_form
XValue
ExactX
ExactXSum
radical
NoRadical
numeric_value
```
