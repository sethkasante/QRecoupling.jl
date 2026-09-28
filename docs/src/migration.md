# Migrating to v0.4

## Numerical defaults now include the product helpers

Symbol functions, `qeval`, `qint`, `qfact`, and `qbinomial` evaluate classically when no q or
k is supplied. Code that used their return values as symbolic objects must request
`Symbolic()` explicitly.

```@example migration
using QRecoupling
q6j(1,1,1,1,1,1)                      # approximately 1/6
q6j(Exact(),1,1,1,1,1,1)          # exact classical coefficient
s = q6j(Symbolic(),1,1,1,1,1,1)       # rule-backed symbolic x-form
@assert qint(5) == 5 # hide
qint(Symbolic(),5)
```

| Earlier pattern | v0.4 replacement |
|:--|:--|
| `q6j(js...)` to construct a DCR | `q6j(Symbolic(), js...).dcr` |
| `qfact(n)` as a symbolic expression | `qfact(Symbolic(), n)` |
| `qfact(n)` inside a monomial callback | Prefer a factorial rule; otherwise `QRecoupling.qfact_mono(n)` |
| Exact level cyclotomic wrapper by default | `Exact(k)` returns `ExactX` |
| `eager=true` | Remove the keyword |

`qseries` remains a constructor: triples produce a `FactorialSum`, callbacks produce a DCR.
Callback functions must return cyclotomic monomials; numeric product defaults cannot be
substituted there. `SymbolicValue` is a rule-backed view, not a replacement monomial with
the same multiplication interface. Build product rules using factorial prefactors or
expand with `x_form` when algebra is needed.

## Exact classical evaluation

Prefer `q6j(Exact(), js...)` over `q6j(js...; exact=true)`, and
`qeval(Exact(), rule)` for a classical factorial sum. `Exact()` is an alias for
`Classical(exact=true)` and preserves the existing exact result types and fast evaluator.
It also works with the other recoupling symbols and product helpers.

`exact=true` and `Classical(exact=true)` remain supported without deprecation warnings
in v0.4. A future release may deprecate the keyword after a migration period; this
release does not remove it. `Exact()` takes no level or form option: use `Exact(k)` for
level arithmetic, or `Symbolic()` for a generic formula.

## Exact levels use x

```@example migration
v = q6j(Exact(10),1,1,1,1,1,1)
@assert v isa ExactX # hide
v
```

This also applies to `k=10, exact=true` for the recoupling and product functions. Read
coefficients with `v.x_value`, `xpolynomial`, and `radicand`, request nested square roots
with `radical(v)` (or `v.rad`), and
convert numerically with `Float64`, `BigFloat`, or `numeric_value`.
Printing shows the stored x-form; radicals and numerical approximations are explicit requests.
`radical` returns `RadExpr` or a reasoned `NoRadical` result.

`Exact(k; form=:canonical)` retains the deprecated cyclotomic representation for
compatibility. `form=:deferred`, `ExactLevelWorkspace`, `canonicalize_exact`, and dedicated
Biedenharn–Elliott verifier functions from development drafts are not current APIs. Write
identities from the public exact arithmetic as in the [identity tutorial](tutorials/identities.md).
Raw DCR/monomial projection remains a separate compatibility path; it is not automatically
converted into an x-form by `qeval(Symbolic(), object)`.

## Symbolic output is bounded and lazy

Use `x_form(value)` to expand a symbolic rule. The earlier `xvalue(value)` spelling is
deprecated; `value.x_value` on an exact level result instead exposes its stored polynomials.

Symbolic expressions display their finite factorial rule in x without carrying out the sum,
for both small and large labels. Radicals use ψ factors. `x_form(value)` explicitly expands and caches
the generic formula, returning an owned copy. `.dcr` constructs a compatibility DCR only
when accessed; printing that DCR does no φ-form arithmetic.

Dimensions and theta graphs use this interface too. R matrices and twists retain `QPhase`
objects with an explicit branch over x. A classical exact R value is an integer sign;
a symbolic or exact-level R value remains a phase.

## Targets and scratch storage

`Level(k)`, `Exact()`, `Exact(k)`, `At(q)`, and `Classical()` select evaluation. Conflicting targets
and keywords are rejected. `exact=true` requires a classical or level target; use
`Symbolic()` for generic formulas.

`EvaluationWorkspace()` can be reused sequentially by supported scalar calls. Batches
manage worker workspaces automatically. `eager=true` is deprecated and delegates to the
standard rule evaluator; it does not enable a different fast path.
