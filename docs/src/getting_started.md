# Getting started

## Labels and conventions

Public symbol functions take physical spins, not doubled integers. Write `1//2` for a half-integer; integer and exactly half-integral floating-point inputs are also accepted. Labels that violate the symbol's admissibility conditions give zero. Invalid fractional spins are rejected.

The symmetric q-integer is `[n] = (qⁿ − q⁻ⁿ)/(q − q⁻¹)` for positive n, with classical limit n. At level k, `q = exp(iπ/(k+2))`. Allowed representation labels are `0, 1//2, …, k//2`; admissible triples also satisfy triangle, parity, and level-sum bounds.

The product helper `qint(0)` uses the package's convention **1**, not the usual algebraic q-integer zero. `qfact(0)` is also 1. Take this convention into account in custom formulas.

```@example start
using QRecoupling
j = 1//2
v = q3j(j, j, 0, j, -j, 0)
@assert isapprox(abs(v), 1/sqrt(2)) # hide
(v, qdim(j), qint(5), qfact(5), qbinomial(5, 2))
```

A Clebsch–Gordan coefficient ⟨j₁m₁; j₂m₂|j m⟩ takes its labels pair by pair, either flat or as `(j, m)` tuples; `m` defaults to `m₁ + m₂` in the flat form. `q3j` takes the rows of the 3j symbol, `(j₁ j₂ j₃; m₁ m₂ m₃)`.

```@example start
a = qcg(1//2, 1//2, 1//2, -1//2, 1, 0)          # ⟨½ ½; ½ −½|1 0⟩
b = qcg((1//2, 1//2), (1//2, -1//2), (1, 0))    # the same, as (j, m) pairs
@assert a === b && a ≈ 1/sqrt(2) # hide
(a, qcg(At(0.8), (1//2, 1//2), (1//2, -1//2), (1, 0)), q3j(1//2, 1//2, 1, 1//2, -1//2, 0))
```

## Numerical, exact, or symbolic?

```@example start
js = (1, 1, 1, 1, 1, 1)
a = q6j(js...)                          # classical, Float64
b = q6j(Exact(), js...)              # exact classical rational/radical
c = q6j(Level(10), js...)               # equivalent to k=10
d = q6j(Exact(10), js...)              # equivalent to k=10, exact=true
real_q = q6j(At(0.8), js...)
complex_q = q6j(At(0.8 + 0.2im), js...)
@assert isapprox(Float64(d), c) # hide
(a, real_q, complex_q)
```

Supply either `k` or `q`. Do not combine an explicit target with competing keywords. `Exact()` selects exact classical evaluation; `Exact(k)` selects exact level evaluation. The keyword `exact=true` and `Classical(exact=true)` remain supported for compatibility. For a generic formula, request `Symbolic()` rather than `exact=true` at a numerical q.

```@example start
s = q6j(Symbolic(), js...)
g = x_form(s)                           # explicit rational x-form with ψ radicals
@assert isapprox(qeval(s; q=0.8), real_q) # hide
g
```

`Symbolic()` also applies to `qint`, `qfact`, `qbinomial`, dimensions, and the other recoupling symbols. It retains and displays the rule without carrying out the sum. `value.dcr` constructs a compatibility DCR lazily; ordinary numerical use needs no DCR.

## Batches and level sweeps

```@example start
labels = [(1,1,1,1,1,1), (2,2,2,2,2,2)]
values = q6j(Level(10), labels)
grid = q6j(Level(6:10), labels)          # rows: labels; columns: levels
@assert size(grid) == (2,5) # hide
@assert values ≈ grid[:,end] # hide
(size(grid), values)
```

Use `all_6j(k=4, canonical=true)` to enumerate representatives under the supported 6j symmetries. Complete enumeration grows rapidly with k; start with small levels or restrict `jmax`.

For a complete change of fusion basis, use `fmatrix` instead of a loop over scalar F symbols. See [Tensor networks and modular data](tutorials/tensor_networks.md).
