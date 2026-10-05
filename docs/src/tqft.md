# Recoupling symbols and conventions

## Symbols and label order

| Function | Meaning |
|:--|:--|
| `q6j(a,b,e,c,d,f)` | `{a b e; c d f}` |
| `q3j(j1,j2,j3,m1,m2,m3)` | Quantum 3j symbol `(j1 j2 j3; m1 m2 m3)`; m labels sum to zero |
| `qcg(j1,m1,j2,m2,j,m)` | Quantum Clebsch–Gordan coefficient `⟨j1 m1; j2 m2∣j m⟩`; also `qcg((j1,m1),(j2,m2),(j,m))` |
| `fsymbol(a,b,e,c,d,f)` | Fusion-basis transformation coefficient |
| `gsymbol(a,b,e,c,d,f)` | Package's dimension-weighted tetrahedral coefficient |
| `qdim(j)` | Quantum dimension `[2j+1]` |
| `rmatrix(a,b,c)` | Braiding eigenvalue in fusion channel c |
| `QRecoupling.tetrahedron(...)` | Tetrahedron evaluation |
| `QRecoupling.theta_value(a,b,c)` | Theta-graph evaluation |

These functions use classical values by default. `Level(k)`, `Exact()`, `Exact(k)`, `At(q)`, and `Symbolic()` select the other supported regimes. The graph helpers are currently qualified names rather than exports.

```@example symbols
using QRecoupling
k = 5
(qdim(1; k=k), q3j(1,1,1,1,-1,0; k=k),
 fsymbol(1,1,1,1,1,1; k=k), gsymbol(1,1,1,1,1,1; k=k))
```

The F convention is

```math
[F^{abc}_d]_{ef}=(-1)^{a+b+c+d}
\sqrt{[2e+1][2f+1]}\begin{Bmatrix}a&b&e\\c&d&f\end{Bmatrix}_q.
```

Keep the returned channel labels when using `fmatrix`: they specify the matrix ordering. At a unitary level the matrix is real orthogonal. At generic complex q the corresponding identity uses **transpose**, not conjugate transpose. Do not assume that a complex recoupling matrix is unitary.

## Exact values and phases

```@example symbols
v = q6j(Exact(10), 1,1,1,1,1,1)
@assert v isa ExactX # hide
v
```

Exact level recoupling coefficients use `ExactX`, in the real algebraic variable `x = 2cos(π/(k+2))`, with square-root factors as needed. They do not default to the old cyclotomic-field wrapper. See [Exact forms](tutorials/exact_forms.md). `q3j` and `qcg` have no `Exact(k)` or `Symbolic()` form yet: their level values are complex, in ℚ(ζ) up to a square root. Use `Level(k)` for them.

Braiding and twists carry q phases. `rmatrix(Symbolic(), ...)` and exact level phase calls retain a `QPhase`; these are not ordinary real x-polynomials.

```@example symbols
phase = rmatrix(Symbolic(), 1//2, 1//2, 1)
phase
```

The display makes the algebraic extension explicit: `u² − xu + 1 = 0`, with `u = q`. It retains the branch needed to distinguish q from its inverse. At a numerical target:

```@example symbols
(rmatrix(1//2, 1//2, 1; k=5), twist(1//2; k=5))
```

Generic complex square roots follow the package's balanced branch convention. Replacing products of roots by a principal root of their product can change a sign. Use the supplied symbol evaluators to preserve that convention.

Negative real `q` uses the same branch as `complex(q)` throughout the recoupling symbols, `qcg`, and `fmatrix`; the result type is complex even when the value is real. At a point on the negative real axis, `arg(q)=π` and the square root of a negative balanced factor is the positive imaginary root. Passing a nonzero imaginary part evaluates that complex parameter without snapping it to the axis.

## Admissibility and model boundaries

For recoupling symbols, level-inadmissible labels return zero by convention. A custom factorial series is a different object: a denominator can vanish at that level and produce a pole. `QRecoupling.level_admissible(QRecoupling.SixJ(), k, labels...)` checks representation constraints without evaluating a coefficient.

A local 6j coefficient is an ingredient of a state sum, not an entire Turaev–Viro or Ponzano–Regge invariant. Edge weights, summation ranges, triangulation data, normalization, and any boundary observables must be supplied by the application.
