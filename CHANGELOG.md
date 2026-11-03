# Changelog

## v0.6.2

- Added TensorKit sectors. `SU2k(k)` is SU(2) at level `k` and `SU2q(q)` is SU(2) at real `q`, with F- and R-symbols from the package's kernels, through an extension that needs only TensorKitSectors and no Nemo. F-symbols are real, so planar tensors can stay real.

## v0.6.1

- Nemo 0.55, 0.56 and 0.57 are allowed, as well as 0.54.
- `radical`, `xform`, `xpolynomial` and `radicand` on the result of arithmetic, or on a value at q = 1, say why they cannot answer. They used to say that Nemo was missing, even with Nemo loaded.
- CI runs on Node 24 actions.

## v0.6.0

Version 0.5.2 was never registered, and its changes are listed here.

### Breaking

- `Exact(k)` and `Symbolic()` return a `QForm`, the value written in quantum integers, and need no backend. `Exact()` stays the exact value at q = 1. `SymbolicValue` is gone, and `ExactX` and `ExactXSum` are no longer exported.
- Nemo is optional, needed only for `radical`, `xform` and `prove_identity`.
- The old cyclotomic layer is gone: `CyclotomicMonomial`, `DCR`, the `project_*` functions, the callback form of `qseries`, `phi_form` and `splits_completely`. `qseries(factors, range)` is unchanged.
- `x_form` is now `xform`. The old name still works, with a deprecation warning.
- A sum of unlike roots at q = 1 is now an exact `ClassicalSum`, where it was a `BigFloat`.
- Removed `eager`, one-argument `xvalue` and the canonical exact paths. See [Migration](docs/src/migration.md).

### New

- `==` and `iszero` on forms are exact and expand nothing, at a level and at every q, and they cover forms with powers of q.
- Exact arithmetic at a level covers every symbol without powers of q, and needs no backend.
- Forms cover `qcg` and `q3j`, printed with their powers of q.
- `ClassicalSum` keeps sums at q = 1 exact: `1//2 + √15/6`, with exact `+`, `-`, `*`, `/`, `==` and `iszero`.
- With Nemo, `xform` on an `Exact(k)` form prints a polynomial in x = 2cos(π/(k+2)), such as `(2x² − 7)/3`.
- `fmatrix`, `bmatrix`, `qcg_matrix`, `qcg_row` and `fmatrix_labels` take `Level(k)`, `At(q)` and `Classical()`, and `smatrix` and `tmatrix` take `Level(k)`.
- Numerical calls infer their result type at a level and at q = 1.
- `numeric_value`, `float` and `big` accept classical exact values.
- `benchmark/minimal.jl` is a short timed and checked benchmark. The article describing the package is now the preferred citation.

### Faster

- A first call takes about 1 ms, where it compiled for seconds, at the cost of about 30 s more when the package is installed.
- Sweeps over levels and real q no longer use a `BigFloat` per table entry, and real q falls back to arbitrary precision far less often. `q6j(Level(10i), …)` over i ≤ 2500 went from 41 s to 2.5 s, and `q6j(At(exp(-π/10i)), …)` from 72 s to 5 s. Table entries are also more accurate at large levels.
- `Exact(k)` takes about 0.15 µs at any size, because a form writes what it prints when first shown. `float` of a form takes about 0.3 µs on small forms, where it took 18 µs, with the same correctly rounded value on every input tested.
- Different level forms compare in about 5 µs at spins 100, k = 401, where it took 270 ms.
- Exact zero decisions no longer use Nemo. Long sums at a level, sums with several segments and weighted sums are decided in big integers, up to about 100 times faster than through the exact field.

### Fixed

- `Exact(k)` printed the wrong sign for labels past the level: `qint(Exact(3), 7)` is −[2] and printed `[2]`. Symbols with admissible labels never reach this.
- Equal forms had different hashes, so `unique` and `Dict` kept duplicates, and hashing a `ClassicalResult` threw.
- `qdim(j; k)` above k/2 gives zero, where it returned [2j+1]. `T = Float32` is respected at real and complex q.
- `==` between a form and a float is exact, `Exact` rejects invalid levels when it is built, and `twist(Exact(), j)` is 1.
- `qint`, `qfact`, `qbinomial`, `qdim` and `theta_value` at q = ±i now always raise an error naming `Level(0)`, unless the product is identically 0 or 1. They used to answer for some labels only.
- Exact values at q = 1 with large labels threw `InexactError` on Windows, where GMP takes 32-bit words.
- Fixed loading of GMP on Julia 1.14 under Windows (thanks to Kristoffer Carlsson). CI also runs on Windows and macOS.

## v0.5.1

- Fixed large-spin level `qcg`/`q3j` accuracy and inconsistent phases at purely imaginary q.
- Improved small-spin evaluation, compensated coupling passes, recurrence fallbacks and exact zero checks.
- Reduced phase-table setup costs and limited table-cache retention during parameter sweeps.

## v0.5.0

- Added `qcg`, `qcg_matrix` and `qcg_row`, with recurrence-based evaluation of coupling families.
- **Breaking:** `q3j` now includes the quantum deformation weights; classical values are unchanged.
- **Breaking:** fewer names are exported; qualified access remains. See [Migration](docs/src/migration.md) for target support and convention changes.
- Improved numerical accuracy, zero confirmation and precision handling; added `level_spectrum(...; prove=true)`.

## v0.4.0

- **Breaking:** omitted q and k mean classical evaluation, including `qint`, `qfact` and `qbinomial`. `Exact(k)` returns `ExactX` in x = 2cos(π/(k+2)); generic formulas use `Symbolic()`. Conflicting targets and keywords are rejected. See the [migration guide](docs/src/migration.md).
- **Added:** `Exact()` for exact classical values; x-form arithmetic (`x_form`, `prove_identity`) and explicit `radical(v)`; direct factorial-rule kernels at q = 1, real and complex q, with adaptive precision; workspaces, batches, level grids and column recurrences; F/braid matrices and modular S/T data.
- **Fixed:** exact equality uses algebraic arithmetic rather than a tolerance; threaded `iszero_at` no longer loses flags; theta-graph normalisation (`theta_value(j,j,0) = qdim(j)`); `all_6j(canonical = true)` merges all 144 symmetries.
- **Deprecated:** `eager = true` and `Exact(k; form = :canonical)`.

## v0.3.4

Correctness release; upgrading from v0.3.3 or earlier is recommended.

- **Fixed:** the sign of exact and analytic 6j symbols (wrong in about half of cases), the phase of `q3j`, exact phases at q = e^{iπ/(k+2)}, several crashes, and silent rounding of invalid spins.
- Documenter and Test are no longer runtime dependencies; randomized cross-checks against independent references were added.
