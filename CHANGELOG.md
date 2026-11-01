# Changelog

## v0.6.0 (in development)

- **Breaking:** `Exact()` and `Exact(k)` return a `QForm`, the value written in quantum integers from its formula. `Exact()` needs no q; used as a number it means its value at q = 1, and `classical` gives that exactly. `Exact(k)` needs no backend, with exact `==`, `iszero` and arithmetic at the level. `qcg` and `q3j` now have both, with powers of q. See [Migration](docs/src/migration.md).
- `==` between `Exact()` values decides equality at every q, and a sum that vanishes at every q is shown as 0. Nothing is expanded: an integer polynomial cannot vanish at an integer above the sum of its coefficients, so one such point is read through residues at word primes.
- `==` at a level expands nothing and covers forms with powers of q. It compares squares at places modulo word primes, from the rule's factorials: different values take 5 µs at spins 100, k = 401, where they took 270 ms, and a proof of equality takes 6.5 ms. The number of places comes from a bound on the norm, which cut them 7-fold at that size and also applies to zero tests. Sign comes from the certified number.
- Nemo is optional, and needed only for nested radicals and expanded polynomials.
- Added `qform` and quantum-integer displays for exact values, including native level zero checks. Polynomial views remain explicit.
- `qform` now covers `qcg` and `q3j`: powers of q are printed in front and on each term, reduced by q^h = −1 at a level, where these values are complex.
- The relations at a level are built from their distinct generators and kept reduced, with checked arithmetic. Levels such as k = 2000 no longer stall, and tables that took 384 ms at k = 1000 take 9 ms. At very smooth levels long rows are left out, so a form can stay longer there, and values are never affected.
- Square classes use the generators directly, which keeps them complete at every level with small entries.
- Certified native bounds, bad-prime rejection, stable `qform` evaluation, and quadratic-character nonzero certificates. Faster modular table setup and bounded cache retention.
- **Breaking:** removed the deprecated canonical/DCR exact paths, `eager`, and one-argument `xvalue`. Polynomial result types now parameterize their payload. See [Migration](docs/src/migration.md).

## Unreleased

- Exact zero decisions in numerical evaluation and zero queries no longer use Nemo. Long sums at a level, sums with several segments and weighted sums are decided in big integers, up to about 100 times faster than through the exact field.
- `numeric_value`, `float` and `big` accept classical exact values.
- CI also runs on Windows and macOS.
- Added `benchmark/minimal.jl`, a short timed and checked benchmark. The article describing the package is now the preferred citation.

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
