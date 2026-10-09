# Changelog

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
