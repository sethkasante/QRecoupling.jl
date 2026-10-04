# Changelog

## v0.5.0

- **Breaking:** `q3j` is now the quantum 3j symbol of U_q(sl₂), orthogonal at q ≠ 1; classical values are unchanged bit for bit. The previous factorial-substitution symbol is `QRecoupling.q3j_factorial`. Negative real q uses the branch of `complex(q)`, and `fmatrix(...; T)` treats `T` as a precision floor.
- **Breaking:** the export list shrinks from about 90 to 45 names; the others remain available as `QRecoupling.name` (see the migration page).
- **Added:** `qcg`, `qcg_matrix` and `qcg_row`: quantum Clebsch–Gordan coefficients (labels flat or as `(j, m)` pairs), with whole matrices and rows from three-term recurrences. `level_spectrum(...; prove = true)` confirms `:cancels` entries exactly.
- **Fixed:** zeros are proved before they are returned (exact confirmation of modular candidates at levels, classically and at generic q). Real-q tables are correctly rounded at large spin (worst errors ≤ 1e−14, from 2.3e−12), and values near |q| = 1 are about 5× faster.
- **Deprecations:** removal of `eager = true`, `Exact(k; form = :canonical)` and `xvalue` will be done in later versions. `Exact(k)` and `Symbolic()` for `q3j` and `qcg` are not supported yet.

## v0.4.0

- **Breaking:** omitted q and k mean classical evaluation, including `qint`, `qfact` and `qbinomial`. `Exact(k)` returns `ExactX` in x = 2cos(π/(k+2)); generic formulas use `Symbolic()`. Conflicting targets and keywords are rejected. See the [migration guide](docs/src/migration.md).
- **Added:** `Exact()` for exact classical values; x-form arithmetic (`x_form`, `prove_identity`) and explicit `radical(v)`; direct factorial-rule kernels at q = 1, real and complex q, with adaptive precision; workspaces, batches, level grids and column recurrences; F/braid matrices and modular S/T data.
- **Fixed:** exact equality uses algebraic arithmetic rather than a tolerance; threaded `iszero_at` no longer loses flags; theta-graph normalisation (`theta_value(j,j,0) = qdim(j)`); `all_6j(canonical = true)` merges all 144 symmetries.
- **Deprecated:** `eager = true` and `Exact(k; form = :canonical)`.

## v0.3.4

Correctness release; upgrading from v0.3.3 or earlier is recommended.

- **Fixed:** the sign of exact and analytic 6j symbols (wrong in about half of cases), the phase of `q3j`, exact phases at q = e^{iπ/(k+2)}, several crashes, and silent rounding of invalid spins.
- Documenter and Test are no longer runtime dependencies; randomized cross-checks against independent references were added.
