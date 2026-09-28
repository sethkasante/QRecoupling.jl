# Changelog

## v0.4.0 (unreleased)

This section describes the current release-candidate implementation, not a registered release.

### Breaking changes and migration

- Symbol functions and `qeval` evaluate classically when q and k are omitted. The product
  functions `qint`, `qfact`, and `qbinomial` now follow the same default. Use `Symbolic()`
  for generic expressions; callback DCR series require explicit monomial constructors.
- `Exact(k)` and rule-backed symbol/product calls with `k=k, exact=true` return the real
  x-basis representation `ExactX`, with `x=2cos(π/(k+2))`. Generic symbolic results use
  `SymbolicValue`; sums of exact level radical classes use `ExactXSum`.
- Conflicting targets and evaluation keywords are rejected. Generic numerical q with
  `exact=true` is rejected; request `Symbolic()` for a generic formula.
- Classical exact R values are integer phases; symbolic and exact-level braiding retain `QPhase`.
- Development-only deferred exact carriers/workspaces and dedicated identity-verifier APIs
  have been removed. Compose identities using the public exact arithmetic instead.

### Added and improved

- `Exact()` requests exact classical evaluation, as shorthand for `Classical(exact=true)`.
  It reuses existing evaluators and result types for symbols, products, and factorial sums.
  `exact=true` remains supported without deprecation warnings in v0.4.

- Generic x-form arithmetic through `x_form`, `XValue`, `generic_sixj`, and `prove_identity`.
  Exact level coefficient access and numerical conversion. Default display shows the stored
  x-form without radical extraction or approximation; both are explicit requests.
- `radical(v)` / `v.rad` returns nested square roots or a reasoned `NoRadical` result
  (`:none`, `:untried`, `:long`, or `:failed`). `v.x_value` exposes stored polynomial pairs.
  Numerical approximation in display is opt-in with `:approximate=>true`.
- Symbolic display retains finite factorial rules in x without summation; expanded formulas
  have a bounded display. ψ factors identify unresolved radicals. Explicit expansion is cached and
  returns an owned copy. Compatibility DCRs are constructed lazily and display structurally.
- Direct ψ multiplicity counts from factorial indices, shared denominator products in exact
  cleared sums, and a rational-denominator shortcut for reduced exact arithmetic.
- Direct exact classical evaluation from factorial rules using integer Horner summation and
  prime-exponent cancellation, with BigInt ratios when machine integers would overflow.
- Direct real-q and complex-q factorial-rule evaluation, scaled products, compensated sums,
  adaptive precision, and balanced complex prefactors. Near-classical inputs retain their
  supplied q rather than being snapped to 1. `T` is a precision floor at generic q.
- Reusable `EvaluationWorkspace` storage, shared fixed-target tables, batches and level grids.
  Eligible families and F matrices use recurrences; guarded near-edge evaluation can accelerate
  selected single symbols, with fallback to summation.
- F/braid matrices, modular S/T data, twists, total quantum dimension, Gauss sums, monodromy,
  and numerical Verlinde coefficients. A common `QSymbol` interface connects rules and queries.
- Exact numerical conversion honors requested precision for single `ExactX` values and rejects
  negative real radicands rather than clamping them to zero.
- Factorial rules also construct the 3j/6j/F/G/tetrahedron compatibility DCRs. General-slope
  DCR ratio construction accumulates exponents without expanding every q-integer separately.

### Fixed

- Corrected the theta-graph denominator so the unsigned normalization satisfies
  `theta_value(j,j,0) = qdim(j)` instead of returning 1 for every admissible triple.
- Restored expanded x-form rendering after symbolic-display cleanup.

### Deprecated

- `eager=true` delegates to the standard factorial-rule evaluator and warns.
- `Exact(k; form=:canonical)` retains legacy cyclotomic output and warns; the default is `:x`.
  DCR and low-level projection APIs remain compatibility interfaces in v0.4.

### Documentation and validation

- Reorganized the documentation around current targets, factorial rules, x-form arithmetic,
  and explicit label/phase conventions. Added getting-started, applications, accuracy/performance,
  and tensor-network pages.
- Replaced stale exact-output transcripts and an unfinished identity tutorial with executable
  examples and assertions for numerical matrices, exact residuals, and generic identities.
- Rewrote the README and migration guide, including the changed product-function defaults.

### Known limitations before registration

- `ExactXSum` equality can use a numerical fallback for dependent specialized radical classes;
  arbitrary multi-class zero tests are not universal exact certificates.
- Structural cancellation queries still use modular screening; batch `iszero_at` has a packed-bit
  concurrency issue. Use serial queries and exact confirmation for proof-oriented work.
- F-matrix columns currently use machine-precision storage before conversion to `T`; requesting
  BigFloat matrix elements alone does not increase their precision.
- Cache-key/concurrency coverage and arbitrary exact-sum numerical cancellation remain open.
  See the accuracy guide and release review before claiming broader guarantees.

See [the migration guide](docs/src/migration.md) for replacements and examples.

## v0.3.4 (released)

Correctness release. Several results were wrong without any warning; upgrading is recommended.

### Fixed
- **Exact and analytic 6j symbols had the wrong sign in about half of cases.** `evaluate_exact` and
  analytic evaluation took the principal square root of a radical that carries a phase. Both now use
  the balanced branch √(q^P ΠΨ_d) = q^{P/2} √(ΠΨ_d) with Ψ_d(q) = q^{-φ(d)} Φ_d(q²), the branch that is
  continuous from q = 1.
- **`q3j` through the DCR path (including `q = 1`) missed the phase (−1)^{j₁−j₂−m₃}**; it now agrees
  with the eager path and with the classical Wigner 3j symbol.
- **`q3j` with |m| > j or mismatched parity** returned nonzero values; it now returns zero.
- **Discrete projection dropped phases and signs.** Values at q = e^{iπ/(k+2)} now keep exact integer
  phases and the sign of every cyclotomic factor, returning `Complex` results when the value is not real
  (for example Σ q^z), and `qdim` above the level has the correct sign.
- **Φ_{mh} (m ≥ 2) was treated as zero** in the discrete and exact tables. Its value is now
  Λ̃(m) Π_{e | mh, h ∤ e} (q^{2e} − 1)^{μ(mh/e)}.
- **Terms vanishing at Φ_h** are skipped by valuation instead of producing `NaN` or spurious poles.
- `empty_caches!()` threw `UndefVarError`; it now clears every cache and is exported.
- `fuse_root` threw a `MethodError`.
- `QPhase * CompositeExactResult` threw a `FieldError`; integer powers of q now multiply exactly and
  half-integer powers raise an `ArgumentError`.
- `rmatrix_mono` returned a `String` for half-integer powers; it now returns a `QPhase`.
- `qseries`/`build_series` silently truncated a series at an interior zero term; it now raises an
  `ArgumentError` (trailing zeros are still allowed).
- `build_dcr!` returned `ZERO_DCR` when handed a reused buffer whose sign was zero.
- Spins that are not multiples of 1/2 are rejected with an `ArgumentError` instead of being rounded.
  Spin arguments now accept any `Integer`, `Rational` or `AbstractFloat`.
- The root-of-unity table cache is read under its lock.

### Added
- `^` with a non-negative integer exponent for `CompositeExactResult`.

### Packaging
- Documenter and Test are no longer runtime dependencies.
- Removed the stray `src/Project.toml`; `docs/build/` and `test/Manifest.toml` are no longer tracked.

### Tests
- Randomized cross-backend tests against independent BigFloat/BigInt Racah formulas (6j and 3j,
  level k and q = 1), cyclotomic values checked against Nemo, the Biedenharn–Elliott and orthogonality
  identities, and regression tests for every fix above.
