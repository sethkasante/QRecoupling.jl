# Changelog

## v0.4.0

### Breaking changes and migration

- Symbol functions and `qeval` evaluate classically when q and k are omitted. The product functions `qint`, `qfact`, and `qbinomial` now follow the same default. Use `Symbolic()` for generic expressions; callback DCR series require explicit monomial constructors.
- `Exact(k)` and rule-backed symbol/product calls with `k=k, exact=true` return the real x-basis representation `ExactX`, with `x=2cos(π/(k+2))`. Generic symbolic results use `SymbolicValue`; sums of exact level radical classes use `ExactXSum`.
- Conflicting targets and evaluation keywords are rejected. Generic numerical q with `exact=true` is rejected; request `Symbolic()` for a generic formula.
- Classical exact R values are integer phases; symbolic and exact-level braiding retain `QPhase`.
- Development-only deferred exact carriers/workspaces and dedicated identity-verifier APIs have been removed. Compose identities using the public exact arithmetic instead.

### Added and improved

- `Exact()` requests exact classical evaluation, as shorthand for `Classical(exact=true)`. It reuses existing evaluators and result types for symbols, products, and factorial sums. `exact=true` remains supported without deprecation warnings in v0.4.

- Generic x-form arithmetic through `x_form`, `XValue`, `generic_sixj`, and `prove_identity`. Exact level coefficient access and numerical conversion. Default display shows the stored x-form without radical extraction or approximation; both are explicit requests.
- `radical(v)` / `v.rad` returns nested square roots or a reasoned `NoRadical` result (`:none`, `:untried`, `:long`, or `:failed`). `v.x_value` exposes stored polynomial pairs. Numerical approximation in display is opt-in with `:approximate=>true`.
- Symbolic display retains finite factorial rules in x without summation; expanded formulas have a bounded display. ψ factors identify unresolved radicals. Explicit expansion is cached and returns an owned copy. Compatibility DCRs are constructed lazily and display structurally.
- Direct ψ multiplicity counts from factorial indices, shared denominator products in exact cleared sums, and a rational-denominator shortcut for reduced exact arithmetic.
- Direct exact classical evaluation from factorial rules using integer Horner summation and prime-exponent cancellation, with BigInt ratios when machine integers would overflow. Accumulators are presized and prime powers batched in machine words: `q6j(Exact(), …)` is 1.3–1.5× faster with 3–8× fewer allocations for spins 20–300, with bit-identical results.
- Direct real-q and complex-q factorial-rule evaluation, scaled products, compensated sums, adaptive precision, and balanced complex prefactors. Near-classical inputs retain their supplied q rather than being snapped to 1. `T` is a precision floor at generic q.
- Reusable `EvaluationWorkspace` storage, shared fixed-target tables, batches and level grids. Eligible families and F matrices use recurrences; guarded near-edge evaluation can accelerate selected single symbols, with fallback to summation.
- F/braid matrices, modular S/T data, twists, total quantum dimension, Gauss sums, monodromy, and numerical Verlinde coefficients. A common `QSymbol` interface connects rules and queries.
- Exact numerical conversion honors requested precision for single `ExactX` values and rejects negative real radicands rather than clamping them to zero.
- Factorial rules also construct the 3j/6j/F/G/tetrahedron compatibility DCRs. General-slope DCR ratio construction accumulates exponents without expanding every q-integer separately.

### Fixed issues

- Exact-sum equality resolves dependent radical classes with exact algebraic arithmetic, rather than a numerical tolerance that could turn small nonzero values into zeros.
- Threaded `iszero_at` batches use independently writable Boolean storage before packing the returned `BitVector`, preventing lost flags at worker boundaries.
- Modular sine-table cache keys include the numeric type as well as precision.
- BigFloat F matrices use scalar `fsymbol` evaluation at the requested precision; machine-precision matrices retain the shared recurrence.
- Corrected the theta-graph denominator so the unsigned normalization satisfies `theta_value(j,j,0) = qdim(j)` instead of returning 1 for every admissible triple.
- Restored expanded x-form rendering after symbolic-display cleanup.
- `all_6j(canonical = true)` keeps one label set per class of all 144 6j symmetries (tetrahedral and Regge), as documented; previously only the 24 tetrahedral relabellings were merged. Canonical sweeps are about 6× faster.

### Deprecated

- `eager=true` delegates to the standard factorial-rule evaluator and warns.
- `Exact(k; form=:canonical)` retains legacy cyclotomic output and warns; the default is `:x`. DCR and low-level projection APIs remain compatibility interfaces in v0.4.

### Documentation and validation

- Reorganized the documentation around current targets, factorial rules, x-form arithmetic, and explicit label/phase conventions. Added getting-started, applications, accuracy/performance, and tensor-network pages.
- Replaced stale exact-output transcripts and an unfinished identity tutorial with executable examples and assertions for numerical matrices, exact residuals, and generic identities.
- Rewrote the README and migration guide, including the changed product-function defaults.
- Added targeted release regression tests while keeping the default package suite small. The full independent-reference and identity suites remain local development checks.

See [the migration guide](docs/src/migration.md) for replacements and examples.

## v0.3.4

Correctness release; upgrading from v0.3.3 or earlier is recommended.

### Fixed
- **Sign of exact and analytic 6j symbols** (wrong in about half of cases, including at |q| = 1): the radical now uses the balanced branch √(q^P ΠΨ_d) = q^{P/2} √(ΠΨ_d), continuous from q = 1.
- **`q3j`**: the DCR path (including q = 1) now includes the phase (−1)^{j₁−j₂−m₃}, and |m| > j or mismatched parity gives zero.
- **Values at q = e^{iπ/(k+2)}** keep exact phases and signs (complex results where the value is not real); Φ_{mh} (m ≥ 2) is no longer treated as zero; terms vanishing at Φ_h no longer give `NaN`.
- Crashes in `empty_caches!` (now exported), `fuse_root`, `QPhase * CompositeExactResult` and `rmatrix_mono`; `qseries` no longer truncates silently at an interior zero term.
- Spins that are not multiples of 1/2 raise an `ArgumentError` instead of being rounded.

### Packaging and tests
- Documenter and Test are no longer runtime dependencies.
- Randomized cross-checks against independent BigFloat/BigInt Racah formulas, Nemo, and the Biedenharn–Elliott and orthogonality identities, with regression tests for every fix.
