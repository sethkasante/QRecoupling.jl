# Changelog

## v0.4.0 (unreleased)

### Breaking API changes
- Symbol functions and `qeval` evaluate classically when `q` and `k` are omitted.
  Use `Symbolic()` as the first argument for the previous symbolic output.
- Conflicting targets (`k` together with `q`, or a target object together with
  competing evaluation keywords) now raise an error. Generic analytic requests
  with `exact=true` are rejected rather than silently returning numeric results.
- `rmatrix(...; exact=true)` without a level returns the classical integer phase;
  explicit symbolic requests and exact level requests retain `QPhase`.

### Changed
- `eager=true` warns and delegates to the standard factorial-rule evaluator,
  including its standard exact result type instead of the old eager result wrapper.
- 3j, 6j, F, G, and tetrahedron DCRs are constructed from their factorial rules.
- General DCR ratios with non-unit slopes accumulate cyclotomic exponents in one
  pass over indices instead of expanding each q-integer.
- Symbolic batch construction and explicit lowering with `qeval(Symbolic(), rule)`.

See `docs/src/migration.md` for examples. This is a development version, not a
registered release. Zero-certification and concurrency follow-ups remain open.


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
