# Accuracy and performance

## Choose the computation before the precision

For one coefficient, use a scalar symbol call. For many coefficients at one target, use a batch. For a whole fusion transformation, use `fmatrix`: its recurrence shares work across entries. Request exact or generic symbolic algebra when that algebra is needed; it is not a prerequisite for numerical evaluation.

```@example perf
using QRecoupling
work = EvaluationWorkspace()
v = q6j(10,10,10,10,10,10; k=50, workspace=work)
again = q6j(10,10,10,10,10,10; k=50, workspace=work)
@assert v == again # hide
labels = [(j,j,j,j,j,j) for j in 1:8]
q6j(labels; k=30, threads=1)
```

One workspace can be reused sequentially. Do not share it across concurrent tasks. For repeated `qcg` or `q3j` calls, it also retains the CG recurrence coefficients for the current `(j1,j2,m1+m2)` sector and reuses scratch arrays. Changing the target or sector refreshes the coefficients automatically. For an entire row or sector, prefer `qcg_row` or `qcg_matrix` so that work is shared across entries. A call without a workspace borrows column scratch kept per thread, so a warm 6j evaluated along its column need not allocate. Batches own worker scratch and reuse tables. On Julia before 1.12, batches run serially because BigFloat precision scopes are shared. Exact batches also remain serial. Avoid clearing caches while evaluations are running.

Warm calls may avoid particular buffer allocations; this is **not** a package-wide zero-allocation guarantee. Table construction, BigFloat fallbacks, exact arithmetic, returned arrays, and symbolic expansion all allocate. Report both cold and warmed timings when benchmarking, and include labels, target, precision, Julia version, and thread count.

## Sweeping levels or deformation parameters

Repeated calls at one target reuse tables. A sweep that changes `k` or `q` at every point must build new numerical data, even when it reuses a workspace. Classical spin sweeps can reuse the same growing factorial tables, so they are a different workload.

Large Float64 level tables now use reflection and addition formulas with a small set of wide-precision seeds, rather than a wide number for every entry. Positive-real-q recurrence tables use the same approach within their exponent range. Small tables and inputs outside that range retain wide-precision construction. This reduces setup cost without removing compensated passes or higher-precision fallbacks.

```@example perf
function scaled_sweep(n)
    work = EvaluationWorkspace()
    [i^1.5 * q6j(Level(10i), (i .* (3,3,3,2,2,2))...; workspace=work)
     for i in 1:n]
end
samples = scaled_sweep(20)
@assert all(isfinite, samples) # hide
```

For real q, replace `Level(10i)` with `At(exp(-pi/(10i)))`. Compute the samples before plotting so that evaluation time is separate from plotting and editor communication. A workspace reuses scratch; it does not make tables for different targets interchangeable. Cache limits bound retention within individual caches, not the total process memory.

On an Apple M2 Pro, this sweep through `n = 2500` takes about 2.5 s at a level and about 5 s at real q. A real-q value very close to a node of its recurrence column cannot be certified from the recurrence, so it still falls back to arbitrary precision, at up to about 500 ms per call at these spins. Sweeps over complex q still build a wide-precision table at every point.

## Precision and supplied inputs

Numerical kernels use error estimates, compensated arithmetic, and higher-precision recomputation where needed. Level evaluation can move from Float64 through multiword arithmetic to BigFloat. Classical exact evaluation uses integers rather than trying to recover exact values from floating-point answers.

```@example perf
high = setprecision(BigFloat, 256) do
    q = BigFloat("0.8")
    q6j(1,1,1,1,1,1; q=q, T=BigFloat)
end
@assert high isa BigFloat # hide
high
```

Construct q at the intended precision: converting an existing Float64 to BigFloat does not recover missing input digits. At generic q, `T` sets a floor on the working/output precision; a higher-precision q is not narrowed. Classically and at a level, `T` is the result type: `T = Float32` returns `Float32`. For level calls, use `Level(k; T=BigFloat)` inside the desired precision scope.

For phase and braid data, use `T=Complex{BigFloat}` inside a precision scope with `rmatrix`, `twist`, `tmatrix` and `bmatrix`. Both the phases and the basis transformations are computed at the requested precision. `Level(k; T=BigFloat)` also works for `rmatrix` and `twist`: the level phase determines the complex result type.

A floating approximation to a root of unity is evaluated as the number supplied; it is not snapped to a level. Use `Level(k)` or `Exact(k)` when the root is known exactly. Near a singularity even a one-ulp input change can matter. Adaptive agreement is not a universal interval certificate, and a small numerical result is not proof of an exact zero. Final conversion can still overflow or underflow the output type.

When generic-q evaluation needs an exact zero decision, it first evaluates the sum modulo a prime at the exact stored real or complex parameter, including any integer q-power weights. A nonzero residue proves the sum is nonzero and avoids growing rational arithmetic. A zero residue or an unusable modular image still requires exact confirmation. This speeds up exceptional zero checks; ordinary scalar calls may never need the filter.

## Current limits to account for

- **F-matrix precision:** Float64/ComplexF64 matrices use the shared column recurrence. BigFloat inputs or `T=BigFloat`/`T=Complex{BigFloat}` use scalar `fsymbol` evaluations to preserve the requested precision. This higher-precision path costs more per entry.
- **Exact values:** construction retains a compact rule. First display materializes its quantum-integer presentation, which can cost substantially more than construction even when the printed output is truncated. `==` and `iszero` can reject candidates using modular residues; proving equality or zero may require more work. Above h = 1100, display omits distribution-relation reduction while preserving the value. Arithmetic on unweighted level forms uses the native exact engine; weighted `qcg` and `q3j` forms do not yet support it. `float` rounds from the compensated pass when its bound leaves almost nothing beyond the final rounding (within 0.501 ulp), and otherwise from 128 bits, so it is more accurate than a `Level(k)` value and can cost more. Forms with powers of q always take the 128-bit route. Numerical conversion of sums from level arithmetic adapts precision to cancellation. See [Exact forms](tutorials/exact_forms.md) and [Identity checks](tutorials/identities.md).
- **Tiny values at real q:** away from |q| = 1, coefficients can shrink rapidly with the spins. Subnormal values can still be represented below `floatmin(Float64)`; sufficiently small nonzero values round to `0.0`. Use `T = BigFloat` when the floating-point exponent range is insufficient. This can also affect `tetrahedron` at large spins. A rounded numerical zero is not an exact-zero certificate; use the exact form's `iszero` test.
- **Zero queries:** `iszero_at` uses structural proofs and modular filters, then confirms unresolved candidates exactly. Exact confirmation can cost more at high degree and runs serially after threaded screening. Workers use independently writable Boolean storage. The `:cancels` entries of `level_spectrum` are screening results unless you pass `prove = true`.
- **Custom series at singular targets:** cancellation of poles between separate summands is not generally regularized. Keep denominator factorials away from roots when possible.
- **Caches:** modular sine tables distinguish numeric type and precision. Broader concurrent cache-clearing workflows are not supported; clear caches between evaluations.

These limits distinguish the numerical coefficients, exact algebraic representation, and proof facilities; none should be inferred solely from another.

## Inspecting memory use

`empty_caches!()` clears the numerical and legacy projection caches covered by that API; it is not a promise to release every x-polynomial cache or Julia allocation. Use it between benchmark phases, not inside a hot loop. Workspace capacity grows to fit requests and can be discarded when no longer needed.
