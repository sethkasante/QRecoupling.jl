# Accuracy and performance

## Choose the computation before the precision

For one coefficient, use a scalar symbol call. For many coefficients at one target, use
a batch. For a whole fusion transformation, use `fmatrix`: its recurrence shares work
across entries. Request exact or generic symbolic algebra when that algebra is needed;
it is not a prerequisite for numerical evaluation.

```@example perf
using QRecoupling
work = EvaluationWorkspace()
v = q6j(10,10,10,10,10,10; k=50, workspace=work)
again = q6j(10,10,10,10,10,10; k=50, workspace=work)
@assert v == again # hide
labels = [(j,j,j,j,j,j) for j in 1:8]
q6j(labels; k=30, threads=1)
```

One workspace can be reused sequentially. Do not share it across concurrent tasks.
Batches own worker scratch and reuse tables. On Julia before 1.12, batches run serially
because BigFloat precision scopes are shared. Exact batches also remain serial.
Avoid clearing caches while evaluations are running.

Warm calls may avoid particular buffer allocations; this is **not** a package-wide
zero-allocation guarantee. Table construction, BigFloat fallbacks, exact arithmetic,
returned arrays, and symbolic expansion all allocate. Report both cold and warmed timings
when benchmarking, and include labels, target, precision, Julia version, and thread count.

## Precision and supplied inputs

Numerical kernels use error estimates, compensated arithmetic, and higher-precision
recomputation where needed. Level evaluation can move from Float64 through multiword
arithmetic to BigFloat. Classical exact evaluation uses integers rather than trying to
recover exact values from floating-point answers.

```@example perf
high = setprecision(BigFloat, 256) do
    q = BigFloat("0.8")
    q6j(1,1,1,1,1,1; q=q, T=BigFloat)
end
@assert high isa BigFloat # hide
high
```

Construct q at the intended precision: converting an existing Float64 to BigFloat does
not recover missing input digits. At generic q, `T` sets a floor on the working/output
precision; a higher-precision q is not narrowed. For level calls, use `Level(k; T=BigFloat)`
inside the desired precision scope.

A floating approximation to a root of unity is evaluated as the number supplied; it is
not snapped to a level. Use `Level(k)` or `Exact(k)` when the root is known exactly.
Near a singularity even a one-ulp input change can matter. Adaptive agreement is not a
universal interval certificate, and a small numerical result is not proof of an exact zero.
Final conversion can still overflow or underflow the output type.

## Current limits to account for

- **F-matrix precision:** the current recurrence forms Float64/ComplexF64 columns before
  storing the requested matrix element type. `T=BigFloat` does not provide arbitrary-precision
  F-matrix entries. For a high-precision matrix, evaluate scalar `fsymbol` entries explicitly.
- **Exact sums:** coefficient arithmetic is exact, but `ExactXSum` equality can use a numerical
  fallback when different radical classes specialize to dependent roots. An empty exact
  residual or a generic polynomial identity is a stronger certificate; see
  [Identity checks](tutorials/identities.md).
- **Zero screening:** `iszero_at` and cancellation entries of `level_spectrum` can rely on
  modular screening. Treat cancellation zeros as candidates for exact confirmation rather
  than as a general proof. For batch `iszero_at`, use `threads=1` while the packed-bit output
  concurrency issue remains open.
- **Custom series at singular targets:** cancellation of poles between separate summands is
  not generally regularized. Keep denominator factorials away from roots when possible.
- **Caches:** concurrency and cache-key coverage remain release-review items. In particular,
  do not rely on mixed Float64/64-bit-BigFloat modular-table requests sharing a cache safely.

These limits distinguish the numerical coefficients, exact algebraic representation, and
proof facilities; none should be inferred solely from another. The release review tracks
open issues separately from the documented supported workflows.

## Inspecting memory use

`empty_caches!()` clears the numerical and legacy projection caches covered by that API;
it is not a promise to release every x-polynomial cache or Julia allocation. Use it between
benchmark phases, not inside a hot loop. Workspace capacity grows to fit requests and can
be discarded when no longer needed.
