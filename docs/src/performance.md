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

- **F-matrix precision:** Float64/ComplexF64 matrices use the shared column recurrence.
  BigFloat inputs or `T=BigFloat`/`T=Complex{BigFloat}` use scalar `fsymbol` evaluations to
  preserve the requested precision. This higher-precision path costs more per entry.
- **Exact sums:** `ExactXSum` equality uses deterministic algebraic arithmetic, including
  dependent radical classes. The algebraic fallback can be expensive at high degree; see
  [Identity checks](tutorials/identities.md). Numerical conversion of a sum is a separate
  operation and does not yet adapt precision to cancellation between its terms.
- **Zero screening:** `iszero_at` and cancellation entries of `level_spectrum` can rely on
  modular screening. Treat cancellation zeros as candidates for exact confirmation rather
  than as a general proof. Threaded batch queries pack their Boolean results after the
  workers finish, so individual flags do not share writable storage.
- **Custom series at singular targets:** cancellation of poles between separate summands is
  not generally regularized. Keep denominator factorials away from roots when possible.
- **Caches:** modular sine tables distinguish numeric type and precision. Broader concurrent
  cache-clearing workflows are not supported; clear caches between evaluations.

These limits distinguish the numerical coefficients, exact algebraic representation, and
proof facilities; none should be inferred solely from another. The release review tracks
open issues separately from the documented supported workflows.

## Inspecting memory use

`empty_caches!()` clears the numerical and legacy projection caches covered by that API;
it is not a promise to release every x-polynomial cache or Julia allocation. Use it between
benchmark phases, not inside a hot loop. Workspace capacity grows to fit requests and can
be discarded when no longer needed.
