# Minimal benchmark of QRecoupling.jl: a few representative calls, each timed and checked.
#
#   julia --project=. benchmark/minimal.jl
#
# It needs only the package. Every case compares its value with an exact result, a higher-precision
# evaluation or an identity, and the script fails if a comparison does. Times depend on the machine and are
# printed for comparison between versions; they are not tested.
using QRecoupling, Printf

"""
Median time of one call of `f` in nanoseconds. A sample is the mean over `evals` consecutive calls, because
a single call at small spins is shorter than the clock resolution of some machines.
"""
function time_call(f; evals = 100, samples = 41)
    f()                                            # compile and build tables outside the timing
    t = Vector{Float64}(undef, samples)
    for s in 1:samples
        t0 = time_ns()
        for _ in 1:evals
            Base.donotdelete(f())
        end
        t[s] = (time_ns() - t0) / evals
    end
    return sort!(t)[(samples+1)÷2]
end

relerr(v, r) = Float64(abs(big(v) - big(r)) / abs(big(r)))
show_time(ns) = ns < 1e3 ? @sprintf("%7.1f ns", ns) : ns < 1e6 ? @sprintf("%7.2f µs", ns / 1e3) : @sprintf("%7.2f ms", ns / 1e6)

const RESULTS = Tuple{String,Float64,Bool,String}[]
"Record one case: its label, the time of `f`, whether `ok` holds, and a note shown beside the result."
function case(label, f, ok, note; kw...)
    push!(RESULTS, (label, time_call(f; kw...), ok, note))
end

# Labels are read from a `Ref` inside the timed call so that the compiler cannot evaluate it in advance.
const ONES = Ref((1, 1, 1, 1, 1, 1))
const J50 = Ref(ntuple(_ -> 50, 6))
const J20 = Ref(ntuple(_ -> 20, 6))
const J100 = Ref(ntuple(_ -> 100, 6))
const M3 = Ref((2, 2, 0, 1, -1, 0))
const CG = Ref((1, 1, 1, -1, 1, 0))
const W3 = Ref((10, 10, 10, 1, -1, 0))
const ZERO = Ref((1, 5 // 2, 5 // 2, 5 // 2, 3, 3))
const TOL = 2.0^-40                               # the package's acceptance threshold for numerical values
const Q09 = At(0.9)                               # targets are built once, outside the timed calls
const Q08 = At(0.8)

# ---- classical ------------------------------------------------------------------------------------------
case("6j   classical, spins 1", () -> q6j(ONES[]...), q6j(ONES[]...) ≈ 1 / 6, "equals 1/6")
let e = relerr(q6j(J50[]...), BigFloat(q6j(Exact(), J50[]...)))
    case("6j   classical, spins 50", () -> q6j(J50[]...), e <= TOL, @sprintf("error %.1e against the exact value", e))
end
case("3j   classical", () -> q3j(M3[]...), q3j(M3[]...) ≈ -1 / sqrt(5) && q3j(Exact(), M3[]...)^2 == 1 // 5,
     "equals -1/√5")

# ---- unitary level --------------------------------------------------------------------------------------
case("6j   level 10, spins 1", () -> q6j(Level(10), ONES[]...),
     numeric_value(q6j(Exact(10), ONES[]...)) ≈ q6j(Level(10), ONES[]...), "agrees with the exact level value")
let e = relerr(q6j(Level(401), J100[]...), q6j(Level(401; T = BigFloat), J100[]...))
    case("6j   level 401, spins 100", () -> q6j(Level(401), J100[]...), e <= TOL,
         @sprintf("error %.1e against BigFloat (strong cancellation)", e))
end
let e = relerr(q3j(Level(50), W3[]...), q3j(Level(50; T = BigFloat), W3[]...))
    case("3j   level 50, spins 10", () -> q3j(Level(50), W3[]...), e <= TOL, @sprintf("error %.1e against BigFloat", e))
end
case("6j   level 10, exact zero", () -> q6j(Level(10), ZERO[]...),
     q6j(Level(10), ZERO[]...) == 0 && iszero_at(10, ZERO[]...) && iszero(q6j(Exact(10), ZERO[]...)),
     "numerical, query and exact outputs agree on zero")

# ---- real deformation -----------------------------------------------------------------------------------
let e = relerr(q6j(Q09, J20[]...), q6j(At(big(0.9)), J20[]...))
    case("6j   q = 0.9, spins 20", () -> q6j(Q09, J20[]...), e <= TOL, @sprintf("error %.1e against BigFloat", e))
end
case("CG   q = 0.8", () -> qcg(Q08, CG[]...),
     qcg(Q08, CG[]...) ≈ sqrt(qdim(Q08, 1)) * q3j(Q08, 1, 1, 1, 1, -1, 0), "matches the 3j normalization")

# ---- families and exact output --------------------------------------------------------------------------
let (F, rows, cols) = fmatrix(40, 35, 30, 45; k = 300)
    defect = maximum(abs, transpose(F) * F - one(F))
    case("F    matrix 61 × 61, level 300", () -> fmatrix(40, 35, 30, 45; k = 300), defect < 1e-12,
         @sprintf("orthogonality defect %.1e", defect); evals = 5)
end
case("6j   Exact(10), spins 1", () -> q6j(Exact(10), ONES[]...), !iszero(q6j(Exact(10), ONES[]...)),
     "algebraic value"; evals = 10)

# ---- report ---------------------------------------------------------------------------------------------
println("QRecoupling ", pkgversion(QRecoupling), ", Julia ", VERSION, ", ", Sys.cpu_info()[1].model)
println(rpad("case", 34), rpad("time per call", 16), "check")
for (label, ns, ok, note) in RESULTS
    println(rpad(label, 34), rpad(show_time(ns), 16), ok ? "ok    " : "FAILED", "  ", note)
end
failed = count(r -> !r[3], RESULTS)
failed == 0 || error("$failed of $(length(RESULTS)) checks failed")
println("All ", length(RESULTS), " checks passed.")
