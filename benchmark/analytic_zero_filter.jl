# From the repository root:
#   julia --project=. benchmark/analytic_zero_filter.jl
#   julia --project=. benchmark/analytic_zero_filter.jl --exact
#
# Default: warmed zero-check timings and allocations, plus ordinary scalar calls.
# --exact also measures the unfiltered rational fallback up to j=50 (several seconds).
# These internal zero-check costs are conditional: ordinary calls may never reach them.
# Edit SPINS, TARGETS and REPEATS below to change the experiment. No data files are overwritten.

using QRecoupling, Printf, Test
const QR = QRecoupling
const SPINS = (20,50,100)
const TARGETS = (0.8,0.7+0.4im)
const REPEATS = 100

function measure(f; repeats=REPEATS)
    for _ in 1:3; f(); end
    bytes = @allocated f()
    times = [(@elapsed begin
        for _ in 1:repeats; f(); end
    end)/repeats for _ in 1:7]
    sort!(times)
    return times[4], bytes # median of seven warmed samples
end

function run_benchmark()
    println("Julia ",VERSION,"; ",Sys.CPU_NAME,"; threads=",Threads.nthreads())
    println("Zero-check costs (not complete q6j calls):")
    @printf("%5s %-14s %12s %10s %14s\n","j","q","filtered μs","bytes","raw exact ms")
    for j in SPINS, q in TARGETS
        s = QR.sixj_sum(ntuple(_->2j,6)...)
        @test !QR._analytic_sum_iszero(s,q,0)
        seconds,bytes = measure(()->QR._analytic_sum_iszero(s,q,0))
        raw = if "--exact" in ARGS && j <= 50
            @test !QR._analytic_sum_iszero_exact(s,q,0)
            @sprintf("%.3f",1000*(@elapsed QR._analytic_sum_iszero_exact(s,q,0)))
        else
            "--"
        end
        @printf("%5d %-14s %12.3f %10d %14s\n",j,string(q),1e6*seconds,bytes,raw)
    end
    println("\nOrdinary warmed scalar calls (separate measurement):")
    @printf("%5s %-14s %12s %10s\n","j","q","q6j μs","bytes")
    for j in SPINS, q in TARGETS
        seconds,bytes = measure(()->q6j(ntuple(_->j,6)...;q))
        @printf("%5d %-14s %12.3f %10d\n",j,string(q),1e6*seconds,bytes)
    end
    # A weighted identity must still take an exact/proved-zero route, unlike nonzero candidates.
    identity = FactorialSum(0:2;factors=((1,0,-1),(-1,2,-1)),alternating=true)
    @test QR._analytic_sum_iszero(identity,0.7+0.4im,1)
    println("\nWeighted cancellation identity confirmed exactly.")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_benchmark()
