@testset "correctness regressions" begin
    # A modular collision must not turn a nonzero exact sum into zero.
    s = FactorialSum(0:1; factors=((1,1,120),), alternating=true)
    tabs = (QR.ClassicalModTable(UInt64(5),2), QR.ClassicalModTable(UInt64(13),2))
    @test QR._sum_vanishes_mod(s,tabs) && !QR.is_classical_zero(s,tabs)
    @test qeval(Exact(), s) == 1-big(2)^120

    # Weights change zero identities; tiny nonzero results must survive cancellation.
    s = FactorialSum(0:2; factors=((1,0,-1),(-1,2,-1)), alternating=true)
    @test !QR._analytic_sum_iszero(s,0.8,0)
    @test iszero(QR.analytic_value(s,0.8; weight=(1,0)))
    q = 1+big(1)//big(2)^500
    @test qeval(s; q) ≈ BigFloat(2/(q+inv(q))-1) rtol=big"1e-60"

    # An actual level cancellation, including writes across a packed-bit boundary.
    cancelled = (1,3//2,3//2,2,3//2,3//2)
    @test iszero(q6j(Exact(10), cancelled...))
    labels = [isodd(i) ? cancelled : ONES for i in 1:65]
    @test iszero_at(10,labels; threads=4) == isodd.(1:65)
    @test level_spectrum(cancelled...; k=[3,10],prove=true) == [:inadmissible,:cancels]
    @test level_spectrum(ONES...; k=Int[],threads=4) == Symbol[]
    @test_throws DomainError iszero_at(-1,ONES...)

    # Dependent radical classes need exact equality, even after tiny rescaling.
    a = q6j(Exact(4), 0,0,0,1//2,1//2,1//2)
    b = q6j(Exact(4), 0,0,0,3//2,3//2,3//2)
    tiny = big(1)//big(10)^100
    @test a == b && iszero(tiny*(a-b)) && !iszero(tiny*(a+b))

    # Retain the sign on the negative-real branch and precision in matrix output.
    l = (1//2,1//2,0,1//2,1//2,1)
    @test fsymbol(l...; q=-0.8) == fsymbol(l...; q=complex(-0.8))
    @test real(fsymbol(l...; q=-0.8)) < 0
    smatrix(7) # populate a Float64 table first
    setprecision(BigFloat,128) do
        F, es, fs = fmatrix(1,1,1,1; q=big"0.8",T=Float64)
        @test eltype(F) === BigFloat
        @test F ≈ [fsymbol(1,1,e,1,1,f; q=big"0.8") for e in es, f in fs] rtol=big"1e-30"
        S, = smatrix(7; T=BigFloat)
        @test eltype(S) === BigFloat && precision(S[1,1]) == 128
        @test S[1,1] ≈ sqrt(big(2)/9)*sinpi(big(1)/9) rtol=eps(BigFloat)*4
    end

    # Scaled intermediates must handle both ends of the Float64 exponent range.
    for q in (nextfloat(0.0),floatmax(Float64))
        expected = Float64(-inv(BigFloat(q)+inv(BigFloat(q))))
        @test q6j(1//2,1//2,0,1//2,1//2,0; q) ≈ expected rtol=1e-13 atol=nextfloat(0.0)
    end

    # The imaginary-axis branch must agree across precision and scalar/family routes.
    q = 0.999im
    ref = setprecision(() -> q6j(ONES...; q=Complex{BigFloat}(q), T=BigFloat), BigFloat, 256)
    @test q6j(ONES...; q) ≈ ref rtol=1e-12
    F, es, fs = fmatrix(2, 3//2, 2, 3//2; q)
    @test F ≈ [fsymbol(2, 3//2, e, 2, 3//2, f; q) for e in es, f in fs] rtol=1e-12

    # F/G recurrence fallbacks must retain their dimension factors.
    l = (381,296,365,393,260,389) .// 2
    for f in (fsymbol, gsymbol)
        ref = setprecision(() -> f(l...; k=800, T=BigFloat), BigFloat, 256)
        @test f(l...; k=800) ≈ ref rtol=1e-13
    end
    @test q6j(Level(10), cancelled...) === 0.0
    @test q6j(2,2,2,3//2,3//2,3//2) === 0.0

    # Eviction must preserve results when a level table is rebuilt.
    cache = QR.QINT_F64_TABLES; limit = QR.LEVEL_CACHE_LIMIT[]
    try
        empty!(cache); QR.LEVEL_CACHE_LIMIT[] = 200_000
        values = [q6j(Level(k), ONES...) for k in (700,800,900)]
        @test cache.bytes <= 200_000
        @test values == [q6j(Level(k), ONES...) for k in (700,800,900)]
    finally
        QR.LEVEL_CACHE_LIMIT[] = limit; empty_caches!()
    end
end

@testset "zero decisions in integers" begin
    # cyclotomic polynomials, including prime powers and several primes
    @test QR._cyclotomic_coeffs(1) == [-1, 1]
    @test QR._cyclotomic_coeffs(12) == [1, 0, -1, 0, 1]
    @test QR._cyclotomic_coeffs(30) == [1, 1, 0, -1, -1, -1, 0, 1, 1]
    for n in (7, 16, 36, 105, 210)
        p = QR.phi_q(n)
        @test QR._cyclotomic_coeffs(n) == [Int(QR.coeff(p, i)) for i in 0:QR.degree(p)]
    end

    # the big-integer test against the exact field, on one segment and across gaps
    rule(l...) = QR._rule_for(q6j, l...)
    for (l, k) in (((1, 5//2, 5//2, 5//2, 3, 3), 10), ((1, 1, 1, 1, 1, 1), 10), ((5, 5, 5, 5, 5, 5), 20),
                   ((5, 5, 5, 5, 5, 5), 21), ((4, 4, 4, 4, 4, 4), 7), ((6, 6, 6, 6, 6, 6), 5), ((3, 3, 3, 4, 4, 4), 4))
        s = rule(l...)
        st, segs = QR.classify_at_level(s, k)
        st === :finite || continue
        @test QR._level_ring_zero(s, segs, k) == iszero(QR.exact_x(s, k))
    end

    # a rule with slopes and exponents the machine-word test does not take
    s = FactorialSum(0:3; factors = [(2, 1, 1), (-1, 3, -2), (1, 0, -1), (1, 2, -1)], alternating = true)
    for k in 1:12
        st, segs = QR.classify_at_level(s, k)
        st === :finite || continue
        @test QR._level_ring_zero(s, segs, k) == iszero(QR.exact_x(s, k))
    end

    # long sums decided by the theorem: six equal odd spins cancel at k = 4j, not at k = 4j + 1
    for j in (17, 41)
        s = rule(ntuple(_ -> j, 6)...)
        @test QR._level_ring_zero(s, last(QR.classify_at_level(s, 4j)), 4j)
        @test !QR._level_ring_zero(s, last(QR.classify_at_level(s, 4j + 1)), 4j + 1)
    end

    # a phase can carry the zero: without alternation, w = h is the alternating sum
    s = rule(1, 5//2, 5//2, 5//2, 3, 3)
    t = FactorialSum(s.zlo:s.zhi; factors = [(f.a, f.b, f.c) for f in s.fac], alternating = false)
    @test QR._level_ring_zero(t, (s.zlo:s.zhi,), 10, 12)
    @test !QR._level_ring_zero(t, (s.zlo:s.zhi,), 10, 0)
    @test QR._weighted_level_zero(t, 12, 10)

    # queries take the same route
    @test iszero_at(10, 1, 5//2, 5//2, 5//2, 3, 3)
    @test iszero_at(10, [(1, 5//2, 5//2, 5//2, 3, 3), (1, 1, 1, 1, 1, 1)]) == [true, false]
    @test level_spectrum(1, 5//2, 5//2, 5//2, 3, 3; k = 8:12, prove = true)[3] === :cancels
end

@testset "conversions of classical exact values" begin
    v = q3j(Exact(), 2, 2, 0, 1, -1, 0)                  # −√(1/5)
    @test float(v) === Float64(v) === q3j(2, 2, 0, 1, -1, 0)
    @test big(v) == BigFloat(v)
    @test float(typeof(v)) === Float64
    @test numeric_value(v; bits = 256) ≈ -sqrt(big(1) / 5) atol = big(2.0)^-250
    @test precision(numeric_value(v; bits = 256)) == 256
    @test iszero(numeric_value(zero(v)))
end
