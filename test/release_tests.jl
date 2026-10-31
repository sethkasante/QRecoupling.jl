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

@testset "removed in v0.6" begin
    @test_throws MethodError Exact(10; form = :canonical)
    @test_throws MethodError q6j(1, 1, 1, 1, 1, 1; eager = true)
    @test_throws MethodError QR.xvalue(q6j(Symbolic(), 1, 1, 1, 1, 1, 1))
    @test !isdefined(QR, :CompositeExactResult) && !isdefined(QR, :project_exact)
    m = QR.qint_mono(3)
    @test_throws ArgumentError qeval(m; k = 10, exact = true)
    @test_throws ArgumentError qeval(QR.q6j_dcr(2, 2, 2, 2, 2, 2); k = 10, exact = true)
    # what stays: numerical and classical evaluation of the same objects
    @test qeval(m; k = 10) ≈ qint(3; k = 10)
    @test qeval(m; q = 1, exact = true) == 3
    @test string(Exact(10)) == "Exact(10)"
end

@testset "values in quantum integers" begin
    sshow(f) = string(f)
    @test sshow(qform(q6j, 1, 1, 1, 1, 1, 1)) == "([5] − 1)/[4]!"
    @test sshow(qform(q6j, 1, 1, 1, 1, 1, 1; k = 3)) == "−1/[2]²   (k = 3)"
    @test sshow(qform(q6j, 1, 1, 1, 2, 1, 1; k = 5)) == "[2]/[3]²   (k = 5)"
    @test sshow(qform(q6j, 1, 1, 1, 1, 1, 1; k = 5)) == "([2] − 1)/([2][3]²)   (k = 5)"
    @test sshow(qform(fsymbol, 1, 1, 1, 1, 1, 1; k = 3)) == "−1/[2]   (k = 3)"

    # the same form at every generic q, with q = 1 the classical value
    f = qform(q6j(Symbolic(), 2, 2, 2, 2, 2, 2))
    @test f(1) ≈ q6j(2, 2, 2, 2, 2, 2)
    for q in (0.7, 0.95, 1.3)
        @test f(q) ≈ q6j(At(q), 2, 2, 2, 2, 2, 2) rtol = 1e-12
    end

    # values at levels, composite ones included, with and without the relations
    for k in (5, 10, 12, 16, 30), l in ((1, 1, 1, 1, 1, 1), (2, 2, 1, 2, 2, 1), (2, 2, 2, 2, 2, 2), (3//2, 3//2, 1, 3//2, 3//2, 2))
        QR._admissible_at(q6j, k, l...) || continue
        @test float(qform(q6j, l...; k = k)) ≈ q6j(Level(k), l...) atol = 1e-13
        @test float(qform(q6j, l...; k = k, reduce = false)) ≈ q6j(Level(k), l...) atol = 1e-13
    end
    @test numeric_value(qform(q6j, 1, 1, 1, 1, 1, 1; k = 5); bits = 256) ≈ numeric_value(q6j(Exact(5), 1, 1, 1, 1, 1, 1); bits = 256)

    # the relations reveal zeros that no single term shows: reflection and a sporadic one
    @test iszero(qform(q6j, 5, 5, 5, 5, 5, 5; k = 20))
    @test iszero(qform(q6j, 1, 3//2, 3//2, 2, 3//2, 3//2; k = 10))
    @test iszero(qform(q6j, 1, 5//2, 5//2, 5//2, 3, 3; k = 10))
    @test !iszero(qform(q6j, 5, 5, 5, 5, 5, 5; k = 21))
    @test iszero(qform(q6j, 3, 3, 9, 3, 3, 9; k = 4))                       # inadmissible at the level
    @test string(qform(q6j, 3, 3, 9, 3, 3, 9; k = 4)) == "0   (k = 4)"

    # the distribution relation itself: [pb] ∏[jh/p] = p ∏[b + jh/p]
    qn(n, h) = sinpi(n / h) / sinpi(1 / h)
    for h in (12, 15, 30), p in (2, 3, 5)
        h % p == 0 || continue
        for b in 1:h÷2
            (any(j -> (b + j * h ÷ p) % h == 0, 0:p-1) || (p * b) % h == 0) && continue
            @test qn(p * b, h) * prod(qn(j * h ÷ p, h) for j in 1:p-1) ≈ p * prod(qn(b + j * h ÷ p, h) for j in 0:p-1)
        end
    end
    # each stored row is a true relation between positive numbers
    for h in (12, 30), (_, row, con) in QR._relations_at(h)
        @test prod(qn(n, h)^row[n] for n in 2:h÷2) ≈ Float64(con)
    end

    @test_throws ArgumentError qform(q6j, 1, 1, 1, 1, 1, 1; k = 3)(0.8)
    @test_throws DomainError qform(q6j, 1, 1, 1, 1, 1, 1)(-0.8)
    @test_throws ArgumentError numeric_value(qform(q6j, 1, 1, 1, 1, 1, 1))
end

@testset "exact level values without Nemo" begin
    QR = QRecoupling
    # the value in q-integers agrees with the exact field, at prime and composite h, with and without a root
    for (k, l) in ((0, ntuple(_ -> 0, 6)), (1, (1//2, 1//2, 0, 1//2, 1//2, 0)), (3, ntuple(_ -> 1, 6)), (5, (1, 1, 1, 2, 1, 1)),
                   (10, ntuple(_ -> 2, 6)), (11, (1//2, 1, 3//2, 2, 3//2, 1)), (16, (3, 3, 2, 3, 3, 2)),
                   (30, ntuple(_ -> 3, 6)), (41, ntuple(_ -> 10, 6)), (121, ntuple(_ -> 30, 6)))
        v = QR.level_value(QR._rule_for(q6j, l...), k)
        @test float(v) ≈ Float64(numeric_value(q6j(Exact(k), l...))) rtol = 1e-12
        @test length(v.coeffs) == QR.euler_phi(2(k + 2)) ÷ 2
    end
    # [1], …, [d] carry the arithmetic: at k = 10 the value (2√3 − 3)/3 is (2[3] − 5)/3
    @test QR.level_value(QR._rule_for(q6j, ntuple(_ -> 1, 6)...), 10).coeffs == [-5 // 3, 0 // 1, 2 // 3, 0 // 1]
    # zeros that only additive relations show, proved without the exact field
    for (k, l) in ((20, (5//2, 9//2, 5, 4, 5, 11//2)), (22, (2, 3, 3, 5//2, 7//2, 7//2)), (10, (1, 5//2, 5//2, 5//2, 3, 3)))
        @test iszero(qform(q6j, l...; k = k))
        @test iszero(QR.level_value(QR._rule_for(q6j, l...), k))
        @test iszero_at(k, l...)
    end
    @test !iszero(qform(q6j, 2, 3, 3, 5//2, 7//2, 7//2; k = 23))
    # the valuations behind the proof: the norm of [n] is the matching power of its prime
    for h in (4, 6, 12, 16, 30, 49, 60)
        F = QR.LevelField(h); as = [a for a in 1:h-1 if gcd(a, 2h) == 1]
        for n in 2:h÷2
            l, w = F.valuation[n]
            @test sum(log(abs(sinpi(a * n / h) / sinpi(a / h))) for a in as) ≈
                  (l == 0 ? 0.0 : F.d * w * log(l)) atol = 1e-9
        end
    end
end

@testset "display of exact level values" begin
    plain(v) = sprint(show, MIME("text/plain"), v)
    v = q6j(Exact(5), 1, 1, 1, 2, 1, 1)
    @test string(v) == "ExactX(k = 5, [2]/[3]²)"
    @test occursin("≈ 0.3568958", plain(v)) && occursin("= [2]/[3]²", plain(v))
    @test !occursin("x²", plain(v))
    # a radical when the level has one
    w = q6j(Exact(3), 1, 1, 1, 1, 1, 1)
    @test string(w) == "ExactX(k = 3, −1/[2]²)"
    @test occursin("(√5 − 3)/2", plain(w))
    # the form shown is the value held
    for (k, l) in ((5, (1, 1, 1, 1, 1, 1)), (10, (2, 2, 2, 2, 2, 2)), (12, (2, 2, 1, 2, 2, 1)), (30, (3, 3, 3, 3, 3, 3)))
        u = q6j(Exact(k), l...)
        @test float(qform(u.rule; k = k)) ≈ Float64(numeric_value(u))
    end
    # a radical beside the form when one exists and is short, and none when it does not exist
    @test occursin("(2√3 − 3)/3", plain(q6j(Exact(10), 1, 1, 1, 1, 1, 1)))
    @test !occursin("√", plain(q6j(Exact(5), 1, 1, 1, 1, 1, 1)))
    # the result of arithmetic has no rule, and still reads in q-integers
    a = q6j(Exact(5), 1, 1, 1, 1, 1, 1) + v
    @test a.rule === nothing && string(a) == "ExactX(k = 5, [3] − [2])"
    @test !occursin("x²", plain(v * v))
    @test (v * v).rule === nothing && v * v == q6j(Exact(5), 1, 1, 1, 2, 1, 1)^2
    # a long sum is cut, never replaced by a polynomial
    long = replace(string(q6j(Exact(40), ntuple(_ -> 10, 6)...)), "ExactX" => "")
    @test occursin("terms)", long) && !occursin("x", long)
    # factorials combine across runs of labels
    @test string(qform(q6j, 1, 1, 1, 1, 1, 1)) == "([5] − 1)/[4]!"
    @test !occursin(r"\\[15\\]!.*\\[15\\]!", string(qform(q6j, ntuple(_ -> 30, 6)...; k = 200)))
    @test string(q6j(Exact(10), 1, 5//2, 5//2, 5//2, 3, 3)) == "ExactX(k = 10, 0)"
    # equality does not look at the display
    @test v == q6j(Exact(5), 1, 1, 1, 2, 1, 1) && v.rule !== nothing
end

