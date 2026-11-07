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
        p = Nemo.cyclotomic(n, Nemo.polynomial_ring(Nemo.ZZ, "q")[2])
        @test QR._cyclotomic_coeffs(n) == [Int(Nemo.coeff(p, i)) for i in 0:Nemo.degree(p)]
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
    f = q3j(Symbolic(), 2, 2, 0, 1, -1, 0)                  # −√(1/5)
    v = f(1)
    @test v == q3j(Exact(), 2, 2, 0, 1, -1, 0)
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
    @test string(Exact(10)) == "Exact(10)"
end

# `float(f)` and `f(q)` evaluate the rule behind a form, so they cannot see a mistake in the form. These two read the
# form itself: its fields, and the text a reader gets.
"The expanded value from the Nemo backend, independent of the forms and of the plain-Julia kernel."
oracle(k, l...) = QR.exact_x(QR._rule_for(q6j, l...), k)

"The value of a form's own fields, at its level or at a real q."
form_value(f; q = nothing) = setprecision(BigFloat, 1024) do
    v = q === nothing ? QR._qform_level(f) :
        QR._qform_value(f, isone(q) ? (n -> BigFloat(n)) : (n -> (big(q)^n - big(q)^(-n)) / (big(q) - 1 / big(q))),
                        a -> sqrt(big(q))^a)
    v isa Complex ? ComplexF64(v) : Float64(v)
end
function text_value(str::AbstractString; k = nothing, q = nothing)
    sup = Dict('⁰' => '0', '¹' => '1', '²' => '2', '³' => '3', '⁴' => '4', '⁵' => '5', '⁶' => '6', '⁷' => '7', '⁸' => '8', '⁹' => '9', '⁻' => '-')
    t = replace(first(split(str, "   (k")), "−" => "-", "√" => "*sqrt")
    t = replace(t, r"[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+" => m -> "^" * join(sup[c] for c in m))
    t = replace(t, r"\[(\d+)\]!" => s"F(\1)", r"\[(\d+)\]" => s"Q(\1)")
    t = replace(t, r"(q\^\([^)]*\)|q\^-?\d+|q)(?=[FQ(*\s]*[FQ(])" => s"\1*", r"(?<=[)\d])\s*(?=q)" => "*")
    t = replace(t, r"\)\s*(?=[FQ(\d])" => ")*", r"(?<=\d)(?=[FQ(])" => "*", r"\s\*sqrt" => "*sqrt", r"^\*sqrt" => "sqrt", r"^-\*sqrt" => "-sqrt")
    t = replace(t, r"\*+" => "*")
    setprecision(BigFloat, 1024) do
        Q(n) = k !== nothing ? sinpi(big(n) / (k + 2)) / sinpi(1 / big(k + 2)) :
               isone(q) ? BigFloat(n) : (big(q)^n - big(q)^(-n)) / (big(q) - 1 / big(q))
        F(n) = prod(Q(i) for i in 2:n; init = BigFloat(1))
        v = Base.invokelatest(eval(:((Q, F, q) -> $(Meta.parse(t)))), Q, F, k !== nothing ? cispi(1 / big(k + 2)) : big(q))
        v isa Complex ? ComplexF64(v) : Float64(v)
    end
end
@testset "forms as written" begin
    # the helpers have to notice a wrong form, or the tests below prove nothing
    f = QR.qform(q6j, 1, 1, 1, 2, 1, 1; k = 5)
    w = QR._written(f)
    wrong = QR.QForm(f.k, f.rule, f.weight, true, QR.QWritten(w.k, w.sign, 7 * w.coef, w.out, w.radcoef, w.rad, w.terms, w.qpow, w.powers))
    @test float(wrong) ≈ float(f)
    @test form_value(wrong) ≈ 7 * float(f)
    @test text_value(string(wrong); k = 5) ≈ 7 * float(f)
    # at levels, prime and composite, with and without the relations, including roots and factorials
    for k in (3, 5, 10, 12, 16, 21, 30), l in collect(all_6j(k = min(k, 12)))[1:97:end]
        QR._admissible_at(q6j, k, l...) || continue
        reference = q6j(Level(k), l...)
        for reduce in (true, false)
            g = QR.qform(q6j, l...; k = k, reduce = reduce)
            iszero(g) && continue
            @test form_value(g) ≈ reference atol = 1e-12
            length(QR._written(g).terms) <= QR.QFORM_SHOWN_TERMS && @test text_value(string(g); k = k) ≈ reference atol = 1e-12
        end
    end
    # at generic q, q = 1 included
    for l in collect(all_6j(k = 10))[1:211:end], q in (0.7, 1, 1.3)
        g = QR.qform(q6j, l...)
        iszero(g) && continue
        reference = isone(q) ? q6j(l...) : q6j(At(q), l...)
        @test form_value(g; q = q) ≈ reference rtol = 1e-10
        length(QR._written(g).terms) <= QR.QFORM_SHOWN_TERMS && @test text_value(string(g); q = q) ≈ reference rtol = 1e-10
    end
end

@testset "forms with powers of q" begin
    # a wrong power of q has to show, in the form and in its text
    f = QR.qform(qcg, 3//2, 1//2, 1, 0, 3//2, 1//2; k = 9)
    w = QR._written(f)
    wrong = QR.QForm(f.k, f.rule, f.weight, true, QR.QWritten(w.k, w.sign, w.coef, w.out, w.radcoef, w.rad, w.terms, w.qpow + 2, w.powers))
    @test float(wrong) ≈ float(f)
    @test !(form_value(wrong) ≈ form_value(f))
    @test !(text_value(string(wrong); k = 9) ≈ form_value(f))
    labels = [(a, m1, b, m2, j, m1 + m2) for a in (0:6) .// 2 for b in (0:6) .// 2 for j in abs(a - b):1:a+b
              for m1 in -a:1:a for m2 in -b:1:b if abs(m1 + m2) <= j]
    for l in labels[1:7:end]
        g = QR.qform(qcg, l...)
        @test !iszero(g)
        for q in (0.7, 1.3)
            @test form_value(g; q = q) ≈ qcg(l...; q = q) rtol = 1e-11
            @test g(q) ≈ qcg(l...; q = q) rtol = 1e-11
            length(QR._written(g).terms) <= QR.QFORM_SHOWN_TERMS && @test text_value(string(g); q = q) ≈ qcg(l...; q = q) rtol = 1e-11
        end
        @test form_value(g; q = 1) ≈ qcg(l...) rtol = 1e-12
        @test g(1) ≈ qcg(l...) rtol = 1e-12
        for k in (1, 2, 3, 4, 6, 10, 15)
            v = qcg(l...; k = k)
            for reduce in (true, false)
                h = QR.qform(qcg, l...; k = k, reduce = reduce)
                @test iszero(h) == iszero(v)
                iszero(h) && continue
                @test form_value(h) ≈ v atol = 1e-11
                @test ComplexF64(numeric_value(h)) ≈ v atol = 1e-11
                @test float(h) ≈ v atol = 1e-11
                length(QR._written(h).terms) <= QR.QFORM_SHOWN_TERMS && @test text_value(string(h); k = k) ≈ v atol = 1e-11
            end
        end
    end
    for (a, b, c, m1, m2) in ((1, 1, 1, 1, -1), (2, 2, 2, 0, 0), (3//2, 1, 1//2, 1//2, 0), (5//2, 2, 3//2, -1//2, 1))
        m3 = -m1 - m2
        @test form_value(QR.qform(q3j, a, b, c, m1, m2, m3); q = 0.8) ≈ q3j(a, b, c, m1, m2, m3; q = 0.8) rtol = 1e-11
        k = 2Int(a + b + c)
        @test form_value(QR.qform(q3j, a, b, c, m1, m2, m3; k = k)) ≈ q3j(a, b, c, m1, m2, m3; k = k) atol = 1e-11
    end
    # the printed powers: a form with one term, at generic q and reduced by q^h = −1 at a level
    @test string(QR.qform(qcg, 1//2, 1//2, 1//2, -1//2, 1, 0)) == "q^(−1/2) √(1/[2])"
    @test string(QR.qform(qcg, 1//2, 1//2, 1//2, -1//2, 0, 0)) == "q^(1/2) √(1/[2])"
    @test endswith(string(QR.qform(qcg, 1//2, 1//2, 1//2, -1//2, 1, 0; k = 3)), "(k = 3, q = e^{iπ/5})")
    # outside the fusion rule
    @test iszero(QR.qform(qcg, 1, 0, 1, 0, 2, 0; k = 3))
    @test iszero(numeric_value(QR.qform(qcg, 1, 0, 1, 0, 2, 0; k = 3)))
    # terms with one monomial and different powers do not merge at generic q, and do merge when q^h = −1 aligns them
    h = 9; s = FactorialSum(0:2; factors = [(0, 2, 1)], alternating = false)
    @test length(QR._written(QR.qform(s; weight = (2h ÷ 3, 0))).terms) == 3
    # 1 + q^(2h/3) + q^(4h/3) is a sum of cube roots of unity: no two terms merge, and it vanishes
    @test length(QR._written(QR.qform(s; k = h - 2, reduce = false, weight = (2h ÷ 3, 0))).terms) == 3
    @test iszero(QR.qform(s; k = h - 2, weight = (2h ÷ 3, 0)))
    @test !iszero(QR.qform(s; k = h - 2, weight = (2h ÷ 3 + 2, 0)))
    @test iszero(QR.qform(FactorialSum(0:1; factors = [(0, 2, 1)], alternating = true); k = h - 2, weight = (2h, 0)))
    @test length(QR._written(QR.qform(FactorialSum(0:1; factors = [(0, 2, 1)], alternating = true); k = h - 2, weight = (h, 0))).terms) == 1
end

@testset "values in quantum integers" begin
    sshow(f) = string(f)
    @test sshow(QR.qform(q6j, 1, 1, 1, 1, 1, 1)) == "([5] − 1)/[4]!"
    @test sshow(QR.qform(q6j, 1, 1, 1, 1, 1, 1; k = 3)) == "−1/[2]²   (k = 3)"
    @test sshow(QR.qform(q6j, 1, 1, 1, 2, 1, 1; k = 5)) == "[2]/[3]²   (k = 5)"
    @test sshow(QR.qform(q6j, 1, 1, 1, 1, 1, 1; k = 5)) == "([2] − 1)/([2][3]²)   (k = 5)"
    @test sshow(QR.qform(fsymbol, 1, 1, 1, 1, 1, 1; k = 3)) == "−1/[2]   (k = 3)"

    # the same form at every generic q, with q = 1 the classical value
    f = q6j(Symbolic(), 2, 2, 2, 2, 2, 2)
    @test f(1) ≈ q6j(2, 2, 2, 2, 2, 2)
    for q in (0.7, 0.95, 1.3)
        @test f(q) ≈ q6j(At(q), 2, 2, 2, 2, 2, 2) rtol = 1e-12
    end

    # values at levels, composite ones included, with and without the relations
    for k in (5, 10, 12, 16, 30), l in ((1, 1, 1, 1, 1, 1), (2, 2, 1, 2, 2, 1), (2, 2, 2, 2, 2, 2), (3//2, 3//2, 1, 3//2, 3//2, 2))
        QR._admissible_at(q6j, k, l...) || continue
        @test float(QR.qform(q6j, l...; k = k)) ≈ q6j(Level(k), l...) atol = 1e-13
        @test float(QR.qform(q6j, l...; k = k, reduce = false)) ≈ q6j(Level(k), l...) atol = 1e-13
    end
    @test numeric_value(QR.qform(q6j, 1, 1, 1, 1, 1, 1; k = 5); bits = 256) ≈ numeric_value(oracle(5, 1, 1, 1, 1, 1, 1); bits = 256)
    @test setprecision(() -> QR._qform_level(QR.qform(q6j, 1, 1, 1, 1, 1, 1; k = 5)), BigFloat, 288) ≈
          numeric_value(oracle(5, 1, 1, 1, 1, 1, 1); bits = 288) rtol = big(2.0)^-256

    # the relations reveal zeros that no single term shows: reflection and a sporadic one
    @test iszero(QR.qform(q6j, 5, 5, 5, 5, 5, 5; k = 20))
    @test iszero(QR.qform(q6j, 1, 3//2, 3//2, 2, 3//2, 3//2; k = 10))
    @test iszero(QR.qform(q6j, 1, 5//2, 5//2, 5//2, 3, 3; k = 10))
    @test !iszero(QR.qform(q6j, 5, 5, 5, 5, 5, 5; k = 21))
    @test iszero(QR.qform(q6j, 3, 3, 9, 3, 3, 9; k = 4))                       # inadmissible at the level
    @test string(QR.qform(q6j, 3, 3, 9, 3, 3, 9; k = 4)) == "0   (k = 4)"

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
    for h in (12, 30), (_, row, con, _) in QR._relations_at(h)
        @test prod(qn(n, h)^row[n] for n in 2:h÷2) ≈ Float64(con)
    end
    # generators are distinct, and they span the relations: the rank is what the units leave free
    for h in (12, 30, 60, 77, 121, 210, 420, 714, 1002)
        @test length(QR._level_relations(h)) == (h ÷ 2 - 1) - (count(n -> gcd(n, h) == 1, 1:h) ÷ 2 - 1)
    end
    @test isempty(QR._level_relations(101))
    # smooth levels: finite in time and memory, every row true
    logq(n, h) = log(sinpi(big(n) / h)) - log(sinpi(1 / big(h)))
    for h in (910, 1050)
        t = @elapsed rows = QR._level_relations(h)
        @test t < 60 && !isempty(rows)
        @test all(r -> maximum(abs, r[2]) <= QR.RELATION_EXPONENT_LIMIT, rows)
        for (piv, row, con, _) in rows[1:29:end]
            @test all(iszero, row[piv+1:end]) && row[piv] > 0
            @test abs(sum(row[n] * logq(n, h) for n in 2:h÷2 if !iszero(row[n])) - log(BigFloat(con))) < 1e-30
        end
    end
    # the classes do not read the table, so they stay complete above it, with small entries
    for h in (1050, 2002, 2310)
        cl = QR._class_rows(QR._level_field(h))
        @test count(r -> r.pivot > 0, cl) == length(unique(r.pivot for r in cl if r.pivot > 0)) > 0
        @test maximum(maximum(abs, r.e) for r in cl) < 10^4
    end
    v = QR.level_value(QR._rule_for(q6j, 3, 3, 3, 3, 3, 3), 2000)
    @test float(v) ≈ q6j(Level(2000), 3, 3, 3, 3, 3, 3) rtol = 1e-12
    @test float(v * v) ≈ float(v)^2 rtol = 1e-12
    @test float(QR.qform(q6j, 3, 3, 3, 3, 3, 3; k = 2000)) ≈ q6j(Level(2000), 3, 3, 3, 3, 3, 3) rtol = 1e-12
    # above the limit there is no table to wait for, and a form still reads and evaluates right
    @test isempty(QR._level_relations(QR.RELATION_TABLE_LIMIT + 2))
    string(q6j(Exact(5000), 5, 5, 5, 5, 5, 5))
    @test haskey(QR.LEVEL_RELATIONS, 5002) && isempty(QR.LEVEL_RELATIONS[5002])   # printing asked, and nothing was built
    @test string(q6j(Exact(5000), 5, 5, 5, 5, 5, 5)) == string(QR.qform(q6j, 5, 5, 5, 5, 5, 5; k = 5000, reduce = false))
    @test QR._qform_level(q6j(Exact(5000), 5, 5, 5, 5, 5, 5)) ≈ numeric_value(q6j(Exact(5000), 5, 5, 5, 5, 5, 5)) rtol = 1e-12
    # an exponent past the range of Int stops a reduction, and never wraps
    @test_throws OverflowError QR._reduce_exponents!([0, typemin(Int) + 10, 100], [(3, [0, 4, 1], 1 // big(1), [2, 3])])

    @test_throws ArgumentError QR.qform(q6j, 1, 1, 1, 1, 1, 1; k = 3)(0.8)
    @test QR.qform(q6j, 1, 1, 1, 1, 1, 1)(-0.8) ≈ q6j(At(-0.8), 1, 1, 1, 1, 1, 1)
    @test_throws ArgumentError numeric_value(QR.qform(q6j, 1, 1, 1, 1, 1, 1))
end

@testset "exact level values without Nemo" begin
    QR = QRecoupling
    # the value in q-integers agrees with the exact field, at prime and composite h, with and without a root
    for (k, l) in ((0, ntuple(_ -> 0, 6)), (1, (1//2, 1//2, 0, 1//2, 1//2, 0)), (3, ntuple(_ -> 1, 6)), (5, (1, 1, 1, 2, 1, 1)),
                   (10, ntuple(_ -> 2, 6)), (11, (1//2, 1, 3//2, 2, 3//2, 1)), (16, (3, 3, 2, 3, 3, 2)),
                   (30, ntuple(_ -> 3, 6)), (41, ntuple(_ -> 10, 6)), (121, ntuple(_ -> 30, 6)))
        v = QR.level_value(QR._rule_for(q6j, l...), k)
        @test float(v) ≈ Float64(numeric_value(oracle(k, l...))) rtol = 1e-12
        @test length(v.num) == QR._totient(2(k + 2)) ÷ 2
    end
    # [1], …, [d] carry the arithmetic: at k = 10 the value (2√3 − 3)/3 is (2[3] − 5)/3
    @test QR.coordinates(QR.level_value(QR._rule_for(q6j, ntuple(_ -> 1, 6)...), 10)) == [-5 // 3, 0 // 1, 2 // 3, 0 // 1]
    # zeros that only additive relations show, proved without the exact field
    for (k, l) in ((20, (5//2, 9//2, 5, 4, 5, 11//2)), (22, (2, 3, 3, 5//2, 7//2, 7//2)), (10, (1, 5//2, 5//2, 5//2, 3, 3)))
        @test iszero(QR.qform(q6j, l...; k = k))
        @test iszero(QR.level_value(QR._rule_for(q6j, l...), k))
        @test iszero_at(k, l...)
    end
    @test !iszero(QR.qform(q6j, 2, 3, 3, 5//2, 7//2, 7//2; k = 23))
    # arithmetic agrees with the exact field, with and without a shared root
    let k = 12
        a = oracle(k, 1, 1, 1, 1, 1, 1); b = oracle(k, 3//2, 1, 3//2, 1, 3//2, 1)
        u = QR.level_value(QR._rule_for(q6j, 1, 1, 1, 1, 1, 1), k)
        w = QR.level_value(QR._rule_for(q6j, 3//2, 1, 3//2, 1, 3//2, 1), k)
        @test float(u * w) ≈ Float64(numeric_value(a * b))
        @test float(inv(w)) ≈ Float64(numeric_value(inv(b)))
        @test float(u + w) ≈ Float64(numeric_value(a + b))
        @test w * inv(w) == one(w) && isone(w * inv(w))
        @test (u + w) - w == u && u^3 * inv(u)^3 == 1
    end
    # the quadratic subfields: √5 at k = 3 is a field element, whose key can be merged
    let F = QR._level_field(5), G = only(QR._quadratic_roots(F))
        @test G[1] == 5 && QR._label_mul(G[2], G[2], F) == [big(5), big(0)]
    end
    @test QR.LevelValue(3, [big(1), big(0)], big(1), big(1), Int[]) * 2 == 2 * one(QR.LevelValue, 3)
    # a zero between two keys that the relations do not merge is still proved, by the separation bound
    let v = QR.LevelValue(3, [big(1), big(0)], big(1), big(5), Int[]),
        g = QR.LevelValue(3, QR._quadratic_roots(QR._level_field(5))[1][2], big(1), big(1), Int[])
        s = QR.LevelValueSum(v) - QR.LevelValueSum(g)
        @test length(s.terms) == 2 && iszero(s)
        @test !iszero(QR.LevelValueSum(v) - QR.LevelValueSum(2 * g))
    end
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

@testset "display of the expanded level values of the Nemo backend" begin
    plain(v) = sprint(show, MIME("text/plain"), v)
    v = oracle(5, 1, 1, 1, 2, 1, 1)
    @test string(v) == "ExactX(k = 5, [2]/[3]²)"
    @test occursin("≈ 0.3568958", plain(v)) && occursin("= [2]/[3]²", plain(v))
    @test !occursin("x²", plain(v))
    # a radical when the level has one
    w = oracle(3, 1, 1, 1, 1, 1, 1)
    @test string(w) == "ExactX(k = 3, −1/[2]²)"
    @test occursin("(√5 − 3)/2", plain(w))
    # the form shown is the value held
    for (k, l) in ((5, (1, 1, 1, 1, 1, 1)), (10, (2, 2, 2, 2, 2, 2)), (12, (2, 2, 1, 2, 2, 1)), (30, (3, 3, 3, 3, 3, 3)))
        u = oracle(k, l...)
        @test form_value(QR.qform(u.rule; k = k)) ≈ Float64(numeric_value(u))
        shown = replace(string(u), "ExactX(k = $k, " => "")[1:end-1]
        occursin("terms)", shown) || @test text_value(shown; k = k) ≈ Float64(numeric_value(u))
    end
    # a radical beside the form when one exists and is short, and none when it does not exist
    @test occursin("(2√3 − 3)/3", plain(oracle(10, 1, 1, 1, 1, 1, 1)))
    @test !occursin("√", plain(oracle(5, 1, 1, 1, 1, 1, 1)))
    # the result of arithmetic has no rule, and still reads in q-integers
    a = oracle(5, 1, 1, 1, 1, 1, 1) + v
    @test a.rule === nothing && string(a) == "ExactX(k = 5, [3] − [2])"
    @test !occursin("x²", plain(v * v))
    @test (v * v).rule === nothing && v * v == oracle(5, 1, 1, 1, 2, 1, 1)^2
    # a long sum is cut, never replaced by a polynomial
    long = replace(string(oracle(40, ntuple(_ -> 10, 6)...)), "ExactX" => "")
    @test occursin("terms)", long) && !occursin("x", long)
    # factorials combine across runs of labels
    @test string(QR.qform(q6j, 1, 1, 1, 1, 1, 1)) == "([5] − 1)/[4]!"
    @test !occursin(r"\[(\d+)\]![⁰¹²³⁴⁵⁶⁷⁸⁹]*\[\1\]!", string(QR.qform(q6j, ntuple(_ -> 30, 6)...; k = 200)))
    @test occursin(r"\[(\d+)\]![⁰¹²³⁴⁵⁶⁷⁸⁹]*\[\1\]!", "[15]!²[15]!")
    @test string(oracle(10, 1, 5//2, 5//2, 5//2, 3, 3)) == "ExactX(k = 10, 0)"
    # equality does not look at the display
    @test v == oracle(5, 1, 1, 1, 2, 1, 1) && v.rule !== nothing
end


@testset "v0.6 exact certification regressions" begin
    QR = QRecoupling
    F = QR._level_field(5)
    p = big(QR._places(F, 1, 5).m.p)
    # The first prime divides the whole element: reject it before updating CRT.
    a = [p, p]
    w, D = QR._label_inverse(a, F)
    @test QR._label_mul(a, w, F) == [D, big(0)]
    # A label outside the coordinate basis still has a cached inverse.
    F12 = QR._level_field(12)
    w, D = QR._label_inverse(5, F12)
    @test QR._label_mul(QR._label_coordinates(5, F12), w, F12) == [D; zeros(BigInt, F12.d - 1)]
    f = QR.qform(q6j, 2, 2, 2, 2, 2, 2)
    q = 1 + 1e-8
    ref = setprecision(() -> f(BigFloat(q)), BigFloat, 256)
    @test abs((f(q) - ref) / ref) < 1e-12
    # The enclosure must include the value, not merely agree to a tolerance.
    v = QR.level_value(QR._rule_for(q6j, 2, 2, 2, 2, 2, 2), 12)
    mid, err = setprecision(() -> QR._value_ball(v, 14), BigFloat, 128)
    exactref = numeric_value(oracle(12, 2, 2, 2, 2, 2, 2); bits = 512)
    @test abs(mid - exactref) <= err
    # Very small nonzero differences must not be classified as zeros.
    v = QR.LevelValue(3, [big(1), big(0)], big(1), big(5), Int[])
    g = QR.LevelValue(3, QR._quadratic_roots(F)[1][2], big(1), big(1), Int[])
    @test !iszero((QR.LevelValueSum(v) - QR.LevelValueSum(g)) + QR.LevelValueSum(one(g) / (big(1) << 600)))
end

@testset "prime-count bound" begin
    # the integer bound dominates the true largest conjugate of the cleared sum, and stays close to it
    for (k, l) in ((1, (1//2, 1//2, 0, 1//2, 1//2, 0)), (3, (1, 1, 1, 1, 1, 1)), (10, (2, 2, 2, 2, 2, 2)), (12, (3//2, 1, 3//2, 1, 3//2, 1)), (21, (5, 5, 5, 5, 5, 5)),
                   (30, (3, 3, 3, 3, 3, 3)), (41, (10, 10, 10, 10, 10, 10)), (59, (12, 12, 12, 12, 12, 12)), (121, (30, 30, 30, 30, 30, 30)))
        s = QR._rule_for(q6j, l...); h = k + 2
        zs, pre, _, es, cs = QR._level_exponents(s, k, Rational{BigInt})
        front = [fld(a, 2) for a in pre]
        L = QR._level_rule(s, zs, front, h)
        bits, J, C = QR._level_sum_bound(L, front, es, cs, h)
        truth = setprecision(BigFloat, 256) do
            largest = BigFloat(0)
            for a in 1:h-1
                gcd(a, 2h) == 1 || continue
                total = BigFloat(0)
                for (e, c) in zip(es, cs)
                    total += abs(BigFloat(c)) * prod((abs(sinpi(big(a * n) / h) / sinpi(big(a) / h))^(e[n] + front[n]) for n in 2:length(front)); init = BigFloat(1))
                end
                largest = max(largest, total)
            end
            Float64(log2(largest)) + J + C
        end
        @test truth <= bits <= truth + 16
    end
    # the table encloses the sines, and an exact 1 stays near zero
    let (up, lo) = QR._log_sine_table(60)
        for r in 1:30
            @test lo[r] / 2.0^20 <= log2(sinpi(r / 60)) <= up[r] / 2.0^20
        end
        @test up[30] - lo[30] <= 2
    end
    # integers only after the table: enclosing every conjugate of every term cost gigabytes per call here
    let s = QR._rule_for(q6j, ntuple(_ -> 30, 6)...), h = 202
        zs, pre, _, es, cs = QR._level_exponents(s, 200, Rational{BigInt})
        front = [fld(a, 2) for a in pre]
        L = QR._level_rule(s, zs, front, h)
        QR._level_sum_bound(L, front, es, cs, h)
        @test (@allocated QR._level_sum_bound(L, front, es, cs, h)) < 10^6
    end
    # a nonzero sum and a vanishing one both still finish with a proof
    @test !iszero(QR.level_value(QR._rule_for(q6j, ntuple(_ -> 30, 6)...), 121))
    @test iszero(QR.level_value(QR._rule_for(q6j, 5//2, 9//2, 5, 4, 5, 11//2), 20))
end

@testset "a form is written when it is shown" begin
    f = q6j(Exact(12), 2, 2, 2, 2, 2, 2); g = q6j(Exact(12), 2, 2, 2, 2, 2, 1)
    u = q6j(Symbolic(), 2, 2, 2, 2, 2, 2); v = q6j(Symbolic(), 2, 2, 2, 2, 2, 1)
    # deciding and evaluating read the rule
    @test !iszero(f) && f != g && f == f && f != 1 && u != v && !iszero(u) && u(1) == q6j(Exact(), 2, 2, 2, 2, 2, 2)
    float(f); numeric_value(f; bits = 80); u(0.8); u(1); -f; QR._form_value(f)
    @test all(x -> x.written === nothing, (f, g, u, v))
    # showing writes it once, and the text is what a form written at once would print
    text = string(f)
    w = f.written
    @test w isa QR.QWritten && string(f) == text && f.written === w
    @test text == QR._qform_body(QR._write(f.rule, f.k, true, f.weight)) * "   (k = 12)"
    # what writing alone can settle is still zero, and a rule that is singular says so at once
    @test iszero(q6j(Exact(10), 1, 3//2, 3//2, 2, 3//2, 3//2)) && string(q6j(Exact(10), 1, 3//2, 3//2, 2, 3//2, 3//2)) == "0   (k = 10)"
    @test iszero(qcg(Exact(9), 3//2, 1//2, 1, 0, 9//2, 1//2)) && -q6j(Exact(4), 1, 1, 1, 5, 5, 5) == 0
    @test_throws ArgumentError QR.qform(QR._qfact_rule([5 => -1]); k = 3)
end

@testset "values in containers, and comparing with floats" begin
    # equal values hash alike, whichever rule and printed form they came from
    for mk in (() -> q6j(Exact(10), 1, 1, 1, 1, 1, 1), () -> q6j(Symbolic(), 1, 1, 1, 1, 1, 1),
               () -> qcg(Exact(9), 3//2, 1//2, 1, 0, 3//2, 1//2), () -> q6j(Exact(4), 1, 1, 1, 5, 5, 5),
               () -> QR.level_value(QR._rule_for(q6j, 1, 1, 1, 1, 1, 1), 10), () -> oracle(10, 1, 1, 1, 1, 1, 1),
               () -> q6j(1, 1, 1, 1, 1, 1; exact = true), () -> q6j(1, 1, 1, 1, 1, 2; exact = true),
               () -> q6j(Symbolic(), 1, 1, 1, 1, 1, 1))
        a, b = mk(), mk()
        @test a == b && isequal(a, b) && hash(a) == hash(b) && length(unique([a, b])) == 1
    end
    fs = [q6j(Exact(18), a, b, c, d, e, x) for a in 0:1//2:2 for b in 0:1//2:2 for c in abs(a - b):a+b
          for d in 0:1//2:2 for e in abs(d - c):d+c for x in abs(a - e):a+e if abs(d - b) <= x <= d + b]
    fs = filter(!iszero, fs)
    @test length(unique(fs)) == length(unique(round.(float.(fs); sigdigits = 11)))
    different_rules = 0
    for f in fs, g in fs
        if f !== g && abs(float(f) - float(g)) < 1e-9 && QR._rule_key(f) != QR._rule_key(g)
            different_rules += 1
            @test (f == g) && hash(f) == hash(g)
        end
        different_rules > 200 && break
    end
    @test different_rules > 0
    @test length(Dict(q6j(Exact(10), 1, 1, 1, 1, 1, 1) => 1, q6j(Exact(10), 1, 1, 1, 1, 1, 1) => 2)) == 1
    @test hash(q6j(1, 1, 1, 1, 1, 1; exact = true)) == hash(1//6)

    # a form compares with a float as exactly as with a rational
    @test qint(Exact(3), 1) == 1.0 && 1.0 == qint(Exact(3), 1) && qint(Exact(3), 1) != NaN && qint(Exact(3), 1) != Inf
    f = q6j(Exact(3), 1, 1, 1, 1, 1, 1)
    @test f != float(f) && q6j(Symbolic(), 1, 1, 1, 1, 1, 1) != 1/6 && qint(Symbolic(), 1) == 1.0
    # a form and a number have no hash in common
    @test !isequal(qint(Exact(3), 1), 1) && !isequal(1, qint(Exact(3), 1))
    @test !isequal(q6j(Exact(3), 1, 1, 1, 1, 1, 1), q6j(Exact(4), 1, 1, 1, 1, 1, 1))
    @test Float32(f) ≈ float(f) && ComplexF64(f) ≈ float(f) && Float64[f, f] == fill(float(f), 2)
    @test convert(ComplexF64, qcg(Exact(9), 3//2, 1//2, 1, 0, 3//2, 1//2)) ≈ float(qcg(Exact(9), 3//2, 1//2, 1, 0, 3//2, 1//2))

    # a target says what is wrong with its level at once
    @test_throws DomainError Exact(-1)
    @test_throws DomainError Exact([3, -1])
    @test_throws MethodError Exact(1.5)
    @test Exact(3).k == 3 && Exact() == Classical(exact = true) && Exact(2:4).k == 2:4
end

@testset "Exact is one number, Symbolic is the formula" begin
    js = (1, 1, 1, 1, 1, 1); z = (1, 3//2, 3//2, 2, 3//2, 3//2)
    # Exact() is the released meaning, by the classical route
    @test Exact() == Classical(exact = true)
    @test q6j(Exact(), js...) isa QR.ClassicalResult && q6j(Exact(), js...) == 1//6 == q6j(js...; exact = true)
    @test float(q6j(Exact(), js...)) ≈ 1/6 && q6j(Exact(), js...)^2 == 1//36 && 6 * q6j(Exact(), js...) == 1
    @test qint(Exact(), 6) == 6 == qfact(Exact(), 3) && qdim(Exact(), 1//2) == 2 && qcg(Exact(), 1//2, 1//2, 1//2, -1//2, 1)^2 == 1//2
    @test q6j(Exact(), [js, (2, 2, 2, 2, 2, 2)]) == [q6j(Exact(), js...), q6j(Exact(), 2, 2, 2, 2, 2, 2)]
    @test rmatrix(Exact(), 1, 1, 1) == -1 && QR.twist(Exact(), 1) == 1 && QR.twist(1; exact = true) == 1
    @test Exact(1).k == 1 && iszero(q6j(Exact(1), js...))              # Exact(1) is level 1
    # Symbolic() is a function of q, wherever At is defined
    for (f, args) in ((q6j, js), (q6j, (2, 2, 2, 2, 2, 1)), (fsymbol, js), (qdim, (1,)), (qcg, (3//2, 1//2, 1, 0, 3//2, 1//2)),
                      (q3j, (2, 2, 0, 1, -1, 0)), (qint, (5,)), (qfact, (4,)), (qbinomial, (6, 2)))
        g = f(Symbolic(), args...)
        @test g isa QForm && g.k === nothing
        for q in (0.8, 1.3, 0.8 + 0.2im, -0.8, 1 + 0im)
            @test g(q) ≈ f(At(q), args...) rtol = 1e-11
        end
        @test g(1) == f(Exact(), args...)
        @test g == g && hash(g) == hash(f(Symbolic(), args...))
    end
    @test q6j(Symbolic(), [js, (2, 2, 2, 2, 2, 2)]) isa Vector{QForm} && QR.qeval(Symbolic(), q6j(Symbolic(), js...)) == q6j(Symbolic(), js...)
    s = FactorialSum(0:4; prefactor = (4 => 2,), factors = ((1, 0, -2), (-1, 4, -2)))
    @test qeval(Symbolic(), s)(1) == 70 == qeval(Exact(), s) && qeval(Exact(12), qeval(Symbolic(), s)) == qeval(Exact(12), s)
    @test qeval(qeval(Symbolic(), s); q = 0.8) ≈ qeval(s; q = 0.8)
    # a phase is still a phase
    @test rmatrix(Symbolic(), 1, 1, 1) isa QPhase && QR.twist(Symbolic(), 1) isa QPhase && QR.twist(Exact(10), 1) isa QPhase
    # zeros live at a level, and the formula shows the equation behind them
    @test iszero(q6j(Exact(10), z...)) && string(q6j(Exact(10), z...)) == "0   (k = 10)" && iszero(q6j(Level(10), z...))
    @test !iszero(q6j(Symbolic(), z...)) && string(q6j(Symbolic(), z...)) == "([7] − [2]²)/([3][4][5])"
    # the expansions read a formula, and refuse what they cannot expand
    @test xform(q6j(Symbolic(), js...)) isa XValue
    @test xform(q6j(Exact(10), js...)) isa LevelXValue
    @test_throws ArgumentError xform(qcg(Symbolic(), 1//2, 1//2, 1//2, -1//2, 1))
    @test_throws ArgumentError q6j(Symbolic(), js...; k = 3)
    @test_throws ArgumentError qcg(Symbolic(), 1//2, 1//2, 1//2, -1//2, 1; k = 3)
    @test_throws ArgumentError q6j(Exact(10), js...)(0.8)
end

@testset "numerical calls say what they return" begin
    # a level and q = 1 infer one concrete type through a level known only at run time
    symbols = ((q6j, (1, 1, 1, 1, 1, 1), Float64), (fsymbol, (1, 1, 1, 1, 1, 1), Float64), (gsymbol, (1, 1, 1, 1, 1, 1), Float64),
               (QR.tetrahedron, (1, 1, 1, 1, 1, 1), Float64), (QR.q3j_factorial, (1, 1, 2, 0, 0, 0), Float64), (QR.theta_value, (1, 1, 1), Float64),
               (qcg, (1, 0, 1, 0, 2, 0), ComplexF64), (q3j, (1, 1, 2, 0, 0, 0), ComplexF64), (qdim, (1,), Float64), (qint, (5,), Float64),
               (qfact, (4,), Float64), (qbinomial, (6, 2), Float64), (rmatrix, (1, 1, 1), ComplexF64), (twist, (1,), ComplexF64))
    for (f, a, E) in symbols
        target(k) = f(Level(k), a...); keyword(k) = f(a...; k = k); classical() = f(a...)
        @test only(Base.return_types(target, (Int,))) === E
        @test only(Base.return_types(keyword, (Int,))) === E
        @test isconcretetype(only(Base.return_types(classical, ())))
        @test f(Level(10), a...) isa E && f(a...; k = 10) === f(Level(10), a...)
        # at real q only a positive value is real, so a union is the honest answer, and the exact types stay out of it
        generic(q) = f(At(q), a...)
        @test only(Base.return_types(generic, (Float64,))) <: Union{Float64,ComplexF64}
    end

    # the element type is respected wherever a number is computed
    for (f, a) in ((q6j, (1, 1, 1, 1, 1, 1)), (fsymbol, (1, 1, 1, 1, 1, 1)), (gsymbol, (1, 1, 1, 1, 1, 1)), (QR.tetrahedron, (1, 1, 1, 1, 1, 1)),
                   (QR.theta_value, (1, 1, 1)), (qdim, (1,)), (qint, (5,)), (qfact, (4,)), (qbinomial, (6, 2)))
        for (kw, R) in (((;), Float32), ((; k = 10), Float32), ((; q = 0.8), Float32), ((; q = 0.8f0), Float32), ((; q = 0.8 + 0.2im), ComplexF32), ((; q = -0.8), ComplexF32))
            v = f(a...; kw..., T = Float32)
            @test v isa R
            @test v ≈ f(a...; kw...) rtol = 1e-6
        end
        @test f(a...; q = 0.8, T = BigFloat) isa BigFloat
        @test f(At(big(0.8)), a...) isa BigFloat                      # a wider q keeps its precision
    end
    for (kw, R) in (((;), Float32), ((; k = 10), ComplexF32), ((; q = 0.8), Float32), ((; q = 0.8 + 0.2im), ComplexF32))
        @test qcg(1, 0, 1, 0, 2, 0; kw..., T = Float32) isa R && q3j(1, 1, 2, 0, 0, 0; kw..., T = Float32) isa R
    end
    @test rmatrix(1, 1, 1; q = 0.8, T = ComplexF32) isa ComplexF32 && twist(1; q = 0.8, T = ComplexF32) isa ComplexF32
    s = QR._rule_for(q6j, 1, 1, 1, 1, 1, 1)
    @test qeval(s; q = 0.8, T = Float32) isa Float32 && qeval(s; k = 10, T = Float32) isa Float32
    rule_level(k) = qeval(s; k = k); rule_target(k) = qeval(Level(k), s)
    @test only(Base.return_types(rule_level, (Int,))) === Float64 && only(Base.return_types(rule_target, (Int,))) === Float64

    # a label above k/2 is not in the theory: every engine and every symbol agree it gives zero
    @test qdim(Level(4), 3) == 0 && qdim(3; k = 4) == 0 && qdim(3; k = 4, T = Float32) isa Float32
    @test iszero(qdim(Exact(4), 3)) && iszero(qdim(3; k = 4, exact = true)) && q6j(Level(4), 3, 3, 0, 3, 3, 0) == 0
    @test qdim(Level(4), 2) ≈ 1 && float(qdim(Exact(4), 2)) ≈ 1 && qdim(Level(4), 1) ≈ float(qdim(Exact(4), 1))
    @test qint(Level(4), 7) ≈ -1                                     # a q-integer has no label to reject
    for k in 0:12, j in 0:1//2:k//2+2
        @test qdim(Level(k), j) ≈ float(qdim(Exact(k), j)) atol = 1e-13
    end
    @test qdim(1; k = 2:4) == [qdim(Level(2), 1), qdim(Level(3), 1), qdim(Level(4), 1)] && qdim(Exact(2:3), 1) isa Vector{QForm}
end

@testset "matrix functions take targets" begin
    # the same matrices as the keywords, for each numerical target
    for (f, a) in ((fmatrix, (1, 1, 1, 1)), (bmatrix, (1, 1, 1, 1)), (qcg_matrix, (1, 1, 0)), (qcg_matrix, (1, 1)), (qcg_row, (1, 1, 0, 0)))
        for (t, kw) in ((Level(10), (; k = 10)), (At(0.8), (; q = 0.8)), (At(0.8 + 0.2im), (; q = 0.8 + 0.2im)), (Classical(), (; q = 1)))
            x = f(t, a...); y = f(a...; kw...)
            @test first(x) == first(y) && x[2:end] == y[2:end]
        end
        @test first(f(Level(10), a...; T = ComplexF64)) == first(f(a...; k = 10, T = ComplexF64))
        for t in (Exact(), Exact(10), Symbolic(), Classical(exact = true))
            @test_throws ArgumentError f(t, a...)
        end
        @test_throws ArgumentError f(Level(10), a...; k = 3)
        @test_throws ArgumentError f(At(0.8), a...; q = 0.7)
    end
    @test eltype(first(fmatrix(Level(10; T = BigFloat), 1, 1, 1, 1))) === BigFloat
    @test eltype(first(fmatrix(Level(10; T = Float32), 1, 1, 1, 1))) === Float32
    @test fmatrix_labels(Level(4), 1, 1, 1, 1) == fmatrix_labels(1, 1, 1, 1; k = 4) && fmatrix_labels(Classical(), 1, 1, 1, 1) == fmatrix_labels(1, 1, 1, 1)
    @test fmatrix_labels(At(0.8), 1, 1, 1, 1) == fmatrix_labels(1, 1, 1, 1)
    @test_throws ArgumentError fmatrix_labels(Exact(4), 1, 1, 1, 1)
    @test smatrix(Level(10)) == smatrix(10) && tmatrix(Level(10)) == tmatrix(10)
    @test eltype(first(smatrix(Level(10; T = Float32)))) === Float32 && tmatrix(Level(10); anomaly = false) == tmatrix(10; anomaly = false)
    @test_throws ArgumentError smatrix(At(0.8))
    @test_throws ArgumentError tmatrix(Exact(10))
end

@testset "forms do not mix with values at q = 1" begin
    for f in (q6j(Exact(3), 1, 1, 1, 1, 1, 1), q6j(Symbolic(), 1, 1, 1, 1, 1, 1)), g in (+, -, *, /)
        @test_throws ArgumentError g(f, q6j(Exact(), 1, 1, 1, 1, 1, 1))
        @test_throws ArgumentError g(q6j(Exact(), 1, 1, 1, 1, 1, 1), f)
    end
    # the multiword type keeps its arithmetic with plain numbers, and has no ambiguous methods left
    a = QR.MW{3}(0.25)
    @test Float64(a + 1) == 1.25 && Float64(1.0 - a) == 0.75 && Float64(2a) == 0.5 && Float64(a / 2.0) == 0.125 && Float64(a * true) == 0.25 && a < 1 && 0 < a && a <= 0.25
    @test isempty(Test.detect_ambiguities(QRecoupling; recursive = true))
end

@testset "tiny values at real q" begin
    # a true value below floatmin rounds to zero in Float64, which is not an exact zero, and BigFloat shows it
    l = (300, 300, 300, 300, 300, 299)
    @test q6j(At(0.1), l...) == 0
    @test !iszero(q6j(At(0.1), l...; T = BigFloat))
    @test !iszero(q6j(Symbolic(), l...))
end

@testset "equality at a level, read at places" begin
    # against the dense kernel, on pairs close enough to be equal, with both signs, and with numbers
    for (k, top) in ((4, 2), (10, 5//2), (16, 5//2), (18, 5//2))
        r = 0:1//2:top
        fs = [q6j(Exact(k), a, b, c, d, e, f) for a in r, b in r, c in r, d in r, e in r, f in r
              if QR._admissible_at(q6j, k, a, b, c, d, e, f)]
        fs = filter(!iszero, fs)
        close = Dict{Int,Vector{eltype(fs)}}()
        for f in fs
            push!(get!(close, round(Int, abs(float(f)) * 10^9), eltype(fs)[]), f)
        end
        bad = 0
        for b in values(close), i in eachindex(b), j in i:min(i + 3, length(b))
            bad += (b[i] == b[j]) != (QR._form_value(b[i]) == QR._form_value(b[j]))
            bad += (b[i] == -b[j]) != (QR._form_value(b[i]) == -QR._form_value(b[j]))
        end
        for f in fs[1:7:end], c in (1, -1, 1//2)
            bad += (f == c) != (QR._form_value(f) == c)
        end
        @test bad == 0
    end
    @test q6j(Exact(3), 1, 1, 1, 1, 1, 1) != 1 && 1 == q6j(Exact(10), 0, 0, 0, 0, 0, 0)
    @test iszero(q6j(Exact(4), 1, 1, 1, 5, 5, 5)) && q6j(Exact(4), 1, 1, 1, 5, 5, 5) == 0

    # the proof with powers of q: q^(2h) = 1 gives the same square and value, q^h = −1 the opposite value
    s, (w, e2) = QR.qcg_rule(3, 1, 2, 0, 3, 1)
    k = 9; h = k + 2
    V = QR.LevelSquare(s, k, (w, e2))
    @test QR._level_squares_equal(V, QR.LevelSquare(s, k, (w + 2h, e2 + 4h)), k)
    @test QR._level_squares_equal(V, QR.LevelSquare(s, k, (w, e2 + 2h)), k)
    @test !QR._level_squares_equal(V, QR.LevelSquare(s, k, (w + 1, e2)), k)
    c = QR.qform(s; k = k, weight = (w, e2))
    @test QR.qform(s; k = k, weight = (w, e2 + 4h)) == c
    @test QR.qform(s; k = k, weight = (w, e2 + 2h)) == -c
    @test QR.qform(s; k = k, weight = (w, e2 + 2h)) != c

    # the sizes behind the proof hold at every conjugate, where cancellation can only help
    for (k, j) in ((10, 1), (20, 4), (41, 10))
        s = QR._rule_for(q6j, j, j, j, j, j, j)
        K, tops = QR._square_sizes(QR.LevelSquare(s, k), k)
        zs, pre, r, es, cs = QR._level_exponents(s, k, Rational{BigInt})
        h = k + 2
        setprecision(BigFloat, 512) do
            for (i, a) in enumerate(a for a in 1:h-1 if gcd(a, 2h) == 1)
                qn(n) = sinpi(big(a * n) / h) / sinpi(big(a) / h)
                mono(e) = prod(qn(n)^x for (n, x) in enumerate(e) if !iszero(x); init = big(1))
                square = (sum(big(numerator(c)) / denominator(c) * mono(e) for (c, e) in zip(cs, es)))^2 *
                         big(numerator(r)) / denominator(r) * mono(pre)
                iszero(square) || @test log2(abs(square)) + K <= tops[i] / 2^20
            end
        end
    end
end

@testset "exact values are forms" begin
    # one target with no q: the form alone, and the value at q = 1 when asked for
    f = q6j(Symbolic(), 1, 1, 1, 1, 1, 1)
    @test f isa QForm && f.k === nothing
    @test string(f) == "([5] − 1)/[4]!"
    @test sprint(show, MIME("text/plain"), f) == "([5] − 1)/[4]!"
    @test f(1) isa QR.ClassicalResult && f(1) == 1//6 && f(1//1) == 1//6 && f(1) == q6j(Exact(), 1, 1, 1, 1, 1, 1) && Rational(f(1)) == 1//6
    @test f(1.0) isa Float64 && f(1.0) ≈ 1/6 && f(0.8) ≈ q6j(At(0.8), 1, 1, 1, 1, 1, 1)
    @test sprint(show, MIME("text/plain"), q3j(Symbolic(), 2, 2, 0, 1, -1, 0)) == "−q √(1/[5])"
    @test sprint(show, MIME("text/plain"), qdim(Symbolic(), 1//2)) == "[2]"
    @test (-f)(1) == -1//6 && q3j(Symbolic(), 2, 2, 0, 1, -1, 0)(1)^2 == 1//5 && q6j(Symbolic(), 1, 1, 1, 1, 1, 5)(1) == 0
    # it holds no number: what needs a value of q says so
    for g in (float, Float64, big, numeric_value, inv, x -> 2x, x -> x * x, x -> x + x, x -> x^2, x -> x ≈ 1/6, x -> convert(Float64, x), x -> x * f(1))
        @test_throws ArgumentError g(f)
    end
    # a number is a constant function of q, so `==` still means every q
    @test f != 1//6 && 1//6 != f && f != 1/6
    @test qint(Symbolic(), 1) == 1 && 1 == qint(Symbolic(), 1) && qint(Symbolic(), 1) == 1.0 && q6j(Symbolic(), 0, 0, 0, 0, 0, 0) == 1
    @test qint(Symbolic(), 6) != 6 && qint(Symbolic(), 6)(1) == 6 && qdim(Symbolic(), 1//2) != 2 && qfact(Symbolic(), 3) != 6
    @test qcg(Symbolic(), 0, 0, 1//2, 1//2, 1//2, 1//2) == 1 && -qint(Symbolic(), 1) == -1 && qint(Symbolic(), 1) != -1
    @test qcg(Symbolic(), 1//2, 1//2, 1//2, -1//2, 1) != 1 && qcg(Symbolic(), 1//2, 1//2, 1//2, -1//2, 1)(1)^2 == 1//2
    @test iszero(q6j(Symbolic(), 1, 1, 1, 1, 1, 5))
    # a classical zero is not a zero at every q, and a form says which it is
    z = first(l for l in all_6j(k = 16) if iszero(q6j(Exact(), l...)) && !iszero(q6j(At(0.8), l...)))
    @test !iszero(q6j(Symbolic(), z...)) && q6j(Symbolic(), z...) != 0 && q6j(Symbolic(), z...)(1) == 0

    # a level: the same form, simplified, and no backend needed for any of it
    g = q6j(Exact(3), 1, 1, 1, 1, 1, 1)
    @test g isa QForm && g.k == 3 && string(g) == "−1/[2]²   (k = 3)"
    @test sprint(show, MIME("text/plain"), g) == "−1/[2]²   (k = 3)"
    # arithmetic gives a value at the level, printed as a form is, and a sum of unlike roots in a fixed order
    x = q6j(Exact(10), 1, 1, 1, 1, 1, 1)
    @test string(x^2) == sprint(show, MIME("text/plain"), x^2) == "(−4[3] + 11)/3   (k = 10)"
    @test string(q6j(Exact(10), 3//2, 1, 1//2, 1, 3//2, 1) + q6j(Exact(10), 1, 1, 2, 1, 1, 1)) ==
          "(−[3] + 4)/6 + (−[3] + 4)/6 √([2])   (k = 10)"
    # roots and half-integer powers: exact where the root is again one, and an error that says why where it is not
    d = qdim(Exact(10), 1); p6 = q6j(Exact(10), 1, 1, 2, 1, 1, 1)
    @test string(sqrt(d)) == "√([3])   (k = 10)" && d^(1//2) == sqrt(d) && d^0.5 == sqrt(d) && d^(1/2) == sqrt(d)
    @test sqrt(d)^2 == d && d^(3//2) == d * sqrt(d) && d^(-1//2) * sqrt(d) == 1 && sqrt(p6)^2 == p6 && x^2.0 == x^2
    @test float(sqrt(p6)) ≈ sqrt(float(p6)) && sqrt(x - x + 9//4) == 3//2
    # a root of one product is a form, and so are products, quotients and powers of such forms: nothing is expanded
    R = sqrt(qint(Exact(1098), 550))
    @test sqrt(d) isa QForm && d * p6 isa QForm && string(R) == "√([550])   (k = 1098)" && string(R * R) == "[550]   (k = 1098)"
    @test R * R == qint(Exact(1098), 550) && float(R) ≈ sqrt(float(qint(Exact(1098), 550)))
    @test string(sqrt(d) * sqrt(qdim(Exact(10), 2))) == "√([3][5])   (k = 10)" && string(d^(-1//2)) == "√(1/[3])   (k = 10)"
    @test d / p6 == QR.level_value(d.rule, 10) / QR.level_value(p6.rule, 10) && inv(d) * d == 1 && d^0 == 1 && p6^-2 == inv(p6 * p6)
    @test string(sqrt(qint(Exact(10), 3) * qint(Exact(10), 9))) == "[3]   (k = 10)"   # [9] = [3] at h = 12
    G = sqrt(qdim(Symbolic(), 1))
    @test string(G) == "√([3])" && G(1)^2 == 3 && G(0.9) ≈ sqrt(qdim(At(0.9), 1)) && G * G == qdim(Symbolic(), 1)
    @test string(qdim(Symbolic(), 1) * qdim(Symbolic(), 2)) == "[3][5]" && inv(G)(0.9) ≈ 1 / G(0.9)
    for bad in (() -> sqrt(x), () -> sqrt(x^2), () -> sqrt(q6j(Exact(10), 3//2, 1, 1//2, 1, 3//2, 1)), () -> x^(1//3), () -> x^0.3,
              () -> sqrt(q6j(Symbolic(), 1, 1, 1, 1, 1, 1)), () -> sqrt(q3j(Exact(10), 1, 1, 1, 0, 0, 0)), () -> sqrt(sqrt(d)),
              () -> sqrt(G), () -> q6j(Symbolic(), 1, 1, 1, 1, 1, 1) * G, () -> G + G)
        @test_throws ArgumentError bad()
    end
    @test_throws DomainError sqrt(g)                      # g < 0, at level 3
    c = q6j(Exact(), 1, 1, 1, 1, 1, 1)
    @test sqrt(c^2) == c && sqrt(c) isa QR.ClassicalResult && sqrt(c)^2 == c && c^(1//2) == sqrt(c) && c^0.5 == sqrt(c)
    @test c^(-3//2) == inv(sqrt(c))^3
    @test_throws ArgumentError sqrt(sqrt(c))
    @test_throws ArgumentError c^(1//3)
    @test_throws DomainError sqrt(-c)
    # sums of unlike roots: integer powers, and no root or inverse
    s = q6j(Exact(10), 3//2, 1, 1//2, 1, 3//2, 1) + p6
    @test s^2 == s * s && s^(2//1) == s * s && s^0 == 1 && float(s^3) ≈ float(s)^3
    for bad in (() -> sqrt(s), () -> inv(s), () -> s^-1, () -> sqrt(c + sqrt(c)))
        @test_throws ArgumentError bad()
    end
    @test string(radical(g)) == "(√5 − 3)/2"
    @test_throws ArgumentError g(1)
    @test_throws ArgumentError f == g
    @test q6j(1, 1, 1, 1, 1, 1; k = 3, exact = true) == g
    @test iszero(q6j(Exact(10), 1, 3//2, 3//2, 2, 3//2, 3//2)) && iszero(q6j(Exact(4), 1, 1, 1, 5, 5, 5))
    @test [q6j(Exact(k), 1, 1, 1, 1, 1, 1) for k in 2:3] == q6j(Exact(2:3), 1, 1, 1, 1, 1, 1)
    @test length(q6j(Exact(10), [(1, 1, 1, 1, 1, 1), (2, 2, 2, 2, 2, 2)])) == 2

    # the keywords of `exact = true`: a level is a form, and q = 1 stays the classical number
    @test q6j(1, 1, 1, 1, 1, 1; exact = true) isa QR.ClassicalResult
    @test q6j(Classical(exact = true), 1, 1, 1, 1, 1, 1) == 1//6
    @test qcg(1, 0, 1, 0, 1; exact = true) == 0
    @test qcg(1//2, 1//2, 1//2, -1//2, 1; k = 4, exact = true) == qcg(Exact(4), 1//2, 1//2, 1//2, -1//2, 1)

    # every symbol, against numbers at a level and at generic q, and against Nemo
    for (sym, args) in ((q6j, (2, 2, 2, 2, 2, 2)), (fsymbol, (1, 1, 1, 1, 1, 1)), (gsymbol, (3//2, 1, 3//2, 1, 3//2, 1)),
                        (QR.tetrahedron, (1, 1, 1, 1, 1, 1)), (QR.theta_value, (1, 1, 1)), (qdim, (3//2,)),
                        (qint, (5,)), (qfact, (5,)), (qbinomial, (5, 2)), (QR.q3j_factorial, (1, 1, 1, 1, -1, 0)),
                        (q3j, (1, 1, 1, 1, -1, 0)), (qcg, (1, 0, 1, 0, 2)))
        for k in (5, 10, 16)
            h = sym(Exact(k), args...)
            @test h isa QForm && h.k == k
            @test float(h) ≈ sym(Level(k), args...) atol = 1e-12
            @test form_value(h) ≈ sym(Level(k), args...) atol = 1e-12
        end
        h = sym(Symbolic(), args...)
        @test h isa QForm && h.k === nothing
        @test h(1.3) ≈ sym(At(1.3), args...) rtol = 1e-12
        @test form_value(h; q = 1.3) ≈ sym(At(1.3), args...) rtol = 1e-12
        @test Float64(h(1)) ≈ sym(args...) rtol = 1e-12
    end
    # phases are not field values
    @test QR.rmatrix(Exact(10), 1//2, 1//2, 1) isa QPhase && QR.rmatrix(Symbolic(), 1//2, 1//2, 1) isa QPhase
    @test QR.twist(Exact(10), 1//2) isa QPhase
    # the expressions a level has no rule for

    # arithmetic at a level is exact, and it agrees with Nemo
    a = q6j(Exact(12), 1, 1, 1, 1, 1, 1); b = q6j(Exact(12), 3//2, 1, 3//2, 1, 3//2, 1)
    @test float(a * b) ≈ Float64(numeric_value(oracle(12, 1, 1, 1, 1, 1, 1) * oracle(12, 3//2, 1, 3//2, 1, 3//2, 1)))
    @test float(a + b) ≈ float(a) + float(b)
    @test (a + b) - b == a && a * inv(a) == 1 && a^2 * inv(a)^2 == 1
    @test zero(a) + a == a && 2a == a + a && a / 2 + a / 2 == a
    @test -a == (-1) * a
    @test sum(q6j(Exact(12), 1, 1, x, 1, 1, 1) * q6j(Exact(12), 1, 1, x, 1, 1, 1) * qdim(Exact(12), x) * qdim(Exact(12), 1) for x in 0:2) == 1
    # different forms of one value, at a level
    u = q6j(Exact(4), 0, 0, 0, 1//2, 1//2, 1//2); w = q6j(Exact(4), 0, 0, 0, 3//2, 3//2, 3//2)
    @test u == w
    @test q6j(Exact(10), 1, 1, 1, 1, 1, 1) != q6j(Exact(10), 2, 2, 2, 2, 2, 2)
    # arithmetic at q = 1 belongs to the value there
    @test qcg(Symbolic(), 1//2, 1//2, 1//2, -1//2, 0)(1)^2 == 1//2
    @test f(1) * f(1) == 1//36 && 6f(1) == 1

    # with powers of q: complex numbers, decided like any other value
    c = qcg(Exact(9), 3//2, 1//2, 1, 0, 3//2, 1//2)
    @test c == qcg(Exact(9), 3//2, 1//2, 1, 0, 3//2, 1//2)
    @test c != qcg(Exact(9), 3//2, 1//2, 1, 0, 1//2, 1//2)
    @test c != -c && -c == -c && c == -(-c)
    @test_throws ArgumentError c * c
    @test_throws ArgumentError radical(c)

    # equality at every q: different at q = 1 or not, and true across forms that differ
    @test qint(Symbolic(), 6) != qfact(Symbolic(), 3)                   # both 6 at q = 1
    @test qint(Symbolic(), 6)(1) == qfact(Symbolic(), 3)(1)
    @test qdim(Symbolic(), 1//2) != qdim(Symbolic(), 1)
    @test q6j(Symbolic(), 1, 1, 1, 1, 1, 1) != q6j(Symbolic(), 2, 2, 2, 2, 2, 2)
    @test qdim(Symbolic(), 1//2) != q6j(Symbolic(), 1//2, 1//2, 1, 1//2, 1//2, 1)
    # the Regge symmetry at every q
    regge(a, b, c, d, e, f) = (s = (b + c + e + f) // 2; (a, s - b, s - c, d, s - e, s - f))
    done = 0
    for l in all_6j(k = 12)
        l2 = regge(l...)
        (all(x -> denominator(2x) == 1 && x >= 0, l2) && l2 != l && length(QR._written(QR.qform(q6j, l...)).terms) > 1) || continue
        @test q6j(Symbolic(), l...) == q6j(Symbolic(), l2...)
        (done += 1) == 20 && break
    end
    @test done == 20
    # forms that share nothing: q-Vandermonde, a sum of 4 terms against one binomial, with a proof and not a match
    sum_form(w, e2) = QR.qform(FactorialSum(0:3; prefactor = (3 => 1, 4 => 1), factors = ((1, 0, -2), (-1, 3, -1), (-1, 4, -1))); weight = (w, e2))
    @test length(QR._written(sum_form(7, -24)).terms) == 4 && length(QR._written(qbinomial(Symbolic(), 7, 3)).terms) == 1
    @test sum_form(7, -24)(1.3) ≈ qbinomial(At(1.3), 7, 3)
    @test sum_form(7, -24) == qbinomial(Symbolic(), 7, 3)
    @test sum_form(7, -22)(1) == qbinomial(Symbolic(), 7, 3)(1) && sum_form(7, -22) != qbinomial(Symbolic(), 7, 3)
    @test sum_form(5, -24) != qbinomial(Symbolic(), 7, 3)
    # the proof itself: [2][3] = [4] + [2], q[2] = q² + 1, and near misses
    T(c, b, p...) = (Rational{BigInt}(c), b, collect(Pair{Int,Int}, p))
    vanishes = QR._generic_sum_vanishes
    @test vanishes([T(1, 0, 2 => 1, 3 => 1), T(-1, 0, 4 => 1), T(-1, 0, 2 => 1)])
    @test !vanishes([T(1, 0, 2 => 1, 3 => 1), T(-1, 0, 4 => 1)])
    @test vanishes([T(1, 1, 2 => 1), T(-1, 2), T(-1, 0)])
    @test !vanishes([T(1, 1, 2 => 1), T(-1, 2), T(-1, 1)])
    @test vanishes([T(1//3, 0, 2 => 2, 5 => -1), T(-1//3, 0, 3 => 1, 5 => -1), T(-1//3, 0, 5 => -1)])  # [2]² = [3] + 1
    @test !vanishes([T(1//3, 0, 2 => 2, 5 => -1), T(-1//3, 0, 3 => 1, 5 => -1)])
    @test !vanishes([T(big(10)^40, 0, 2 => 1), T(-big(10)^40 + 1, 0, 2 => 1)])
end

@testset "optional backend product hooks" begin
    @test Float64(qdim(Exact(10),1)) ≈ qdim(Level(10),1)
    @test Float64(qfact(5;k=10,exact=true)) ≈ qfact(5;k=10)
    @test Float64(qbinomial(5,2;k=10,exact=true)) ≈ qbinomial(5,2;k=10)
end

@testset "real-q F-matrices past the exponent range" begin
    # The labels of LorentzianRecoupling's su2_dual_labels(1//2, 1, 3//2, N). At real q the column's
    # coefficients grow like (1/q)^(spin): at q = 0.5 the plain recurrence gave NaN from N ≈ 100 (every
    # entry at N = 150), at q = 0.7 from N ≈ 200; at q = 0.1, N = 300 some coefficients pass 2^1024.
    dual(N) = (1//2 + (N - 1)//2, 1 + (N - 1)//2, 3//2 + (N - 1)//2, 3 + (N - 3)//2)
    for (q, N) in ((0.5, 150), (0.7, 200), (0.1, 300))
        a, b, c, d = dual(N)
        F, es, fs = fmatrix(a, b, c, d; q)
        @test all(isfinite, F)
        @test transpose(F) * F ≈ [i == j for i in eachindex(fs), j in eachindex(fs)] atol = 1e-14
        for (i, j) in ((1, 1), (N ÷ 2, N ÷ 2 + 1), (N + 1, N + 1), (N ÷ 3, N ÷ 3 + 3), (N ÷ 3, 2N ÷ 3))
            ref = setprecision(() -> fsymbol(a, b, es[i], c, d, fs[j]; q = big(q), T = BigFloat), BigFloat, 256)
            @test F[i, j] ≈ ref rtol = 1e-14 atol = 1e-300
        end
    end

    # The scaled column is taken only for large labels; where the plain one also runs, they agree bit for bit.
    a, b, c, d = dual(40)
    A, B, C, D = QR.doubled(a, b, c, d)
    fs = QR._fmatrix_cols(A, B, C, D, nothing)
    es = QR.sixj_column_range(B, A, first(fs), D, C, nothing)
    Q = QR.real_q_tables(0.5, B, A, first(fs), D, C)
    @test QR._column_q(Q, es, B, A, D, C) isa QR.ScaledRealQ
    @test QR._column_q(QR.real_q_tables(0.95, B, A, first(fs), D, C), es, B, A, D, C) isa QR.RealQ
    plain = zeros(length(es), length(fs))
    QR._fmatrix_fill!(plain, Q, es, fs, A, B, C, D, a, b, c, d, nothing, 0.5)
    @test plain == first(fmatrix(a, b, c, d; q = 0.5))

    # The bare 6j column (no dimension factor) through the same path, at q > 1.
    J2, J3, L1, L2, L3 = 151, 150, 201, 153, 152
    X = QR.sixj_column_range(J2, J3, L1, L2, L3, nothing)
    col = zeros(length(X))
    QR.sixj_column!(col, J2, J3, L1, L2, L3, QR.real_q_tables(2.0, J2, J3, L1, L2, L3))
    @test all(isfinite, col)
    for i in (1, length(X) ÷ 2, length(X))
        ref = q6j(X[i]//2, J2//2, J3//2, L1//2, L2//2, L3//2; q = 2.0)
        @test col[i] ≈ ref rtol = 1e-13 atol = 1e-300
    end
end

@testset "exact values at q = 1 are real numbers to generic code" begin
    # `plot([q6j(Exact(), i*ones(6)...) for i in 1:50])` failed in `isinf`, whose `Real` fallback needs `decompose`.
    v = [q6j(Exact(), i, i, i, i, i, i) for i in 1:20]
    out = similar(v, Float64)
    broadcast!(x -> isinf(x) ? NaN : x, out, v)              # what Plots does with series data
    @test out == Float64.(v) && all(isfinite, v) && !any(isnan, v)

    x = q6j(Exact(), 1, 1, 1, 1, 1, 1)                      # 1//6
    r = QR.ClassicalResult(-1, 2//3)                         # −√(2/3)
    s = QR.ClassicalResult(1, 2//1) + QR.ClassicalResult(1, 3//1)   # √2 + √3
    @test s isa QR.ClassicalSum
    # exact order among exact values and rationals, rounded to Float64 against floats as `==` is
    @test r < x < 1//5 < s < 4 && x <= x && !(x < x) && -s < r
    @test sort([s, x, r, -s, 0]) == [-s, r, 0, x, s] && extrema([x, s, r]) == (r, s)
    @test x < 0.2 && !(x < 1/6) && x == 1/6 && isless(x, NaN) && !isless(NaN, x)
    @test max(x, 1) == 1 && max(x, 1) isa QR.ClassicalResult && min(s, x) === x && max(s, 3) === s
    @test max(x, 0.1) isa BigFloat && x + π ≈ 1/6 + π && promote(x, 1)[2] == 1
    @test [x, 1] isa Vector{QR.ClassicalResult}
    @test signbit(r) && !signbit(x) && !signbit(zero(x))
    # rounding is exact, ties to even only where the value is rational
    @test round(QR.ClassicalResult(1, 25//4)) == 2 && round(QR.ClassicalResult(1, 49//4)) == 4
    @test floor(r) == -1 && ceil(r) == 0 && trunc(Int, r) == 0 && round(Int, s) == 3 && floor(Int, -s) == -4
    @test isinteger(q6j(Exact(), 0, 0, 0, 0, 0, 0)) && !isinteger(x) && !isinteger(s)
    big_root = QR.ClassicalResult(1, big(10)^40 + 1)         # √(10⁴⁰+1), just above 10²⁰
    @test floor(BigInt, big_root) == big(10)^20 && ceil(BigInt, big_root) == big(10)^20 + 1
end

@testset "exact values at a level are ordered real numbers" begin
    v = [q6j(Exact(40), i, i, i, i, i, i) for i in 0:8]
    out = similar(v, Float64)
    broadcast!(x -> isinf(x) ? NaN : x, out, v)              # what Plots does with series data
    @test out == Float64.(v) && float(zero(eltype(v))) === 0.0
    @test sort(v) == v[sortperm(Float64.(v))] && maximum(v) == v[argmax(Float64.(v))]

    a = q6j(Exact(10), 1, 1, 1, 1, 1, 1)                     # 0.1547…
    b = q6j(Exact(10), 1, 1, 1, 1, 1, 2)                     # 0.2113…
    c = QR._form_value(a)                                     # the same value as a LevelValue
    d = c + 1 // big(10)^30                                   # differs from it far below Float64's resolution
    @test a < b && !(a < c) && !(c < a) && a <= c && c < d && !(d < c) && sign(d - c) == 1
    @test 0.154 < a < 0.155 && a < 1 && -1 < -a && !(a < float(c)) == (float(c) <= c)     # exact against floats
    @test isless(a, NaN) && !isless(NaN, a) && sign(-a) == -1 && abs(-a) == a && signbit(-a)
    @test max(2c, 1) == 1 && min(2c, 1) == 2c && max(1//2, 2c) == 1//2 && conj(a) == a && iszero(imag(a))
    @test a + 1.0 ≈ float(a) + 1 && a * π ≈ π * float(a) && convert(Float64, c) == float(c) && big(c) ≈ float(c)
    # complex (powers of q) and generic forms have no order, and say so
    @test_throws ArgumentError q3j(Exact(10), 1, 1, 1, 0, 0, 0) < 1
    @test_throws ArgumentError conj(q3j(Exact(10), 1, 1, 1, 0, 0, 0))
    @test_throws ArgumentError q6j(Symbolic(), 1, 1, 1, 1, 1, 1) < 1
end

@testset "matrices of exact values, and hashing against floats" begin
    # a matrix product starts from zero(T), the zero that holds at every q, which must add at any level
    e, f = fmatrix_labels(1, 1, 1, 1; k = 6)
    F = [fsymbol(Exact(6), 1, 1, ei, 1, 1, fj) for ei in e, fj in f]
    @test transpose(F) * F == [i == j for i in axes(F, 2), j in axes(F, 2)]       # orthogonal, exactly
    a = q6j(Exact(10), 1, 1, 1, 1, 1, 1); b = q6j(Exact(10), 1, 1, 1, 1, 1, 2)
    @test zero(QForm) + a == a && a - zero(QForm) == a && zero(QForm) - a == -a && iszero(zero(QForm) * a)
    @test sum([a, b]; init = zero(QForm)) == a + b && [a, QR._form_value(b)] isa Vector{QR.LevelValue}
    @test_throws ArgumentError a + q6j(Exact(12), 1, 1, 1, 1, 1, 1)               # levels still do not mix
    # equal values hash alike: `isequal` is exact against a float, though `==` rounds there
    x = q6j(Exact(), 1, 1, 1, 1, 1, 1); h = QR.ClassicalResult(1, 1//4)            # 1//6 and exactly 1//2
    @test x == 1/6 && !isequal(x, 1/6) && length(Set([x, 1/6])) == 2
    @test isequal(h, 0.5) && isequal(0.5, h) && hash(h) == hash(0.5) && length(Set([h, 0.5])) == 1
    @test isequal(zero(x), 0.0) && !isequal(zero(x), -0.0)
end
