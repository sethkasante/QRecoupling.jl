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
        p = Base.get_extension(QR, :QRecouplingNemoExt).phi_q(n)
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

# `float(f)` and `f(q)` evaluate the rule behind a form, so they cannot see a mistake in the form. These two read the
# form itself: its fields, and the text a reader gets.
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
    f = qform(q6j, 1, 1, 1, 2, 1, 1; k = 5)
    wrong = typeof(f)(f.k, f.sign, 7 * f.coef, f.out, f.radcoef, f.rad, f.terms, f.rule, f.weight, f.qpow, f.powers)
    @test float(wrong) ≈ float(f)
    @test form_value(wrong) ≈ 7 * float(f)
    @test text_value(string(wrong); k = 5) ≈ 7 * float(f)
    # at levels, prime and composite, with and without the relations, including roots and factorials
    for k in (3, 5, 10, 12, 16, 21, 30), l in collect(all_6j(k = min(k, 12)))[1:97:end]
        QR._admissible_at(q6j, k, l...) || continue
        reference = q6j(Level(k), l...)
        for reduce in (true, false)
            g = qform(q6j, l...; k = k, reduce = reduce)
            iszero(g) && continue
            @test form_value(g) ≈ reference atol = 1e-12
            length(g.terms) <= QR.QFORM_SHOWN_TERMS && @test text_value(string(g); k = k) ≈ reference atol = 1e-12
        end
    end
    # at generic q, q = 1 included
    for l in collect(all_6j(k = 10))[1:211:end], q in (0.7, 1, 1.3)
        g = qform(q6j, l...)
        iszero(g) && continue
        reference = isone(q) ? q6j(l...) : q6j(At(q), l...)
        @test form_value(g; q = q) ≈ reference rtol = 1e-10
        length(g.terms) <= QR.QFORM_SHOWN_TERMS && @test text_value(string(g); q = q) ≈ reference rtol = 1e-10
    end
end

@testset "forms with powers of q" begin
    # a wrong power of q has to show, in the form and in its text
    f = qform(qcg, 3//2, 1//2, 1, 0, 3//2, 1//2; k = 9)
    wrong = typeof(f)(f.k, f.sign, f.coef, f.out, f.radcoef, f.rad, f.terms, f.rule, f.weight, f.qpow + 2, f.powers)
    @test float(wrong) ≈ float(f)
    @test !(form_value(wrong) ≈ form_value(f))
    @test !(text_value(string(wrong); k = 9) ≈ form_value(f))
    labels = [(a, m1, b, m2, j, m1 + m2) for a in (0:6) .// 2 for b in (0:6) .// 2 for j in abs(a - b):1:a+b
              for m1 in -a:1:a for m2 in -b:1:b if abs(m1 + m2) <= j]
    for l in labels[1:7:end]
        g = qform(qcg, l...)
        @test !iszero(g)
        for q in (0.7, 1.3)
            @test form_value(g; q = q) ≈ qcg(l...; q = q) rtol = 1e-11
            @test g(q) ≈ qcg(l...; q = q) rtol = 1e-11
            length(g.terms) <= QR.QFORM_SHOWN_TERMS && @test text_value(string(g); q = q) ≈ qcg(l...; q = q) rtol = 1e-11
        end
        @test form_value(g; q = 1) ≈ qcg(l...) rtol = 1e-12
        @test g(1) ≈ qcg(l...) rtol = 1e-12
        for k in (1, 2, 3, 4, 6, 10, 15)
            v = qcg(l...; k = k)
            for reduce in (true, false)
                h = qform(qcg, l...; k = k, reduce = reduce)
                @test iszero(h) == iszero(v)
                iszero(h) && continue
                @test form_value(h) ≈ v atol = 1e-11
                @test ComplexF64(numeric_value(h)) ≈ v atol = 1e-11
                @test float(h) ≈ v atol = 1e-11
                length(h.terms) <= QR.QFORM_SHOWN_TERMS && @test text_value(string(h); k = k) ≈ v atol = 1e-11
            end
        end
    end
    for (a, b, c, m1, m2) in ((1, 1, 1, 1, -1), (2, 2, 2, 0, 0), (3//2, 1, 1//2, 1//2, 0), (5//2, 2, 3//2, -1//2, 1))
        m3 = -m1 - m2
        @test form_value(qform(q3j, a, b, c, m1, m2, m3); q = 0.8) ≈ q3j(a, b, c, m1, m2, m3; q = 0.8) rtol = 1e-11
        k = 2Int(a + b + c)
        @test form_value(qform(q3j, a, b, c, m1, m2, m3; k = k)) ≈ q3j(a, b, c, m1, m2, m3; k = k) atol = 1e-11
    end
    # the printed powers: a form with one term, at generic q and reduced by q^h = −1 at a level
    @test string(qform(qcg, 1//2, 1//2, 1//2, -1//2, 1, 0)) == "q^(−1/2) √(1/[2])"
    @test string(qform(qcg, 1//2, 1//2, 1//2, -1//2, 0, 0)) == "q^(1/2) √(1/[2])"
    @test endswith(string(qform(qcg, 1//2, 1//2, 1//2, -1//2, 1, 0; k = 3)), "(k = 3, q = e^{iπ/5})")
    # outside the fusion rule
    @test iszero(qform(qcg, 1, 0, 1, 0, 2, 0; k = 3))
    @test iszero(numeric_value(qform(qcg, 1, 0, 1, 0, 2, 0; k = 3)))
    # terms with one monomial and different powers do not merge at generic q, and do merge when q^h = −1 aligns them
    h = 9; s = FactorialSum(0:2; factors = [(0, 2, 1)], alternating = false)
    @test length(qform(s; weight = (2h ÷ 3, 0)).terms) == 3
    # 1 + q^(2h/3) + q^(4h/3) is a sum of cube roots of unity: no two terms merge, and it vanishes
    @test length(qform(s; k = h - 2, reduce = false, weight = (2h ÷ 3, 0)).terms) == 3
    @test iszero(qform(s; k = h - 2, weight = (2h ÷ 3, 0)))
    @test !iszero(qform(s; k = h - 2, weight = (2h ÷ 3 + 2, 0)))
    @test iszero(qform(FactorialSum(0:1; factors = [(0, 2, 1)], alternating = true); k = h - 2, weight = (2h, 0)))
    @test length(qform(FactorialSum(0:1; factors = [(0, 2, 1)], alternating = true); k = h - 2, weight = (h, 0)).terms) == 1
    # exact values with weights are not a target yet
    @test_throws ArgumentError qcg(Exact(5), 1, 0, 1, 0, 1, 0)
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
    @test setprecision(() -> QR._qform_level(qform(q6j, 1, 1, 1, 1, 1, 1; k = 5)), BigFloat, 288) ≈
          numeric_value(q6j(Exact(5), 1, 1, 1, 1, 1, 1); bits = 288) rtol = big(2.0)^-256

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
    for h in (12, 30), (_, row, con, _) in QR._relations_at(h)
        @test prod(qn(n, h)^row[n] for n in 2:h÷2) ≈ Float64(con)
    end
    # generators are distinct, and they span the relations: the rank is what the units leave free
    for h in (12, 30, 60, 77, 121, 210, 420, 714, 1002)
        @test length(QR._level_relations(h)) == (h ÷ 2 - 1) - (count(n -> gcd(n, h) == 1, 1:h) ÷ 2 - 1)
    end
    @test isempty(QR._level_relations(101))
    # very smooth levels: finite in time and memory, every row true, and a form still comes back right
    logq(n, h) = log(sinpi(big(n) / h)) - log(sinpi(1 / big(h)))
    for h in (2002, 2310)
        t = @elapsed rows = QR._level_relations(h)
        @test t < 60 && !isempty(rows)
        @test all(r -> maximum(abs, r[2]) <= QR.RELATION_EXPONENT_LIMIT, rows)
        for (piv, row, con, _) in rows[1:29:end]
            @test all(iszero, row[piv+1:end]) && row[piv] > 0
            @test abs(sum(row[n] * logq(n, h) for n in 2:h÷2 if !iszero(row[n])) - log(BigFloat(con))) < 1e-30
        end
        # the classes stay complete without the table, with small entries
        F = QR._level_field(h)
        cl = QR._class_rows(F)
        @test count(r -> r.pivot > 0, cl) == length(unique(r.pivot for r in cl if r.pivot > 0)) > length(rows)
        @test maximum(maximum(abs, r.e) for r in cl) < 10^4
    end
    v = QR.level_value(QR._rule_for(q6j, 3, 3, 3, 3, 3, 3), 2000)
    @test float(v) ≈ q6j(Level(2000), 3, 3, 3, 3, 3, 3) rtol = 1e-12
    @test float(v * v) ≈ float(v)^2 rtol = 1e-12
    @test float(qform(q6j, 3, 3, 3, 3, 3, 3; k = 2000)) ≈ q6j(Level(2000), 3, 3, 3, 3, 3, 3) rtol = 1e-12
    # the largest levels go without a table, and a form is still correct
    @test isempty(QR._level_relations(QR.RELATION_TABLE_LIMIT + 2))
    # an exponent past the range of Int stops a reduction, and never wraps
    @test_throws OverflowError QR._reduce_exponents!([0, typemin(Int) + 10, 100], [(3, [0, 4, 1], 1 // big(1), [2, 3])])

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
        @test length(v.num) == QR._totient(2(k + 2)) ÷ 2
    end
    # [1], …, [d] carry the arithmetic: at k = 10 the value (2√3 − 3)/3 is (2[3] − 5)/3
    @test QR.coordinates(QR.level_value(QR._rule_for(q6j, ntuple(_ -> 1, 6)...), 10)) == [-5 // 3, 0 // 1, 2 // 3, 0 // 1]
    # zeros that only additive relations show, proved without the exact field
    for (k, l) in ((20, (5//2, 9//2, 5, 4, 5, 11//2)), (22, (2, 3, 3, 5//2, 7//2, 7//2)), (10, (1, 5//2, 5//2, 5//2, 3, 3)))
        @test iszero(qform(q6j, l...; k = k))
        @test iszero(QR.level_value(QR._rule_for(q6j, l...), k))
        @test iszero_at(k, l...)
    end
    @test !iszero(qform(q6j, 2, 3, 3, 5//2, 7//2, 7//2; k = 23))
    # arithmetic agrees with the exact field, with and without a shared root
    let k = 12
        a = q6j(Exact(k), 1, 1, 1, 1, 1, 1); b = q6j(Exact(k), 3//2, 1, 3//2, 1, 3//2, 1)
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
        @test form_value(qform(u.rule; k = k)) ≈ Float64(numeric_value(u))
        shown = replace(string(u), "ExactX(k = $k, " => "")[1:end-1]
        occursin("terms)", shown) || @test text_value(shown; k = k) ≈ Float64(numeric_value(u))
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
    @test !occursin(r"\[(\d+)\]![⁰¹²³⁴⁵⁶⁷⁸⁹]*\[\1\]!", string(qform(q6j, ntuple(_ -> 30, 6)...; k = 200)))
    @test occursin(r"\[(\d+)\]![⁰¹²³⁴⁵⁶⁷⁸⁹]*\[\1\]!", "[15]!²[15]!")
    @test string(q6j(Exact(10), 1, 5//2, 5//2, 5//2, 3, 3)) == "ExactX(k = 10, 0)"
    # equality does not look at the display
    @test v == q6j(Exact(5), 1, 1, 1, 2, 1, 1) && v.rule !== nothing
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
    f = qform(q6j, 2, 2, 2, 2, 2, 2)
    q = 1 + 1e-8
    ref = setprecision(() -> f(BigFloat(q)), BigFloat, 256)
    @test abs((f(q) - ref) / ref) < 1e-12
    # The enclosure must include the value, not merely agree to a tolerance.
    v = QR.level_value(QR._rule_for(q6j, 2, 2, 2, 2, 2, 2), 12)
    mid, err = setprecision(() -> QR._value_ball(v, 14), BigFloat, 128)
    exactref = numeric_value(q6j(Exact(12), 2, 2, 2, 2, 2, 2); bits = 512)
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

@testset "optional backend product hooks" begin
    @test Float64(qdim(Exact(10),1)) ≈ qdim(Level(10),1)
    @test Float64(qfact(5;k=10,exact=true)) ≈ qfact(5;k=10)
    @test Float64(qbinomial(5,2;k=10,exact=true)) ≈ qbinomial(5,2;k=10)
end
