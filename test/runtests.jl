using Test, QRecoupling, Nemo

const QR = QRecoupling
const ONES = (1, 1, 1, 1, 1, 1)

@testset "QRecoupling.jl" begin
    @testset "targets and closed forms" begin
        @test q6j(ONES...) ≈ 1/6
        @test q6j(Exact(), ONES...) == 1//6
        @test q6j(0,0,0,0,0,0) == 1
        @test q6j(1,1,1,1,1,4) == 0
        @test q6j(Level(1), ONES...) == 0

        # {1 1 1; 1 1 1} = (x²-3)/((x²-1)(x²-2)), x=q+1/q.
        formula(x) = (x^2-3)/((x^2-1)*(x^2-2))
        for q in (0.8, 0.7+0.4im)
            @test q6j(At(q), ONES...) ≈ formula(q+inv(q)) rtol=1e-13
        end
        v = q6j(Exact(10), ONES...)
        @test v isa QForm
        @test Float64(v) ≈ formula(2cospi(1/12)) rtol=1e-13
        @test q6j(Level(10), ONES...) ≈ Float64(v) rtol=1e-13
        @test string(radical(v)) == "(2√3 − 3)/3"
        @test iszero(v-v) && isone(v*inv(v))

        s = q6j(Symbolic(), ONES...)
        @test s isa QForm && s.k === nothing && x_form(s) isa XValue
        @test s(1) == 1//6 && qeval(Exact(), QRecoupling._rule_for(q6j, ONES...)) == 1//6
        @test qeval(s; q=0.8) ≈ formula(0.8+inv(0.8)) rtol=1e-13
    end

    @testset "products and factorial rules" begin
        @test qint(0) == qfact(0) == 1 # the package's qint(0) convention
        @test qfact(5) == 120 && qbinomial(5,2) == 10
        @test qint(5; q=0.8) ≈ (0.8^5-0.8^-5)/(0.8-0.8^-1)
        @test qint(25; k=10) ≈ sinpi(25/12)/sinpi(1/12)
        @test qdim(Exact(), 1//2) == 2
        # Sum of squared binomial coefficients: sum_z binomial(4,z)^2 = 70.
        s = FactorialSum(0:4; prefactor=(4=>2,), factors=((1,0,-2),(-1,4,-2)))
        @test qeval(Exact(), s) == 70
        @test qeval(s; k=12) ≈ Float64(qeval(Exact(12), s)) rtol=1e-13
    end

    @testset "fusion and repeated evaluation" begin
        for kw in ((; k=6), (; q=0.8+0.3im))
            F, es, fs = fmatrix(1,1,1,1; kw...)
            identity = [i==j for i in eachindex(fs), j in eachindex(fs)]
            @test transpose(F)*F ≈ identity atol=1e-13
            @test F ≈ [fsymbol(1,1,e,1,1,f; kw...) for e in es, f in fs] atol=1e-13
        end
        @test gsymbol(ONES...; k=5) ≈ 1 atol=1e-13
        @test abs(rmatrix(Level(5),1,1,1)) ≈ 1
        @test abs(twist(Level(5),1//2)) ≈ 1
        work = EvaluationWorkspace()
        q6j(ONES...; q=0.8,workspace=work)
        @test q6j(ONES...; q=0.8,workspace=work) == q6j(ONES...; q=0.8)
        labels = [ONES, (2,2,2,2,2,2)]
        @test q6j(Level(10), labels) ≈ [q6j(l...; k=10) for l in labels]
        @test size(q6j(Level(9:10), labels)) == (2,2)
    end

    @testset "input errors" begin
        @test_throws ArgumentError q6j(1//3,1,1,1,1,1)
        @test_throws DomainError q6j(ONES...; k=-1)
        @test_throws ArgumentError q6j(Level(10), ONES...; q=0.8)
        @test_throws ArgumentError fmatrix(1,1,1,1; k=3.0)
    end

    @testset "tables without a wide number per entry" begin
        QR = QRecoupling
        u2 = 2.0^-106
        for k in (62, 63, 64, 100, 1000, 5000)
            tab = QR.QIntTables(Float64, k)
            h = k + 2; N = k + 1
            setprecision(BigFloat, 400) do
                S = sin(big(pi) / h); f = big(1.0)
                qerr = ferr = big(0.0)
                for n in 1:N
                    v = sin(n * big(pi) / h) / S
                    f *= v
                    qerr = max(qerr, abs((big(tab.q[n]) + tab.ql[n]) / v - 1), abs((big(tab.qi[n]) + tab.qil[n]) * v - 1))
                    ferr = max(ferr, abs(ldexp(big(tab.fm[n+1]) + tab.fml[n+1], tab.fe[n+1]) / f - 1),
                               abs(ldexp(big(tab.gm[n+1]) + tab.gml[n+1], tab.ge[n+1]) * f - 1))
                    @test abs(tab.q[n] - Float64(v)) <= eps(tab.q[n])
                end
                qread, fread = QR._read_errors(tab)
                @test qerr <= max(qread, 6) * u2
                @test ferr <= max(fread, 16) * u2
            end
            N >= QR.WIDE_TABLE_BELOW && @test tab.q[N] == 1.0 && tab.ql[N] == 0.0      # [h − 1] = [1]
            @test all(1 .<= tab.q) && all(0.5 .<= tab.fm .<= 1) && all(0.5 .<= tab.gm .<= 1)
        end
        # above the switch the words agree with the wide build, entry by entry
        a, b = QR.QIntTables(Float64, 300), QR._qint_tables_from(Float64, () -> QR._qints_wide(302, 301))
        @test a.q == b.q && a.qi == b.qi && ldexp.(a.fm, a.fe) == ldexp.(b.fm, b.fe) && ldexp.(a.gm, a.ge) == ldexp.(b.gm, b.ge)
        for (q, N) in ((0.9, 200), (1.1, 200), (exp(-pi / 500), 2000), (1 - 1e-6, 500), (1.5, 1000))
            T = QR.RealQ(q, N)
            setprecision(BigFloat, 800) do
                b = big(q); den = b - inv(b); e = big(0.0)
                for n in 1:N
                    v = (b^n - b^-n) / den
                    e = max(e, abs((big(T.ints[n+1][1]) + T.ints[n+1][2]) / v - 1), abs((big(T.invs[n+1][1]) + T.invs[n+1][2]) * v - 1))
                end
                @test e <= 6 * u2
            end
            @test T.cancel <= 3
        end
        # where the words would leave the exponent range the table stays wide, and stays finite on its inverse side
        T = QR.RealQ(3.0, 2000)
        @test T.invs[2001][1] >= 0 && isinf(T.ints[2001][1]) && T.ints[101][1] > 1e40
        @test q6j(At(3.0), 300, 300, 300, 300, 300, 300) ≈ q6j(At(1 / 3.0), 300, 300, 300, 300, 300, 300) rtol = 1e-12
    end
    @testset "float of a form keeps the 128-bit rounding" begin
        wide(f) = Float64(numeric_value(f; bits = 128))
        for (k, ls) in ((10, ((1, 1, 1, 1, 1, 1), (2, 2, 2, 2, 2, 2), (1, 3 // 2, 5 // 2, 2, 2, 3 // 2))),
                        (5, ((3 // 2, 3 // 2, 1, 3 // 2, 3 // 2, 1),)), (200, ((4, 4, 4, 4, 4, 4), (10, 12, 8, 6, 14, 10))),
                        (1000, ((40, 40, 40, 40, 40, 40),)))
            for l in ls
                f = q6j(Exact(k), l...)
                @test float(f) === Float64(f) === convert(Float64, f) === wide(f)
            end
        end
        z = q6j(Exact(4), 3 // 2, 3 // 2, 1, 3 // 2, 3 // 2, 1)
        @test iszero(z) == iszero(float(z))
        @test float(q3j(Exact(5), 1, 1, 1, 1, -1, 0)) isa ComplexF64
        @test_throws ArgumentError float(q6j(Symbolic(), 1, 1, 1, 1, 1, 1))
    end
    @testset "a column run allocates nothing once warm, on any thread" begin
        l = (200, 210, 190, 205, 195, 215)       # labels whose sum cancels, so the value comes from the column
        call(l) = (q6j(l...), q6j(Level(2001), l...))
        allocations(l) = (call(l); call(l); (@allocated(q6j(l...)), @allocated(q6j(Level(2001), l...))))   # inside a function
        @test allocations(l) == (0, 0)
        want = call(l)
        got = fill((0.0, 0.0), 200)
        Threads.@threads for i in eachindex(got)
            got[i] = call(l)
        end
        @test all(==(want), got)
        @test length(QRecoupling._ratio_slots()) >= Threads.nthreads()
    end
    @testset "x_form at a level" begin
        v = x_form(q6j(Exact(10), 1, 1, 1, 1, 1, 1))
        @test v isa LevelXValue && string(v) == "LevelXValue(k = 10, (2x² − 7)/3)"
        @test occursin("(2x² − 7)/3", sprint(show, MIME"text/plain"(), v)) && occursin("2cos(π/12)", sprint(show, MIME"text/plain"(), v))
        @test string(x_form(q6j(Exact(3), 1, 1, 1, 1, 1, 1))) == "LevelXValue(k = 3, x − 2)"
        @test v == x_form(q6j(Exact(10), 1, 1, 1, 1, 1, 1)) && hash(v) == hash(x_form(q6j(Exact(10), 1, 1, 1, 1, 1, 1)))
        @test iszero(x_form(q6j(Exact(10), 1, 3//2, 3//2, 2, 3//2, 3//2)))
        @test string(x_form(q6j(Exact(10), 1, 3//2, 3//2, 2, 3//2, 3//2))) == "LevelXValue(k = 10, 0)"
        @test_throws ArgumentError x_form(q3j(Exact(5), 1, 1, 1, 1, -1, 0))
        # the printed polynomials are the value: read them back at x = 2cos(π/h), with and without a root
        rooted = 0
        for (f, l) in ((q6j, (1, 1, 1, 1, 1, 1)), (q6j, (2, 2, 2, 2, 2, 2)), (q6j, (3//2, 3//2, 1, 1, 1, 1)), (fsymbol, (1, 1, 1, 1, 1, 1)),
                       (fsymbol, (1, 1, 1, 1, 1, 2)), (gsymbol, (1, 1, 1, 1, 1, 1)), (fsymbol, (3//2, 3//2, 1, 3//2, 3//2, 2))), k in 4:3:40
            form = f(Exact(k), l...)
            iszero(form) && continue
            x = x_form(form)
            setprecision(BigFloat, 256) do
                x0 = 2cospi(inv(big(k + 2)))
                P, R = QRecoupling.xpolynomial(x), QRecoupling.radicand(x)
                rooted += length(R) > 1
                @test isapprox(evalpoly(x0, P) * sqrt(evalpoly(x0, R)), numeric_value(form; bits = 256); rtol = big(2.0)^-200)
            end
            @test occursin("k = $k", string(x)) && float(x) ≈ float(form)
        end
        @test rooted > 0
    end
    @testset "the message when the Nemo backend is missing" begin
        @test occursin("Install Nemo", QRecoupling._algebra_hint(:x_form, false, true))
        @test occursin("`include`", QRecoupling._algebra_hint(:x_form, true, false))
        @test occursin("Restart Julia", QRecoupling._algebra_hint(:x_form, true, true))
    end
    include("qcg_tests.jl")
    include("regge_tests.jl")
    include("release_tests.jl")
end
