using Test
using QRecoupling

const QR = QRecoupling

# ---------------------------------------------------------------------------------
#  Everything here is a cross-check between two paths that were written separately,
#  or a value that can be read off by hand, so a failure says which path moved.
# ---------------------------------------------------------------------------------

const LABELS = ((1, 1, 1, 1, 1, 1),
                (2, 2, 2, 2, 2, 2),
                (1//2, 1//2, 1, 1//2, 1//2, 1),
                (3//2, 1, 3//2, 1, 3//2, 1),
                (2, 3//2, 5//2, 2, 5//2, 3//2))

rel(a, b) = abs(a - b) / max(1.0, abs(b))

@testset "QRecoupling.jl" begin

    @testset "classical" begin
        @test q6j(1, 1, 1, 1, 1, 1) ≈ 1 / 6
        @test q6j(1, 1, 1, 1, 1, 1; exact = true) == 1 // 6
        @test q6j(0, 0, 0, 0, 0, 0) == 1.0
        @test q6j(1, 1, 1, 1, 1, 4) == 0.0            # inadmissible is zero, not an error
        # the classical value is the q → 1 limit of the analytic path
        for js in LABELS
            @test rel(q6j(js...; q = 1), q6j(js...)) < 1e-12
            @test rel(q6j(js...; q = 1 + 1e-9), q6j(js...)) < 1e-6
        end
        # column exchange and the (12)(45) swap leave a 6j symbol alone
        a, b, c, d, e, f = 2, 3//2, 5//2, 2, 5//2, 3//2
        @test q6j(a, b, c, d, e, f) ≈ q6j(b, a, c, e, d, f)
        @test q6j(a, b, c, d, e, f) ≈ q6j(d, e, c, a, b, f)
    end

    @testset "exact classical target" begin
        @test Exact() === Classical(exact=true)
        @test q6j(Exact(), 1, 1, 1, 1, 1, 1) == 1//6
        @test q6j(Exact(), 1, 1, 1, 1, 1, 4) == 0
        r = q3j(Exact(), 1//2, 1//2, 0, 1//2, -1//2, 0)
        @test r^2 == 1//2
        for (f, args) in ((q6j, LABELS[4]), (q3j, (1,1,1,1,-1,0)),
                          (fsymbol, LABELS[1]), (gsymbol, LABELS[1]),
                          (QR.tetrahedron, LABELS[1]), (QR.theta_value, (1,1,0)),
                          (qdim, (1//2,)), (rmatrix, (1,1,1)), (twist, (1//2,)),
                          (qint, (5,)), (qfact, (5,)), (qbinomial, (5,2)))
            v = f(Exact(), args...)
            legacy = f(args...; q=1, exact=true)
            @test v == legacy
            @test typeof(v) === typeof(legacy)
        end
        s = q6j(Symbolic(), LABELS[1]...)
        @test qeval(Exact(), s) == 1//6
        @test qeval(Exact(), s.rule) == 1//6
        @test q6j(Exact(), collect(LABELS)) == [q6j(Exact(), js...) for js in LABELS]
        for kw in ((;q=1), (;k=10), (;exact=true), (;exact=false))
            @test_throws ArgumentError q6j(Exact(), LABELS[1]...; kw...)
        end
    end

    @testset "level k" begin
        @test q6j(1, 1, 1, 1, 1, 1; k = 10) ≈ 0.1547005383792515
        for js in LABELS, k in (5, 10, 40)
            QR._qδtet(QR.doubled(js...)..., k) || continue
            v = q6j(js...; k = k)
            @test isfinite(v)
            # three independent routes to the same number: the level kernel, the analytic path at
            # q = exp(iπ/h), and the exact value in ℚ[x]/Ψ_h
            @test rel(v, real(q6j(js...; q = cispi(1 / (k + 2))))) < 1e-10
            @test rel(v, Float64(q6j(Exact(k), js...))) < 1e-12
        end
        # large level and large spin stay finite and agree with the classical limit
        @test rel(q6j(1, 1, 1, 1, 1, 1; k = 100_000), q6j(1, 1, 1, 1, 1, 1)) < 1e-8
        @test isfinite(q6j(30, 30, 30, 30, 30, 30; k = 200))
    end

    @testset "real and complex q" begin
        for js in LABELS
            for q in (0.6, 0.87, 1.4)
                v = q6j(js...; q = q)
                @test isfinite(v)
                @test rel(v, q6j(js...; q = 1 / q)) < 1e-10      # q ↔ 1/q invariance
            end
            for q in (0.87 + 0.28im, 0.3 + 0.9im)
                v = q6j(js...; q = q)
                @test rel(v, conj(q6j(js...; q = conj(q)))) < 1e-10
            end
        end
        @test q6j(1, 1, 1, 1, 1, 1; q = 0.87, T = BigFloat) isa BigFloat
    end

    @testset "exact values in x, and radicals" begin
        v = q6j(Exact(10), 1, 1, 1, 1, 1, 1)
        @test v isa ExactX && QR.level(v) == 10
        @test rel(Float64(v), q6j(1, 1, 1, 1, 1, 1; k = 10)) < 1e-14
        @test string(radical(v)) == "(2√3 − 3)/3"
        @test string(radical(q6j(Exact(3), 1, 1, 1, 1, 1, 1))) == "(√5 − 3)/2"
        # printing is the stored form and nothing computed
        txt = sprint(show, MIME"text/plain"(), v)
        @test occursin("(2x² − 7)/3", txt) && !occursin("≈", txt) && !occursin("√", txt)
        @test occursin("≈", sprint(show, MIME"text/plain"(), v; context = :approximate => true))
        P, R = v.x_value
        @test P === v.p && isone(R)
        # a level whose degree has an odd prime factor has no radical form, and says so
        n = radical(q6j(Exact(5), 1, 1, 1, 1, 1, 1))
        @test n isa NoRadical && n.kind === :none
        @test occursin("no radical form", repr(MIME"text/plain"(), n))
        @test has_radical_form(3) && !has_radical_form(5)
        # arithmetic, including scalars, and the one-class collapse
        w = q6j(Exact(10), 3//2, 1, 3//2, 1, 3//2, 1)
        @test v * w isa ExactX && v + 1 isa ExactX
        @test rel(Float64(v * w), Float64(v) * Float64(w)) < 1e-13
        @test rel(Float64(v + 1), Float64(v) + 1) < 1e-13
        @test iszero(v - v) && isone(v * inv(v))
        @test_throws ArgumentError v + 0.5
        @test_throws ArgumentError q6j(Exact(10), 1, 1, 1, 1, 1, 1) * q6j(Exact(8), 1, 1, 1, 1, 1, 1)
    end

    @testset "symbolic rules" begin
        v = q6j(Symbolic(), 1, 1, 1, 1, 1, 1)
        @test v isa QR.SymbolicValue
        @test occursin("Σ", sprint(show, MIME"text/plain"(), v))
        @test x_form(v) == QR.generic_sixj(QR.doubled(1, 1, 1, 1, 1, 1)...)
        @test occursin("x² − 3", sprint(show, MIME"text/plain"(), x_form(v)))
        @test splits_completely(phi_form(v)) isa Bool
        @test v.dcr isa QR.DCR
        # the rule evaluates to the same numbers the direct calls give
        for k in (5, 10), js in (LABELS[1], LABELS[3])
            QR._qδtet(QR.doubled(js...)..., k) || continue
            @test rel(qeval(q6j(Symbolic(), js...).dcr; k = k), q6j(js...; k = k)) < 1e-12
        end
        @test rel(qeval(v.dcr; q = 0.7), q6j(1, 1, 1, 1, 1, 1; q = 0.7)) < 1e-12
    end

    @testset "q-numbers take the symbols' interface" begin
        @test qint(5) === 5.0 && qint(0) === 1.0        # [0] = 1, the package's convention
        @test abs(qint(25; k = 10) - sinpi(25 / 12) / sinpi(1 / 12)) < 1e-14   # periodic past h
        @test qfact(4) == 24.0 && qbinomial(5, 2) == 10.0
        @test qint(5; exact = true) == 5
        @test rel(qint(5; k = 10), sin(5 * pi / 12) / sin(pi / 12)) < 1e-12
        @test rel(qint(5; q = 0.8), (0.8^5 - 0.8^-5) / (0.8 - 0.8^-1)) < 1e-12
        @test qint(Symbolic(), 5) isa QR.SymbolicValue
        @test qint(Exact(10), 5) isa ExactX
        @test rel(Float64(qint(Exact(10), 5)), qint(5; k = 10)) < 1e-13
        @test rel(qdim(1; k = 5), real(qdim(1; q = cispi(1 / 7)))) < 1e-10
    end

    @testset "the other tensors" begin
        @test q3j(1, 1, 1, 0, 0, 0) isa Float64
        @test fsymbol(1, 1, 1, 1, 1, 1; k = 5) ≈ 0.19806226419516196
        @test gsymbol(1, 1, 1, 1, 1, 1; k = 5) ≈ 1.0 atol = 1e-12
        @test abs(rmatrix(1, 1, 1; k = 5)) ≈ 1.0 atol = 1e-12
        @test qdim(1//2; k = 5) > 0
        @test qdim(1//2; k = 5, exact = true) isa ExactX
        # a q-series through the generic constructor
        # the summand builds monomials: in v0.4 `qfact(z)` is a number, `qfact_mono(z)` the rule
        s = qseries(1:6) do z
            (-1)^z * qfact_mono(z)
        end
        @test isfinite(qeval(s; k = 10))
        @test rel(qeval(s; q = 1), sum((-1)^z * factorial(z) for z in 1:6)) < 1e-10
    end

    @testset "refusals are errors, not wrong answers" begin
        @test_throws DomainError q6j(1, 1, 1, 1, 1, 1; k = -1)
        @test_throws ArgumentError q6j(1//3, 1, 1, 1, 1, 1)
        @test_throws DomainError has_radical_form(-1)
        @test_throws DivideError inv(zero(ExactX, 10))
        # an identically vanishing symbol at generic q is a zero, not an error
        @test q3j(5, 5, 5, 1, -2, 1; q = 0.8) == 0.0
        @test q6j(1, 1, 1, 1, 1, 1; k = 1) == 0.0       # inadmissible at this level
    end

    # # The full suite lives outside the package: twenty files, ~167,000 assertions,
    # # three minutes. Kept out of CI on purpose; run it before a release.
    # if get(ENV, "QRECOUPLING_FULL_TESTS", "") != ""
    #     full = joinpath(@__DIR__, "..", "dev", "tests", "runtests.jl")
    #     isfile(full) ? include(full) :
    #         @warn "QRECOUPLING_FULL_TESTS is set but $full is not present"
    # end
end
