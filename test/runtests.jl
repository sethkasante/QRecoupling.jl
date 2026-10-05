using Test, QRecoupling

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
        @test v isa ExactX
        @test Float64(v) ≈ formula(2cospi(1/12)) rtol=1e-13
        @test q6j(Level(10), ONES...) ≈ Float64(v) rtol=1e-13
        @test string(radical(v)) == "(2√3 − 3)/3"
        @test iszero(v-v) && isone(v*inv(v))

        s = q6j(Symbolic(), ONES...)
        @test s isa SymbolicValue && x_form(s) isa XValue
        @test qeval(Exact(), s) == 1//6
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
        @test abs(rmatrix(1,1,1; k=5)) ≈ 1
        @test abs(twist(1//2; k=5)) ≈ 1
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

    include("qcg_tests.jl")
    include("regge_tests.jl")
    include("release_tests.jl")
end
