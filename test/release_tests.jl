@testset "v0.4 release regressions" begin
    @testset "exact equality is invariant under rational scaling" begin
        # Both values are -sqrt(sqrt(3)/3), in different specialized square classes.
        a = q6j(Exact(4), 0,0,0,1//2,1//2,1//2)
        b = q6j(Exact(4), 0,0,0,3//2,3//2,3//2)
        @test a.sqclass != b.sqclass
        @test a == b
        @test iszero(a - b)
        @test !iszero(a + b)
        for scale in (big(1)//big(10)^100, -big(1)//big(10)^100, big(10)^100)
            @test !iszero(scale * (a + b))
            @test scale * (a + b) != 0
            @test iszero(scale * (a - b))
        end
        tiny = big(1)//big(10)^100
        @test !iszero((a - b) + tiny)
        @test (a - b) + tiny == tiny
        @test (a - b) + 1 == 1
    end

    @testset "parallel zero flags own their storage" begin
        for n in (0, 1, 31, 32, 63, 64, 65, 127, 128, 129)
            allzero = fill((1,1,1,1,1,4), n)
            mixed = [isodd(i) ? (1,1,1,1,1,4) : (1,1,1,1,1,1) for i in 1:n]
            for labels in (allzero, mixed)
                expected = [iszero_at(6, js...) for js in labels]
                for nw in (1, 2, 3, 4), _ in 1:25
                    actual = iszero_at(6, labels; threads=nw)
                    @test actual isa BitVector
                    @test actual == expected
                end
            end
        end
    end

    @testset "modular tables separate numeric types and precisions" begin
        # Exercise both cache insertion orders at fresh levels.
        smatrix(11)
        setprecision(BigFloat, 64) do
            for k in (11, 13)
                S, js = smatrix(k; T=BigFloat)
                @test eltype(S) === BigFloat
                @test all(x -> precision(x) == 64, S)
                h = k + 2
                expected = sqrt(BigFloat(2)/h) * sin(big(pi)/h)
                @test S[1,1] ≈ expected rtol=eps(BigFloat)*4
            end
        end
        @test eltype(first(smatrix(13))) === Float64
        setprecision(BigFloat, 128) do
            S, = smatrix(11; T=BigFloat)
            @test all(x -> precision(x) == 128, S)
            @test S[1,1] ≈ sqrt(BigFloat(2)/13)*sin(big(pi)/13) rtol=eps(BigFloat)*4
        end
    end

    @testset "F matrices retain requested precision" begin
        setprecision(BigFloat, 256) do
            for kw in ((;), (;k=6), (;q=BigFloat("0.8")),
                       (;q=Complex{BigFloat}(BigFloat("0.8"),BigFloat("0.3"))))
                E = haskey(kw, :q) && kw.q isa Complex ? Complex{BigFloat} : BigFloat
                F, es, fs = fmatrix(1,1,1,1; T=E, kw...)
                reference = [fsymbol(1,1,e,1,1,f; T=BigFloat, kw...) for e in es, f in fs]
                @test eltype(F) === E
                @test maximum(abs, F-reference) < big"1e-65"
                identity = [i == j ? one(E) : zero(E) for i in eachindex(fs), j in eachindex(fs)]
                @test maximum(abs, transpose(F)*F-identity) < big"1e-65"
            end
            F, = fmatrix(1,1,1,1; k=6, T=BigFloat)
            @test any(x -> x != BigFloat(Float64(x)), F)
            @test eltype(first(fmatrix(1,1,1,1; q=BigFloat("0.8")))) === BigFloat
            @test eltype(first(fmatrix(1,1,1,1; q=Complex{BigFloat}(0.8,0.3)))) === Complex{BigFloat}
            @test eltype(first(fmatrix(1,1,1,1; k=6, T=Complex{BigFloat}))) === Complex{BigFloat}
        end
    end
end
