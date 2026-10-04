# Quantum Clebsch–Gordan coefficients and the quantum 3j symbol: a closed form, the classical limit,
# orthogonality of whole matrices, and the entry points. The full suite is kept locally.

@testset "qcg and q3j" begin
    @testset "spin-half singlet in closed form" begin
        # |0 0⟩ = (√q |+−⟩ − |−+⟩/√q) / √[2]
        for kw in ((; q = 0.8), (; q = 0.7 + 0.4im), (; k = 8))
            q = haskey(kw, :k) ? cispi(1 / (kw.k + 2)) : complex(kw.q)
            d = sqrt(q + inv(q))
            @test qcg(1//2, 1//2, 1//2, -1//2, 0; kw...) ≈ sqrt(q) / d atol = 1e-14
            @test qcg(1//2, -1//2, 1//2, 1//2, 0; kw...) ≈ -inv(sqrt(q)) / d atol = 1e-14
            @test q3j(1//2, 1//2, 0, 1//2, -1//2, 0; kw...) ≈ sqrt(q) / d atol = 1e-14
        end
    end

    @testset "classical limit" begin
        @test q3j(2, 1, 2, 1, 0, -1) === QR.q3j_factorial(2, 1, 2, 1, 0, -1)
        @test qcg(1, 1, 1, -1, 1) ≈ 1 / sqrt(2)
        @test qcg(Exact(), 1//2, 1//2, 1//2, -1//2, 0)^2 == 1//2
    end

    @testset "qcg_matrix and qcg_row" begin
        for kw in ((;), (; q = 0.8), (; k = 12), (; q = 0.7 + 0.4im))
            C, prod, coup = qcg_matrix(2, 3//2; kw...)
            n = size(C, 2)
            @test maximum(abs.(transpose(C) * C - [i == j for i in 1:n, j in 1:n])) < 1e-13
            @test all(isapprox(C[a, b], m1 + m2 == m ? qcg(2, m1, 3//2, m2, j, m; kw...) : 0; atol = 1e-14)
                      for (a, (m1, m2)) in enumerate(prod), (b, (j, m)) in enumerate(coup))
        end
        c, js = qcg_row(2, 1, 3//2, -1//2; q = 0.8)
        @test c ≈ [qcg(2, 1, 3//2, -1//2, j; q = 0.8) for j in js] atol = 1e-14
    end

    @testset "entry points and refusals" begin
        l = (3//2, 1//2, 1, -1, 3//2, -1//2)
        @test qcg((3//2, 1//2), (1, -1), (3//2, -1//2)) === qcg(l...)
        @test qcg(Level(10), (3//2, 1//2), (1, -1), (3//2, -1//2)) === qcg(l...; k = 10)
        @test qcg(At(0.8), l...) === qcg(l...; q = 0.8)
        @test qcg([l, l]; q = 0.8) == fill(qcg(l...; q = 0.8), 2)
        @test qcg(1, 1, 1, 0, 1, 0) == 0                    # m ≠ m₁ + m₂
        @test_throws ArgumentError q3j(Exact(5), 1, 1, 1, 1, -1, 0)
        @test_throws ArgumentError qcg(Symbolic(), 1, 1, 1, -1, 1)
    end
end
