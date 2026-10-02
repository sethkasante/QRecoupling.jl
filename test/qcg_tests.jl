# ---------------------------------------------------------------------------------
#  Quantum Clebsch–Gordan coefficients and the quantum 3j symbol (src/qcg.jl).
#
#  The reference is the definition, not another code path: the coupled vectors must intertwine the
#  U_q(sl₂) action (K|m⟩ = q^m|m⟩, E|m⟩ = √([j−m][j+m+1])|m+1⟩, Δ(E) = E⊗K + K⁻¹⊗E). Orthogonality,
#  the symmetry relations and the q = 1 limit are checked on top of that.
# ---------------------------------------------------------------------------------
@testset "qcg and q3j" begin
    _qn(n, q) = (q^n - q^(-n)) / (q - 1 / q)
    # the few matrix helpers needed, so the tests take no extra dependency
    function _kron(A, B)
        p, q = size(B)
        return [A[(r - 1) ÷ p + 1, (c - 1) ÷ q + 1] * B[(r - 1) % p + 1, (c - 1) % q + 1]
                for r in 1:size(A, 1) * p, c in 1:size(A, 2) * q]
    end
    _diag(v) = [i == j ? v[i] : zero(eltype(v)) for i in eachindex(v), j in eachindex(v)]
    _eye(n) = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]
    pairs = ((1//2, 1//2), (1, 1), (2, 1), (3//2, 1//2), (5//2, 3//2), (3, 2))

    @testset "intertwiner: the defining property" begin
        function rep(j, q)
            ms = collect(j:-1:-j)
            K = _diag([q^m for m in ms]); E = zeros(typeof(q), length(ms), length(ms))
            for (i, m) in enumerate(ms)
                i > 1 && (E[i-1, i] = sqrt(_qn(j - m, q) * _qn(j + m + 1, q)))
            end
            return ms, K, E
        end
        for q in (0.8, 1.3), (j1, j2) in pairs
            m1s, K1, E1 = rep(j1, q); m2s, K2, E2 = rep(j2, q)
            ΔK = _kron(K1, K2); ΔE = _kron(E1, K2) + _kron(_diag([1 / K1[i, i] for i in axes(K1, 1)]), E2)
            for j in abs(j1 - j2):(j1 + j2)
                ms, Kj, Ej = rep(j, q)
                V = zeros(length(m1s) * length(m2s), length(ms))
                for (c, m) in enumerate(ms), (a, x) in enumerate(m1s), (b, y) in enumerate(m2s)
                    x + y == m && (V[(a - 1) * length(m2s) + b, c] = qcg(j1, x, j2, y, j, m; q = q))
                end
                @test maximum(abs, ΔK * V - V * Kj) < 1e-13
                @test maximum(abs, ΔE * V - V * Ej) < 1e-12
            end
        end
    end

    @testset "q = 1 is the classical symbol, bit for bit" begin
        for (j1, j2) in pairs, j3 in abs(j1 - j2):(j1 + j2), m1 in -j1:j1, m2 in -j2:j2
            m3 = -m1 - m2
            abs(m3) <= j3 || continue
            @test q3j(j1, j2, j3, m1, m2, m3) === q3j_factorial(j1, j2, j3, m1, m2, m3)
            @test q3j(Exact(), j1, j2, j3, m1, m2, m3) == q3j_factorial(Exact(), j1, j2, j3, m1, m2, m3)
            cg = (-1)^Int(j1 - j2 - m3) * sqrt(2j3 + 1) * q3j_factorial(j1, j2, j3, m1, m2, m3)
            @test qcg(j1, m1, j2, m2, j3) ≈ cg atol = 1e-15
        end
        @test qcg(Exact(), 1//2, 1//2, 1//2, -1//2, 0)^2 == 1//2
        @test q3j(1, 1, 1, 1, -1, 0; q = 1) === q3j(1, 1, 1, 1, -1, 0)
    end

    @testset "orthogonality, C Cᵀ = 1" begin
        function defect(f)
            worst = 0.0
            for (j1, j2) in pairs, m in -(j1 + j2):(j1 + j2)
                js = [j for j in abs(j1 - j2):(j1 + j2) if abs(m) <= j]
                m1s = [m1 for m1 in -j1:j1 if abs(m - m1) <= j2]
                C = [f(j1, m1, j2, m - m1, j) for j in js, m1 in m1s]
                worst = max(worst, maximum(abs, C * transpose(C) - _eye(length(js))))
            end
            return worst
        end
        @test defect((a...) -> qcg(a...; q = 0.8)) < 1e-13
        @test defect((a...) -> qcg(a...; q = 1.3)) < 1e-13
        @test defect((a...) -> qcg(a...; q = 0.7 + 0.4im)) < 1e-13
        # at a level, only when the whole block is inside the fusion rule (2(j₁ + j₂) ≤ k for these pairs)
        @test defect((a...) -> qcg(Level(11), a...)) < 1e-13
        # the substituted formula is not orthogonal: this is why it is no longer `q3j`
        old = (j1, m1, j2, m2, j) -> (-1)^Int(j1 - j2 + m1 + m2) * sqrt(qint(Int(2j + 1); q = 0.8)) *
                                     q3j_factorial(j1, j2, j, m1, m2, -m1 - m2; q = 0.8)
        @test defect(old) > 0.5
    end

    @testset "symmetries" begin
        for q in (0.8, 1.3), (j1, j2) in pairs, j3 in abs(j1 - j2):(j1 + j2), m1 in -j1:j1, m2 in -j2:j2
            m3 = -m1 - m2
            abs(m3) <= j3 || continue
            x = q3j(j1, j2, j3, m1, m2, m3; q = q)
            σ = (-1)^Int(j1 + j2 + j3)
            @test x ≈ σ * q3j(j2, j1, j3, m2, m1, m3; q = 1 / q) atol = 1e-14
            @test x ≈ σ * q3j(j1, j2, j3, -m1, -m2, -m3; q = 1 / q) atol = 1e-14
            @test q3j(j2, j3, j1, m2, m3, m1; q = q) ≈ q^m2 * x atol = 1e-14
        end
    end

    @testset "levels" begin
        # the level kernel takes every number at the root itself; for small labels the analytic path at the
        # rounded root agrees with it
        for k in (3, 7, 12), (j1, j2) in pairs, j in abs(j1 - j2):(j1 + j2), m1 in -j1:j1, m2 in -j2:j2
            abs(m1 + m2) <= j || continue
            v = qcg(Level(k), j1, m1, j2, m2, j)
            @test v isa ComplexF64
            if j1 + j2 + j <= k
                @test v ≈ qcg(j1, m1, j2, m2, j; q = cispi(1 / (k + 2))) atol = 1e-12
            else
                @test v == 0
            end
        end
        b = qcg(Level(10; T = BigFloat), 3, 1, 2, -1, 4)
        @test b isa Complex{BigFloat}
        @test abs(ComplexF64(b) - qcg(Level(10), 3, 1, 2, -1, 4)) < 1e-15
        @test qcg(1, 1, 1, -1, 1; k = 1:3) == [qcg(Level(k), 1, 1, 1, -1, 1) for k in 1:3]
    end

    @testset "conventions and refusals" begin
        @test qcg(1, 1, 1, 0, 1, 0) == 0                   # m ≠ m₁ + m₂
        @test qcg(1, 2, 1, 0, 1) == 0                      # |m₁| > j₁
        @test q3j(1, 1, 3, 0, 0, 0; q = 0.8) == 0          # triangle
        @test qcg(1//2, 1//2, 1, 0, 1//2; q = -0.8) isa Real
        @test qcg(1//2, 1//2, 1//2, -1//2, 0; q = 0.8) isa Real
        @test qcg(At(0.8), 1, 1, 1, -1, 1) === qcg(1, 1, 1, -1, 1; q = 0.8)
        @test q3j(Level(5), 1, 1, 1, 1, -1, 0) === q3j(1, 1, 1, 1, -1, 0; k = 5)
        L = [(1, 1, 1, 1, -1, 0), (2, 2, 2, 0, 0, 0), (2, 1, 1, 1, 0, -1)]
        @test q3j(L; q = 0.8) == [q3j(l...; q = 0.8) for l in L]
        @test qcg([(1, 1, 1, -1, 1, 0), (1, 0, 1, 0, 2, 0)]; k = 5) ==
              [qcg(1, 1, 1, -1, 1, 0; k = 5), qcg(1, 0, 1, 0, 2, 0; k = 5)]
        for f in (q3j, qcg), t in (Exact(5), Symbolic())
            args = f === q3j ? (1, 1, 1, 1, -1, 0) : (1, 1, 1, -1, 1)
            @test_throws ArgumentError f(t, args...)
            @test_throws ArgumentError f(args...; k = 5, exact = true)
        end
        @test_throws ArgumentError qcg(1, 1, 1, -1, 1; q = 0.8, k = 5)
    end
end
