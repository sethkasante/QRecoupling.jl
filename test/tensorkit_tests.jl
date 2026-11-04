using Test, LinearAlgebra, QRecoupling, TensorKit, TensorKitSectors
using TensorKitSectors: pentagon_equation, hexagon_equation, triangle_equation, frobenius_schur_phase, findindex
const QR = QRecoupling

include(joinpath(dirname(dirname(pathof(TensorKitSectors))), "test", "testsuite.jl"))

"Spins as sectors, up to the fusion rules of the type."
spins(I, top) = [I(j) for j in 0:1//2:top]

"The largest entry of F^{abc}_d F^{abc†}_d − 1 over all channels."
function unitarity_defect(a, b, c)
    worst = 0.0
    for d in ⊗(a, b, c)
        es = [e for e in ⊗(a, b) if d in ⊗(e, c)]
        fs = [f for f in ⊗(b, c) if d in ⊗(a, f)]
        F = [Fsymbol(a, b, c, d, e, f) for e in es, f in fs]
        worst = max(worst, maximum(abs, F' * F - I))
    end
    return worst
end

"The spectrum of −Σ P0(i,i+1) on L spins 1/2, with P0 the projector on the trivial channel."
function spin_half_chain(I, L)
    P = Vect[I](I(1//2) => 1); T = sectorscalartype(I)
    h = zeros(T, P ⊗ P ← P ⊗ P)
    for (c, b) in blocks(h)
        b .= isunit(c) ? 1.0 : 0.0
    end
    H = zeros(T, P^L ← P^L)
    for i in 1:L-1
        term = h
        i > 1 && (term = id(P^(i - 1)) ⊗ term)
        i < L - 1 && (term = term ⊗ id(P^(L - i - 1)))
        H = H - term
    end
    return sort!(vcat([eigvals(Hermitian(Matrix(b))) for (_, b) in blocks(H)]...))
end

"The same chain written on (ℂ²)^L, where the Temperley–Lieb generator e = d·P0 is a 4 × 4 matrix (Pasquier–Saleur)."
function xxz_chain(q, L)
    e = ComplexF64[0 0 0 0; 0 q -1 0; 0 -1 inv(q) 0; 0 0 0 0] / (q + inv(q))
    one(n) = Matrix{ComplexF64}(I, 2^n, 2^n)
    return eigvals(-sum(kron(one(i - 1), e, one(L - i - 1)) for i in 1:L-1))
end

@testset "TensorKit sectors" begin
    @test Base.get_extension(QR, :QRecouplingTensorKitExt) !== nothing
    @test Base.get_extension(QR, :QRecouplingNemoExt) === nothing

    # the interface as TensorKit itself tests it, including pentagon, hexagon, triangle, A-, B-symbols and show-and-parse,
    # so that a change in TensorKitSectors shows up here: an abelian level, Fibonacci, a level with half-integer spins and
    # one past the size where spaces switch to dictionaries, and real q on both sides of 1
    for I in (SU2k(1), SU2k(3), SU2k(10), SU2k(40), SU2q(0.8), SU2q(1.3))
        SectorTestSuite.test_sector(I)
    end

    @testset "construction" begin
        I = SU2k(10)
        @test I === SU2k(10) && I(1//2) == I(0.5) && I(2) == I(2.0) && I(3//2).J == 3
        @test convert(I, 1) == I(1) && unit(I) == I(0) && dual(I(5//2)) == I(5//2) && isunit(I(0)) && !isunit(I(1))
        @test isbitstype(I) && isbitstype(SU2q(0.8))
        @test_throws DomainError I(11//2)
        @test_throws DomainError I(-1)
        @test_throws DomainError I(1//3)
        @test_throws DomainError SU2k(-1)
        @test_throws DomainError SU2q(0.0)
        @test_throws DomainError SU2q(-0.8)
        @test_throws DomainError SU2q(Inf)
        @test string(I) == "SU2k(10)" && string(I(3//2)) == "SU2k(10)(3//2)" && string(SU2q(0.8)) == "SU2q(0.8)"
        @test sprint(show, Vect[I](1//2 => 2)) == "Vect[SU2k(10)](1//2 => 2)"
        @test length(values(I)) == 11 && collect(values(I)) == spins(I, 5) && values(I)[3] == I(1)
        @test findindex(values(I), I(1)) == 3
        @test_throws BoundsError values(I)[12]
        @test Base.IteratorSize(values(SU2q(0.8))) isa Base.IsInfinite
        @test sectorscalartype(I) === ComplexF64 && fusionscalartype(I) === Float64 && sectorscalartype(SU2q(0.8)) === Float64
    end

    @testset "scalars are inferred" begin
        a = SU2k(10)(1)
        @test (@inferred Fsymbol(a, a, a, a, a, a)) isa Float64
        @test (@inferred Rsymbol(a, a, a)) isa ComplexF64
        @test (@inferred dim(a)) isa Float64
        b = SU2q(0.8)(1)
        @test (@inferred Fsymbol(b, b, b, b, b, b)) isa Float64
        @test (@inferred dim(b)) isa Float64
    end

    @testset "fusion rules are truncated by the level" begin
        I = SU2k(4)
        @test collect(I(1) ⊗ I(3//2)) == [I(1//2), I(3//2)]
        @test collect(I(2) ⊗ I(2)) == [I(0)]
        @test length(I(1) ⊗ I(1)) == 3
        @test Nsymbol(I(1), I(1), I(2)) && !Nsymbol(I(3//2), I(3//2), I(2)) && !Nsymbol(I(1), I(1), I(1//2))
        J = SU2q(0.8)
        @test collect(J(1) ⊗ J(3//2)) == [J(1//2), J(3//2), J(5//2)] && length(J(1) ⊗ J(1)) == 3
        # a channel fusion forbids has no symbol
        @test Fsymbol(I(1), I(1), I(1), I(1), I(1//2), I(1)) == 0.0 && Rsymbol(I(1), I(1), I(1//2)) == 0
        @test Fsymbol(J(1), J(1), J(1), J(1), J(1//2), J(1)) == 0.0
        # dimensions multiply as the fusion rules say, at every level
        for k in (1, 3, 12), a in spins(SU2k(k), min(k // 2, 3)), b in spins(SU2k(k), min(k // 2, 3))
            @test dim(a) * dim(b) ≈ sum(dim, a ⊗ b)
        end
    end

    @testset "every associativity and braiding equation, small levels" begin
        for k in (2, 4)
            I = SU2k(k); all_ = collect(values(I))
            @test all(triangle_equation(a, b; atol = 1e-13) for a in all_, b in all_)
            @test all(pentagon_equation(a, b, c, d; atol = 1e-12) for a in all_, b in all_, c in all_, d in all_)
            @test all(hexagon_equation(a, b, c; atol = 1e-12) for a in all_, b in all_, c in all_)
            @test all(unitarity_defect(a, b, c) < 1e-13 for a in all_, b in all_, c in all_)
            @test all(abs(abs(Rsymbol(a, b, c)) - 1) < 1e-14 for a in all_, b in all_ for c in a ⊗ b)
        end
    end

    @testset "the equations hold at large levels and spins" begin
        # one large spin among small ones, in each position, since a large spin fused with a large spin has thousands of channels
        for k in (200, 5000)
            K = SU2k(k); small = [K(1//2), K(1)]
            for m in (K((k ÷ 2) // 2), K(k // 2)), a in small, b in small, c in small
                @test pentagon_equation(m, a, b, c; atol = 1e-9) && pentagon_equation(a, m, b, c; atol = 1e-9) &&
                      pentagon_equation(a, b, m, c; atol = 1e-9)
                @test hexagon_equation(m, a, b; atol = 1e-9) && hexagon_equation(a, m, b; atol = 1e-9)
                @test unitarity_defect(m, a, b) < 1e-11 && unitarity_defect(a, m, b) < 1e-11
            end
        end
        for q in (0.3, 2.0)
            I = SU2q(q); s = spins(I, 3//2)
            @test all(pentagon_equation(a, b, c, d; atol = 1e-8) for a in s, b in s, c in s, d in s[1:3])
            @test all(unitarity_defect(a, b, c) < 1e-10 for a in s, b in s, c in s)
        end
    end

    @testset "against independent data" begin
        # q = 1 is SU(2), and large level approaches it
        for a in spins(SU2q(1.0), 1), b in spins(SU2q(1.0), 1), c in spins(SU2q(1.0), 3//2), d in ⊗(a, b, c)
            U = SU2Irrep
            @test dim(a) == dim(U(a.J // 2))
            for e in a ⊗ b, f in b ⊗ c
                @test Fsymbol(a, b, c, d, e, f) ≈ Fsymbol(U(a.J // 2), U(b.J // 2), U(c.J // 2), U(d.J // 2), U(e.J // 2), U(f.J // 2)) atol = 1e-13
                K = SU2k(10000); to(x) = K(x.J // 2)
                @test Fsymbol(to(a), to(b), to(c), to(d), to(e), to(f)) ≈ Fsymbol(a, b, c, d, e, f) atol = 1e-5
            end
        end
        # the integer spins of level 3 are the Fibonacci category, gauge-invariantly
        I = SU2k(3); τ = I(1); T = FibonacciAnyon(:τ); u = FibonacciAnyon(:I)
        @test dim(τ) ≈ dim(T)
        @test Fsymbol(τ, τ, τ, τ, τ, τ) ≈ Fsymbol(T, T, T, T, T, T)
        @test Fsymbol(τ, τ, τ, τ, I(0), I(0)) ≈ Fsymbol(T, T, T, T, u, u)
        F(a, b, c, d, e, f) = Fsymbol(a, b, c, d, e, f) * Fsymbol(a, b, c, d, f, e)     # invariant under rescaling e, f alike
        @test F(τ, τ, τ, τ, I(0), τ) ≈ F(T, T, T, T, u, T)
        # the two chiral conjugates of Fibonacci differ in the sign of the exponent
        @test Rsymbol(τ, τ, I(0)) ≈ cispi(4 / 5) || Rsymbol(τ, τ, I(0)) ≈ cispi(-4 / 5)
        @test Rsymbol(τ, τ, τ) ≈ cispi(3 / 5) || Rsymbol(τ, τ, τ) ≈ cispi(-3 / 5)
        @test abs(TensorKitSectors.twist(τ) - TensorKitSectors.twist(T)) < 1e-12 || abs(TensorKitSectors.twist(τ) - conj(TensorKitSectors.twist(T))) < 1e-12
        # level 2 has the fusion rules and F-matrix magnitudes of Ising, with the opposite Frobenius-Schur sign at 1/2
        J = SU2k(2); σ = J(1//2)
        @test dim(σ) ≈ √2 && abs(Fsymbol(σ, σ, σ, σ, J(0), J(0))) ≈ 1 / √2 && frobenius_schur_phase(σ) == -1 && frobenius_schur_phase(J(1)) == 1
        @test collect(σ ⊗ σ) == [J(0), J(1)] && collect(J(1) ⊗ J(1)) == [J(0)] && collect(J(1) ⊗ σ) == [σ]
        # the S-matrix of the package follows from these twists: S_ab = Σ_c d_c θ_c / (θ_a θ_b) / √(Σ_c d_c²)
        for k in (3, 5)
            K = SU2k(k); S, _ = smatrix(k); θ(x) = TensorKitSectors.twist(x)
            for a in values(K), b in values(K)
                s = sum(dim(c) * θ(c) / (θ(a) * θ(b)) for c in a ⊗ b)
                @test abs(s / sum(dim(c)^2 for c in values(K))^0.5) ≈ abs(S[a.J+1, b.J+1]) atol = 1e-12
            end
        end
        # the F-symbols are those of `fsymbol`, in TensorKit's order of arguments
        K = SU2k(12); a, b, c, d, e, f = K(1), K(3//2), K(1), K(1//2), K(1//2), K(3//2)
        @test Fsymbol(a, b, c, d, e, f) == fsymbol(1, 3//2, 1//2, 1, 1//2, 3//2; k = 12)
    end

    @testset "tensors at a level" begin
        # one level only: each sector type compiles TensorKit's tensor methods anew, which is most of the time here
        for k in (3,)
            I = SU2k(k)
            V = Vect[I](0 => 2, 1//2 => 2, 1 => 1)
            W = Vect[I](0 => 1, 1 => 2)
            @test dim(V) ≈ 2 + 2dim(I(1//2)) + dim(I(1)) && sectortype(V) === I
            t = randn(ComplexF64, V ⊗ V ← V ⊗ W)
            # linear algebra and contraction on the coupled blocks
            @test norm(t)^2 ≈ real(tr(t' * t)) && norm(2t) ≈ 2norm(t) && t + t ≈ 2t
            @test t' * t ≈ (t' * t)' && isapprox(dot(t, t), norm(t)^2)
            # F-moves are unitary: the norm survives any repartition, and it can be undone
            r = repartition(t, 1, 3)
            @test norm(r) ≈ norm(t) && repartition(r, 2, 2) ≈ t
            @test norm(repartition(t, 3, 1)) ≈ norm(t)
            # planar and braided permutations, undone
            s = permute(t, ((2, 1), (3, 4)))
            @test norm(s) ≈ norm(t)
            @test braid(braid(t, ((2, 1), (3, 4)), (1, 2, 3, 4)), ((2, 1), (3, 4)), (2, 1, 3, 4)) ≈ t
            @test TensorKit.twist(TensorKit.twist(t, 1), 1; inv = true) ≈ t
            # Yang-Baxter, as two ways to reverse three strands
            x = randn(ComplexF64, V ⊗ V ⊗ V ← W)
            lv = (1, 2, 3, 4)
            b1 = braid(braid(braid(x, ((2, 1, 3), (4,)), lv), ((1, 3, 2), (4,)), (2, 1, 3, 4)), ((2, 1, 3), (4,)), (2, 3, 1, 4))
            b2 = braid(braid(braid(x, ((1, 3, 2), (4,)), lv), ((2, 1, 3), (4,)), (1, 3, 2, 4)), ((1, 3, 2), (4,)), (3, 1, 2, 4))
            @test b1 ≈ b2
            # factorizations reconstruct
            m = t
            U, Σ, Vh = svd_compact(m)
            @test U * Σ * Vh ≈ m && U' * U ≈ id(domain(U)) && Vh * Vh' ≈ id(codomain(Vh))
            Q, R = qr_compact(m)
            @test Q * R ≈ m && Q' * Q ≈ id(domain(Q))
            h = m' * m
            D, P = eigh_full(h)
            @test P * D * P' ≈ h && all(real(x) >= -1e-10 for (_, b) in blocks(D) for x in diag(b))
            # @planar agrees with composition
            A = randn(ComplexF64, V ⊗ V ← V); B = randn(ComplexF64, V ← V ⊗ V)
            @planar C1[-1 -2; -3 -4] := A[-1 -2; 1] * B[1; -3 -4]
            @test C1 ≈ A * B
        end
        # planar tensors can stay real, since the F-symbols are
        I = SU2k(3); V = Vect[I](0 => 2, 1//2 => 2, 1 => 2)
        t = randn(Float64, V ⊗ V ← V)
        @test scalartype(repartition(t, 1, 2)) === Float64 && norm(repartition(t, 1, 2)) ≈ norm(t)
        # large levels do not put the whole sector in the type of a space
        for k in (1000,)
            L = SU2k(k); V = Vect[L](0 => 1, 1//2 => 1, k // 2 => 1)
            t = randn(ComplexF64, V ⊗ V ← V)
            @test dim(V) ≈ 1 + dim(L(1//2)) + 1 && norm(repartition(t, 1, 2)) ≈ norm(t)
            @test braid(braid(t, ((2, 1), (3,)), (1, 2, 3)), ((2, 1), (3,)), (2, 1, 3)) ≈ t
        end
    end

    @testset "tensors at real q" begin
        I = SU2q(0.8); V = Vect[I](0 => 2, 1//2 => 2, 1 => 1, 5//2 => 1)
        t = randn(Float64, V ⊗ V ← V ⊗ V)
        @test scalartype(repartition(t, 1, 3)) === Float64
        @test norm(repartition(t, 1, 3)) ≈ norm(t) && repartition(repartition(t, 1, 3), 2, 2) ≈ t
        U, Σ, Vh = svd_compact(t)
        @test U * Σ * Vh ≈ t
        @test_throws Exception permute(t, ((2, 1), (3, 4)))
    end

    @testset "threads" begin
        I = SU2k(23); s = [I(j) for j in 0:1//2:5//2]
        labels = [(a, b, c, d, e, f) for a in s for b in s for c in s for d in ⊗(a, b, c) for e in a ⊗ b for f in b ⊗ c if d in e ⊗ c && d in a ⊗ f]
        serial = [Fsymbol(l...) for l in labels]
        QR.empty_caches!()
        parallel = Vector{Float64}(undef, length(labels))
        Threads.@threads for i in eachindex(labels)
            parallel[i] = Fsymbol(labels[i]...)
        end
        @test parallel == serial
    end

    @testset "spectra do not depend on the gauge" begin
        # a chain of τ with the projector on the trivial channel, against TensorKit's own Fibonacci sector
        function chain(τ, L)
            P = Vect[typeof(τ)](τ => 1)
            h = zeros(ComplexF64, P ⊗ P ← P ⊗ P)
            for (c, b) in blocks(h)
                b .= isunit(c) ? 1.0 : 0.0
            end
            H = zeros(ComplexF64, P^L ← P^L)
            for i in 1:L-1
                left = i == 1 ? nothing : id(P^(i - 1)); right = i == L - 1 ? nothing : id(P^(L - i - 1))
                term = h
                left === nothing || (term = left ⊗ term)
                right === nothing || (term = term ⊗ right)
                H = H + term
            end
            return sort!(vcat([eigvals(Hermitian(Matrix(b))) for (_, b) in blocks(H)]...))
        end
        for L in (3, 5)
            @test chain(SU2k(3)(1), L) ≈ chain(FibonacciAnyon(:τ), L)
        end
    end

    @testset "spin-1/2 chains are U_q-invariant XXZ chains" begin
        # every eigenvalue of the anyon chain is one of the XXZ chain, and the ground states agree, at a level and at real q
        for (I, q) in ((SU2k(3), cispi(1 / 5)), (SU2q(0.8), complex(0.8)))
            for L in (5,)
                a = spin_half_chain(I, L); x = xxz_chain(q, L)
                @test all(λ -> minimum(abs.(x .- λ)) < 1e-10, a)
                @test minimum(a) ≈ minimum(real, x) atol = 1e-10
            end
        end
    end

    @testset "the tutorial" begin
        I = SU2k(3)
        a = I(1)
        @test collect(a ⊗ a) == [I(0), I(1)] && dim(a) ≈ (1 + √5) / 2 && Fsymbol(a, a, a, a, a, a) ≈ -(√5 - 1) / 2
        @test abs(Rsymbol(a, a, I(0))) ≈ 1
        V = Vect[I](0 => 2, 1//2 => 2, 1 => 1)
        t = randn(ComplexF64, V ⊗ V ← V ⊗ V)
        @test norm(repartition(t, 1, 3)) ≈ norm(t) && norm(braid(t, ((2, 1), (3, 4)), (1, 2, 3, 4))) ≈ norm(t)
        U, S, Vh = svd_compact(t)
        @test U * S * Vh ≈ t
        @test collect(SU2k(0)(0) ⊗ SU2k(0)(0)) == [SU2k(0)(0)]
    end
end
