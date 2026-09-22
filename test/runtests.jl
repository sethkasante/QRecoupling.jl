using Test
using QRecoupling
using Nemo
using Random

const QR = QRecoupling

# ---------------------------------------------------------------------------
#  Independent references: Racah formulas in BigInt / BigFloat (doubled spins)
# ---------------------------------------------------------------------------

function racah_bounds(J)
    J1, J2, J3, J4, J5, J6 = J
    α = ((J1 + J2 + J3) ÷ 2, (J1 + J5 + J6) ÷ 2, (J2 + J4 + J6) ÷ 2, (J3 + J4 + J5) ÷ 2)
    β = ((J1 + J2 + J4 + J5) ÷ 2, (J1 + J3 + J4 + J6) ÷ 2, (J2 + J3 + J5 + J6) ÷ 2)
    return α, β
end

"Classical 6j symbol from the Racah formula with exact rationals."
function sixj_ref_classical(J)
    QR._δtet(J...) || return 0.0
    f(n) = factorial(big(n))
    J1, J2, J3, J4, J5, J6 = J
    Δ(a, b, c) = f((a + b - c) ÷ 2) * f((a - b + c) ÷ 2) * f((-a + b + c) ÷ 2) // f((a + b + c) ÷ 2 + 1)
    α, β = racah_bounds(J)
    s = sum((-1)^z * f(z + 1) // (prod(f(z - a) for a in α) * prod(f(b - z) for b in β))
            for z in maximum(α):minimum(β))
    pre = Δ(J1, J2, J3) * Δ(J1, J5, J6) * Δ(J2, J4, J6) * Δ(J3, J4, J5)
    return setprecision(() -> Float64(sqrt(BigFloat(pre)) * BigFloat(s)), BigFloat, 256)
end

"q-6j symbol at level k from the Racah formula with signed [n] = sin(nπ/h)/sin(π/h) in BigFloat."
function sixj_ref_level(J, k)
    QR._qδtet(J..., k) || return 0.0
    setprecision(BigFloat, 256) do
        h = k + 2
        qn(n) = sin(n * big(pi) / h) / sin(big(pi) / h)
        qf(n) = n <= 0 ? big(1.0) : prod(qn(i) for i in 1:n)
        J1, J2, J3, J4, J5, J6 = J
        Δ(a, b, c) = qf((a + b - c) ÷ 2) * qf((a - b + c) ÷ 2) * qf((-a + b + c) ÷ 2) / qf((a + b + c) ÷ 2 + 1)
        α, β = racah_bounds(J)
        s = sum((-1)^z * qf(z + 1) / (prod(qf(z - a) for a in α) * prod(qf(b - z) for b in β))
                for z in maximum(α):minimum(β))
        pre = Δ(J1, J2, J3) * Δ(J1, J5, J6) * Δ(J2, J4, J6) * Δ(J3, J4, J5)
        Float64(sqrt(pre) * s)
    end
end

"The same level reference, returned as BigFloat at the caller's precision."
function sixj_ref_level_big(J, k)
    QR._qδtet(J..., k) || return big(0.0)
    h = k + 2
    qn(n) = sin(n * big(pi) / h) / sin(big(pi) / h)
    qf(n) = n <= 0 ? big(1.0) : prod(qn(i) for i in 1:n)
    J1, J2, J3, J4, J5, J6 = J
    Δ(a, b, c) = qf((a + b - c) ÷ 2) * qf((a - b + c) ÷ 2) * qf((-a + b + c) ÷ 2) / qf((a + b + c) ÷ 2 + 1)
    α, β = racah_bounds(J)
    s = sum((-1)^z * qf(z + 1) / (prod(qf(z - a) for a in α) * prod(qf(b - z) for b in β))
            for z in maximum(α):minimum(β))
    return sqrt(Δ(J1, J2, J3) * Δ(J1, J5, J6) * Δ(J2, J4, J6) * Δ(J3, J4, J5)) * s
end

"Classical Wigner 3j symbol (doubled labels)."
function threej_ref_classical(J1, J2, J3, M1, M2, M3)
    (!QR._δ(J1, J2, J3) || !QR._mproj_ok(J1, J2, J3, M1, M2, M3)) && return 0.0
    f(n) = factorial(big(n))
    a1, a2, a3 = (J1 + J2 - J3) ÷ 2, (J1 - J2 + J3) ÷ 2, (-J1 + J2 + J3) ÷ 2
    Δ = f(a1) * f(a2) * f(a3) // f((J1 + J2 + J3) ÷ 2 + 1)
    pre = Δ * f((J1 + M1) ÷ 2) * f((J1 - M1) ÷ 2) * f((J2 + M2) ÷ 2) * f((J2 - M2) ÷ 2) *
          f((J3 + M3) ÷ 2) * f((J3 - M3) ÷ 2)
    zmin = max(0, (J2 - J3 - M1) ÷ 2, (J1 - J3 + M2) ÷ 2)
    zmax = min(a1, (J1 - M1) ÷ 2, (J2 + M2) ÷ 2)
    s = sum((-1)^z // (f(z) * f(a1 - z) * f((J1 - M1) ÷ 2 - z) * f((J2 + M2) ÷ 2 - z) *
                       f((J3 - J2 + M1) ÷ 2 + z) * f((J3 - J1 - M2) ÷ 2 + z)) for z in zmin:zmax)
    ph = iseven((J1 - J2 - M3) ÷ 2) ? 1 : -1
    return setprecision(() -> Float64(ph * sqrt(BigFloat(pre)) * BigFloat(s)), BigFloat, 256)
end

"Random admissible doubled 6j labels (classical when k === nothing)."
function rand6j(rng, Jmax, k)
    while true
        J = Tuple(rand(rng, 0:Jmax, 6))
        (k === nothing ? QR._δtet(J...) : QR._qδtet(J..., k)) && return J
    end
end

spins(J) = J ./ 2
close(a, b; tol = 1e-11) = isapprox(a, b; atol = tol, rtol = tol)


@testset "QRecoupling.jl" begin

    @testset "Quantum dimensions and monomials" begin
        @test qdim(1/2, k=2) ≈ sqrt(2.0) atol=1e-14
        m = qdim(Symbolic(), 1/2)
        @test m.q_pow == -1
        @test m.phi_exps == [2 => 1]
        @test rmatrix(0.5, 0.5, 0.0, k=1) ≈ -cispi(-1/2)
        # [2j+1] keeps its sign above the level: [n] = sin(nπ/h)/sin(π/h)
        for k in (1, 4, 7), J in 0:3(k+2)
            @test qdim(J/2, k=k) ≈ sin((J + 1) * π / (k + 2)) / sin(π / (k + 2)) atol=1e-12
        end
    end

    @testset "Admissibility and input validation" begin
        @test q6j(1, 1, 3, 1, 1, 1, k=10) == 0.0
        @test q6j(1, 1, 1, 1, 1, 1, k=2) == 0.0
        @test q6j(Symbolic(), 1, 1, 3, 1, 1, 1).base.sign == 0
        @test iszero(q6j(1, 1, 1, 1, 1, 1, k=2, exact=true))
        @test_throws ArgumentError q6j(0.3, 0.3, 0.6, 0.3, 0.3, 0.6, k=10)
        @test_throws ArgumentError qdim(1//3, k=5)
        @test q6j(big(1), Int32(1), 1//1, 1.0f0, 1, 1, k=10) ≈ q6j(1, 1, 1, 1, 1, 1, k=10)
        # |m| > j and parity mismatches vanish on every path
        @test q3j(1, 1, 1, 2, -2, 0, k=10) == 0.0
        @test q3j(1, 1, 1, 2, -2, 0, k=10, eager=true) == 0.0
        @test q3j(1, 1, 1, 2, -2, 0, q=1) == 0.0
        @test q3j(1, 1, 1, 1/2, -1/2, 0, q=1) == 0.0
    end

    @testset "6j symbols agree with independent Racah sums" begin
        rng = MersenneTwister(2026)
        for trial in 1:150
            k = rand(rng, 3:30)
            J = rand6j(rng, min(2k, 16), k)
            js = spins(J)
            ref = sixj_ref_level(J, k)
            @test q6j(js..., k=k) ≈ ref atol=1e-12
            @test q6j(js..., k=k, eager=true) ≈ ref atol=1e-12
            @test real(q6j(js..., q=cispi(1 / (k + 2)))) ≈ ref atol=1e-9
            if trial % 5 == 0
                ex = QR.evaluate_exact(q6j(js..., k=k, exact=true))
                @test real(ex) ≈ ref atol=1e-12
                @test abs(imag(ex)) < 1e-12
            end
        end
        for trial in 1:100
            J = rand6j(rng, 16, nothing)
            js = spins(J)
            ref = sixj_ref_classical(J)
            @test q6j(js..., q=1) ≈ ref atol=1e-12
            @test Float64(q6j(js..., q=1, exact=true)) ≈ ref atol=1e-12
        end
    end

    @testset "F and G symbols" begin
        rng = MersenneTwister(7)
        for trial in 1:40
            k = rand(rng, 3:20)
            J = rand6j(rng, min(2k, 12), k)
            js = spins(J)
            dims = [qdim(j, k=k) for j in js]
            six = sixj_ref_level(J, k)
            g = sqrt(prod(dims)) * six
            @test close(gsymbol(js..., k=k), g)
            @test close(real(QR.evaluate_exact(gsymbol(js..., k=k, exact=true))), g; tol = 1e-10)
            phase = iseven((J[1] + J[2] + J[4] + J[5]) ÷ 2) ? 1 : -1
            @test close(fsymbol(js..., k=k), phase * sqrt(dims[3] * dims[6]) * six)
        end
    end

    @testset "3j symbols" begin
        rng = MersenneTwister(3)
        n = 0
        while n < 200
            J1, J2 = rand(rng, 0:8), rand(rng, 0:8)
            J3 = rand(rng, abs(J1 - J2):2:J1+J2)
            M1 = rand(rng, -J1:2:J1)
            M2 = rand(rng, -J2:2:J2)
            M3 = -M1 - M2
            abs(M3) > J3 && continue
            n += 1
            a = (J1, J2, J3, M1, M2, M3) ./ 2
            @test q3j(a..., q=1) ≈ threej_ref_classical(J1, J2, J3, M1, M2, M3) atol=1e-12
            @test q3j(a..., k=30) ≈ q3j(a..., k=30, eager=true) atol=1e-12
        end
    end

    @testset "Phases and cyclotomic values at roots of unity" begin
        q = cispi(1/6)
        powers(z) = CyclotomicMonomial(1, z, Pair{Int,Int}[], 0)
        s = qseries(powers, 0:2)
        @test qeval(s, k=4) ≈ 1 + q + q^2
        @test qeval(s, q=q) ≈ 1 + q + q^2

        Zx, x = polynomial_ring(ZZ, "x")
        for h in (3, 5, 6, 7, 12), d in 1:3h
            d == h && continue
            k = h - 2
            c = cyclotomic(d, x)
            y = cispi(2 / h)
            tv = sum(Int(coeff(c, i)) * y^i for i in 0:degree(c))
            m = CyclotomicMonomial(1, 0, [d => 1], d)
            @test project_discrete(m, k) ≈ tv atol=1e-10
            ex = QR.CompositeExactResult(k, QR.ONE_MONOMIAL, project_exact(m, k))
            @test QR.evaluate_exact(ex) ≈ tv atol=1e-10
        end
        # Φ_h vanishes, and its inverse is a pole
        @test project_discrete(CyclotomicMonomial(1, 0, [7 => 1], 7), 5) == 0
        @test_throws DomainError project_discrete(CyclotomicMonomial(1, 0, [7 => -1], 7), 5)
    end

    @testset "Exact arithmetic and phases" begin
        r = q6j(1, 1, 1, 1, 1, 1, k=5, exact=true)
        v = QR.evaluate_exact(r)
        @test QR.evaluate_exact(QPhase(Int8(-1), 3//1) * r) ≈ -cispi(3/7) * v atol=1e-14
        @test_throws ArgumentError QPhase(Int8(1), 1//2) * r
        @test QR.rmatrix_mono(2, 2, 2) isa QPhase
        @test QR.rmatrix_mono(2, 2, 2) == rmatrix(Symbolic(), 1, 1, 1)
        d = q6j(Symbolic(), 1, 1, 1, 1, 1, 1)
        @test qeval(QR.fuse_root(d, qint(3)), k=10) ≈ qeval(d, k=10) * qeval(qint(3), k=10)

        # exact orthogonality in ℚ(ζ): Σ_x [2x+1] {1 1 x; 1 1 1}^2 = 1/[3]
        k = 5
        lhs = sum(qdim(x, k=k, exact=true) * q6j(1, 1, x, 1, 1, 1, k=k, exact=true)^2 for x in 0:2)
        @test iszero(lhs - 1 / qdim(1, k=k, exact=true))
    end

    @testset "Series construction" begin
        interior(z) = z == 1 ? QR.ZERO_MONOMIAL : qint(z + 2)
        @test_throws ArgumentError qseries(interior, 0:3)
        tail = qseries(z -> z >= 2 ? QR.ZERO_MONOMIAL : qint(z + 2), 0:3)
        @test qeval(tail, q=0.7) ≈ qeval(qint(2), q=0.7) + qeval(qint(3), q=0.7)
        lead = qseries(z -> z == 0 ? QR.ZERO_MONOMIAL : qint(z + 2), 0:3)
        @test qeval(lead, q=0.7) ≈ sum(qeval(qint(n), q=0.7) for n in 3:5)
        buf = QR.CycloBuffer(8)
        buf.sign = 0
        dcr = build_dcr!(buf, b -> nothing, (b, z) -> add_qfact!(b, z), (b, z) -> add_qint!(b, z + 1), 1, 3)
        @test qeval(dcr, q=1) ≈ 9.0
    end

    @testset "Caches" begin
        v = q6j(2, 2, 2, 2, 2, 2, k=10)
        q6j(2, 2, 2, 2, 2, 2, k=10, exact=true)
        @test empty_caches!() === nothing
        @test q6j(2, 2, 2, 2, 2, 2, k=10) == v
    end

    @testset "Orthogonality" begin
        # Σ_x [2x+1][2p+1] {a b x; c d p}{a b x; c d p'} = δ_{pp'}
        for (A, B, C, D) in ((2, 2, 2, 2), (3, 2, 4, 3), (4, 6, 5, 3)), k in (nothing, 6, 12)
            xs = max(abs(A - B), abs(C - D)):2:min(A + B, C + D)
            ps = max(abs(A - D), abs(C - B)):2:min(A + D, C + B)
            kw = k === nothing ? (q = 1,) : (k = k,)
            for P1 in ps, P2 in ps
                valid = k === nothing || (QR._qδ(A, D, P1, k) && QR._qδ(C, B, P1, k))
                val = sum(qdim(X / 2; kw...) * qdim(P1 / 2; kw...) *
                          q6j(A / 2, B / 2, X / 2, C / 2, D / 2, P1 / 2; kw...) *
                          q6j(A / 2, B / 2, X / 2, C / 2, D / 2, P2 / 2; kw...) for X in xs)
                @test val ≈ (P1 == P2 && valid ? 1.0 : 0.0) atol=1e-12
            end
        end
    end

    @testset "Biedenharn–Elliott identity" begin
        # Σ_x (-1)^{S+x} [2x+1] {a b x; c d p}{c d x; e f q}{e f x; b a r} = {p q r; e a d}{p q r; f b c}
        # labels are drawn so that both 6j symbols on the right-hand side are admissible
        rng = MersenneTwister(5)
        checked = 0
        while checked < 150
            L = rand(rng, 0:6, 9)
            A, B, C, D, E, F, P, Q, R = L
            (QR._δtet(P, Q, R, E, A, D) && QR._δtet(P, Q, R, F, B, C)) || continue
            checked += 1
            for kw in ((q = 1,), (k = 9,))
                rhs = q6j(P / 2, Q / 2, R / 2, E / 2, A / 2, D / 2; kw...) *
                      q6j(P / 2, Q / 2, R / 2, F / 2, B / 2, C / 2; kw...)
                lhs = 0.0
                for X in 0:16
                    t = q6j(A / 2, B / 2, X / 2, C / 2, D / 2, P / 2; kw...) *
                        q6j(C / 2, D / 2, X / 2, E / 2, F / 2, Q / 2; kw...) *
                        q6j(E / 2, F / 2, X / 2, B / 2, A / 2, R / 2; kw...)
                    iszero(t) && continue
                    s2 = sum(L) + X
                    @assert iseven(s2)
                    lhs += (iseven(s2 ÷ 2) ? 1 : -1) * qdim(X / 2; kw...) * t
                end
                @test lhs ≈ rhs atol=1e-12
            end
        end
        @test checked == 150
    end

    @testset "Factorial sums: valuations, zeros and poles" begin
        # closed-form valuations must agree with scanning the expanded DCR exponents
        function dcr_scan(dcr, k)
            dcr.base.sign == 0 && return (:empty, UnitRange{Int}[])
            h = k + 2
            eh(m) = QR._phi_exponent(m, h)
            vp = eh(dcr.radical) + 2eh(dcr.root)
            V = eh(dcr.base)
            z = first(dcr.z_range)
            contrib = Int[]
            vp + 2V < 0 && return (:pole, UnitRange{Int}[])
            vp + 2V == 0 && push!(contrib, z)
            for r in dcr.ratios
                z += 1
                V += eh(r)
                vp + 2V < 0 && return (:pole, UnitRange{Int}[])
                vp + 2V == 0 && push!(contrib, z)
            end
            isempty(contrib) && return (:zero, UnitRange{Int}[])
            return (:finite, QR._merge_ranges([z:z for z in contrib]))
        end
        rng = MersenneTwister(77)
        statuses = Set{Symbol}()
        for _ in 1:600
            J = rand6j(rng, 30, nothing)
            k = rand(rng, 1:30)
            c = QR.classify_at_level(QR.sixj_sum(J...), k)
            push!(statuses, c[1])
            @test c == dcr_scan(QR.q6j_dcr(J...), k)
        end
        @test statuses ⊇ Set([:pole, :zero, :finite])
        n3 = 0
        while n3 < 300
            J1, J2 = rand(rng, 0:20), rand(rng, 0:20)
            J3 = rand(rng, abs(J1 - J2):2:J1+J2)
            M1 = rand(rng, -J1:2:J1); M2 = rand(rng, -J2:2:J2)
            abs(M1 + M2) > J3 && continue
            k = rand(rng, 1:20)
            n3 += 1
            @test QR.classify_at_level(QR.threej_sum(J1, J2, J3, M1, M2), k) == dcr_scan(QR.q3j_dcr(J1, J2, J3, M1, M2), k)
        end

        # word-size modular arithmetic
        p = QR.primes_1mod(44, 1)[1]
        m = QR.Montgomery(p)
        for _ in 1:200
            a, b = rand(rng, UInt64(0):p-1), rand(rng, UInt64(0):p-1)
            @test QR.from_mont(m, QR.mont_mul(m, QR.to_mont(m, a), QR.to_mont(m, b))) == UInt64(mod(big(a) * big(b), big(p)))
        end
        @test all(q -> Nemo.is_prime(ZZ(q)) && (q - 1) % 44 == 0, QR.primes_1mod(44, 5))
        for _ in 1:300
            n = rand(rng, UInt64(1) << 40:UInt64(1) << 62) | UInt64(1)
            @test QR.is_prime_u64(n) == Nemo.is_prime(ZZ(n))
        end
        @test !QR.is_prime_u64(UInt64(3215031751))          # strong pseudoprime to bases 2, 3, 5, 7
        g = QR.root_of_unity(p, 44)
        @test QR._powmod(g, UInt64(44), p) == 1 && all(d -> QR._powmod(g, UInt64(d), p) != 1, (4, 22))

        # symbols that are exactly zero by cancellation come back as zero on both paths
        for (J, k) in (((10, 10, 10, 10, 10, 10), 20), ((30, 30, 30, 30, 30, 30), 60), ((50, 50, 6, 50, 50, 50), 100))
            @test isempty(QR.project_exact(QR.q6j_dcr(J...), k).terms)     # Nemo agrees it is zero
            @test q6j(spins(J)..., k=k) === 0.0
            @test iszero(q6j(spins(J)..., k=k, exact=true))
        end

        # the zero test agrees with Nemo on random level-admissible symbols
        for _ in 1:120
            k = rand(rng, 3:30)
            J = rand6j(rng, min(2k, 20), k)
            truth = isempty(QR.project_exact(QR.q6j_dcr(J...), k).terms)
            @test QR.is_zero_at_level(QR.sixj_sum(J...), k) == truth
        end

        # the level path equals the DCR projection, for 6j, F and G symbols
        for _ in 1:200
            k = rand(rng, 3:40)
            J = rand6j(rng, min(2k, 24), k)
            @test q6j(spins(J)..., k=k) ≈ QR.project_discrete(QR.q6j_dcr(J...), k) atol=1e-12
            @test fsymbol(spins(J)..., k=k) ≈ QR.project_discrete(QR.fsymbol_dcr(J...), k) atol=1e-11
            @test gsymbol(spins(J)..., k=k) ≈ QR.project_discrete(QR.gsymbol_dcr(QR.canonical_spins(spins(J)...)...), k) atol=1e-10 rtol=1e-11
        end

        # heavy cancellation: Float64 loses most digits, so the sum is redone at higher precision
        for (J, k) in (((100, 100, 100, 100, 100, 100), 200), ((200, 200, 200, 200, 200, 200), 500), ((160, 150, 140, 130, 170, 180), 400))
            r = sixj_ref_level(J, k)
            @test isapprox(q6j(spins(J)..., k=k), r; rtol=1e-12)
        end

        # BigFloat
        J = (40, 40, 40, 40, 40, 40)
        v = setprecision(() -> q6j(spins(J)..., k=80, T=BigFloat), BigFloat, 256)
        r = setprecision(() -> QR.project_discrete(QR.q6j_dcr(J...), 80, BigFloat), BigFloat, 256)
        @test v isa BigFloat
        @test abs(v - r) <= big(10)^-60 * abs(r)
    end

    @testset "K-word precision tier" begin
        # arithmetic against 512-bit BigFloat
        setprecision(BigFloat, 512) do
            for K in (3, 4)
                x = QR.MW{K}(big(pi)); y = QR.MW{K}(sqrt(big(2)))
                tol = big(2.0)^(-50K)
                @test abs(BigFloat(x * y) / (big(pi) * sqrt(big(2))) - 1) < tol
                @test abs(BigFloat(x / y) / (big(pi) / sqrt(big(2))) - 1) < tol
                @test abs(BigFloat(sqrt(x)) / sqrt(big(pi)) - 1) < tol
                @test abs(BigFloat(x + y) / (big(pi) + sqrt(big(2))) - 1) < tol
                @test abs(BigFloat((x + y) - y) / big(pi) - 1) < tol       # cancellation keeps a leading word
            end
        end
        # symbols beyond the compensated range (κ ≈ 1e17 … 1e43) at full double precision, no BigFloat
        for (J, k) in (((240, 240, 240, 240, 240, 240), 500), ((400, 400, 400, 400, 400, 400), 900),
                       ((600, 600, 600, 600, 600, 600), 1200))
            r = setprecision(() -> sixj_ref_level_big(J, k), BigFloat, 1024)
            @test abs(q6j(spins(J)...; k = k) - r) <= 2.5e-16 * abs(r)
        end
    end

    @testset "Certified bounds and compensated evaluation" begin
        # the certified bounds dominate the true error, for both passes, on typical and hard symbols
        rng = MersenneTwister(44)
        viol = 0; n = 0
        while n < 120
            k = rand(rng, 20:200)
            J = rand6j(rng, min(2k, 120), k)
            s = QR.sixj_sum(J...)
            st, segs = QR.classify_at_level(s, k)
            (st === :finite && QR._within_level_tables(s, segs, k)) || continue
            r = setprecision(() -> sixj_ref_level_big(J, k), BigFloat, 512)
            tab = QR.qint_tables(Float64, k)
            v, _, B, _ = QR._sum_at_level(s, segs, k, tab)
            vc, Bc = QR._sum_compensated(s, segs, k, tab)
            viol += (abs(big(v) - r) > B) + (abs(big(vc) - r) > Bc)
            n += 1
        end
        @test viol == 0

        # compensation gives full double precision where the plain loop loses 6+ digits
        for (J, k) in (((100, 100, 100, 100, 100, 100), 200), ((200, 200, 200, 200, 200, 200), 500),
                       ((160, 150, 140, 130, 170, 180), 400))
            r = setprecision(() -> sixj_ref_level_big(J, k), BigFloat, 512)
            @test abs(q6j(spins(J)...; k = k) - r) <= 2.5e-16 * abs(r)
        end
        # ...and classically, with exact integer ratios
        for J in ((60, 60, 60, 60, 60, 60), (120, 102, 166, 118, 104, 164), (200, 200, 200, 200, 200, 200))
            s = QR.sixj_sum(J...)
            tab = QR.classical_tables(Float64, QR.max_argument(s))
            @test QR._integer_ratios(s, UnitRange{Int}[s.zlo:s.zhi], tab)
            @test isapprox(q6j(spins(J)...; q = 1), sixj_ref_classical(J); rtol = 2.5e-16)
        end
        # the escalation that used to need BigFloat is now a compensated pass (no BigFloat allocations)
        q6j(50, 50, 50, 50, 50, 50; k = 200)
        @test (@allocated q6j(50, 50, 50, 50, 50, 50; k = 200)) < 4096
    end

    @testset "Split-exponent factorial tables" begin
        # every [n]! and 1/[n]! entry is within one rounding of the value computed in 512 bits
        for k in (3, 10, 100, 1000)
            tab = QR.QIntTables(Float64, k)
            h = k + 2
            worst = setprecision(BigFloat, 512) do
                w = big(0.0); f = big(1.0)
                for n in 0:k+1
                    n > 0 && (f *= sin(n * big(pi) / h) / sin(big(pi) / h))
                    w = max(w, abs(ldexp(big(tab.fm[n+1]), tab.fe[n+1]) / f - 1),
                               abs(ldexp(big(tab.gm[n+1]), tab.ge[n+1]) * f - 1))
                end
                w
            end
            @test worst <= 2.0^-52
        end
        # classical tables against exact integer factorials
        tab = QR.classical_tables(Float64, 3000)
        @test all(n -> abs(ldexp(big(tab.fm[n+1]), tab.fe[n+1]) / factorial(big(n)) - 1) <= 2.0^-52, 0:3000)
        # typical accepted values now carry essentially full double precision
        rng = MersenneTwister(21)
        errs = Float64[]
        while length(errs) < 60
            k = rand(rng, 20:60)
            J = rand6j(rng, 40, k)
            r = sixj_ref_level(J, k)
            abs(r) > 1e-200 && push!(errs, abs(q6j(spins(J)...; k = k) - r) / abs(r))
        end
        @test sort(errs)[30] < 5e-15                      # median
    end

    @testset "Classical limit from the factorial rule" begin
        # large spins, where the floating-point projection of the expanded form lost the value
        for J in ((200, 200, 200, 200, 200, 200), (160, 150, 140, 130, 170, 180),
                  (100, 120, 80, 90, 110, 70), (400, 400, 400, 400, 400, 400))
            QR._δtet(J...) || continue
            @test isapprox(q6j(spins(J)...; q = 1), sixj_ref_classical(J); rtol = 1e-10)
        end
        for J in (100, 200)
            @test isapprox(q3j(J / 2, J / 2, J / 2, 0, 0, 0; q = 1), threej_ref_classical(J, J, J, 0, 0, 0); rtol = 1e-10)
        end

        # non-trivial classical zeros: the modular test agrees with exact integer arithmetic, and the
        # value comes back as exactly zero rather than rounding noise
        f(n) = factorial(big(n))
        function exact_zero(J)
            α, β = racah_bounds(J)
            sum((-1)^z * f(z + 1) // (prod(f(z - a) for a in α) * prod(f(b - z) for b in β))
                for z in maximum(α):minimum(β)) == 0
        end
        rng = MersenneTwister(99)
        nzero = 0; agree = 0; n = 0
        while n < 1500
            J = rand6j(rng, 16, nothing)
            n += 1
            z = exact_zero(J)
            agree += (QR.is_classical_zero(QR.sixj_sum(J...)) == z)
            if z
                nzero += 1
                @test q6j(spins(J)...; q = 1) == 0.0
            end
        end
        @test agree == n
        @test nzero > 0

        # F and G symbols take the same route
        J = (6, 4, 4, 6, 4, 6)
        d = [J[i] + 1 for i in 1:6]
        six = sixj_ref_classical(J)
        @test isapprox(gsymbol(spins(J)...; q = 1), sqrt(prod(d)) * six; rtol = 1e-12)
        ph = iseven((J[1] + J[2] + J[4] + J[5]) ÷ 2) ? 1 : -1
        @test isapprox(fsymbol(spins(J)...; q = 1), ph * sqrt(d[3] * d[6]) * six; rtol = 1e-12)
    end

    @testset "Batched API, targets and queries" begin
        L = all_6j(k = 5)
        @test !isempty(L)
        @test all(l -> QR._qδtet(QR.doubled(l...)..., 5), L)

        # a batch equals the per-symbol loop, serially and on workers
        ref = [q6j(l...; k = 5) for l in L]
        @test q6j(L; k = 5) == ref
        @test q6j(L; k = 5, threads = 1) == ref
        @test q6j(L; k = 5, threads = 4) == ref

        # canonical representatives cover the same symmetry classes
        Lc = all_6j(k = 5, canonical = true)
        @test length(Lc) < length(L)
        @test Set(map(l -> QR.canonical_spins(l...), Lc)) == Set(map(l -> QR.canonical_spins(l...), L))

        # one label tuple is accepted as well as a collection
        @test q6j([(1, 1, 1, 1, 1, 1)]; k = 10)[1] == q6j(1, 1, 1, 1, 1, 1; k = 10)
        @test q6j((1, 1, 1, 1, 1, 1); k = 10)[1] == q6j(1, 1, 1, 1, 1, 1; k = 10)
        @test_throws ArgumentError q6j([(1, 1, 1)]; k = 10)

        # level sweeps and label-by-level grids
        ks = 3:12
        @test q6j(1, 1, 1, 1, 1, 1; k = ks) == [q6j(1, 1, 1, 1, 1, 1; k = kk) for kk in ks]
        @test q6j(1, 1, 1, 1, 1, 1; k = ks, exact = true) == [q6j(1, 1, 1, 1, 1, 1; k = kk, exact = true) for kk in ks]
        g = q6j(L[1:4]; k = 4:6)
        @test size(g) == (4, 3)
        @test g[:, 2] == q6j(L[1:4]; k = 5)

        # the other families batch too, including the defaulted m3 of a 3j
        @test q3j([(1, 1, 1, 1, -1, 0), (1, 1, 2, 1, -1, 0)]; k = 10) ==
              [q3j(1, 1, 1, 1, -1, 0; k = 10), q3j(1, 1, 2, 1, -1, 0; k = 10)]
        @test q3j([(1, 1, 1, 1, -1)]; k = 10)[1] == q3j(1, 1, 1, 1, -1; k = 10)
        @test fsymbol(L[1:5]; k = 5) == [fsymbol(l...; k = 5) for l in L[1:5]]
        @test gsymbol(L[1:5]; k = 5) == [gsymbol(l...; k = 5) for l in L[1:5]]
        @test q6j(L[1:5]; k = 5, T = BigFloat) == [q6j(l...; k = 5, T = BigFloat) for l in L[1:5]]

        # heavy cancellation inside a batch: precision escalation has to work on worker tasks
        heavy = [(25, 25, 25, 25, 25, 25), (50, 50, 50, 50, 50, 50), (30, 25, 20, 35, 40, 45)]
        hb = q6j(heavy; k = 200, threads = 4)
        @test hb == [q6j(l...; k = 200) for l in heavy]
        for (n, l) in enumerate(heavy)
            # Float64 keeps about 11 digits once cancellation is accounted for; see `value_at_level`
            @test isapprox(hb[n], sixj_ref_level(Int.(2 .* l), 200); rtol = 1e-10)
        end

        # targets are the keywords
        @test q6j(Level(7), 1, 1, 1, 1, 1, 1) == q6j(1, 1, 1, 1, 1, 1; k = 7)
        @test q6j(Level(7; T = BigFloat), 1, 1, 1, 1, 1, 1) == q6j(1, 1, 1, 1, 1, 1; k = 7, T = BigFloat)
        @test q6j(Exact(7), 1, 1, 1, 1, 1, 1) == q6j(1, 1, 1, 1, 1, 1; k = 7, exact = true)
        @test q6j(At(0.7), 1, 1, 1, 1, 1, 1) == q6j(1, 1, 1, 1, 1, 1; q = 0.7)
        @test q6j(Classical(), 1, 1, 1, 1, 1, 1) == q6j(1, 1, 1, 1, 1, 1; q = 1)
        @test q6j(Classical(exact = true), 1, 1, 1, 1, 1, 1) == q6j(1, 1, 1, 1, 1, 1; q = 1, exact = true)
        @test qdim(Level(10), 1 // 2) == qdim(1 // 2; k = 10)
        @test fsymbol(Level(5), 1, 1, 1, 1, 1, 1) == fsymbol(1, 1, 1, 1, 1, 1; k = 5)
        @test q3j(Level(10), 1, 1, 1, 1, -1, 0) == q3j(1, 1, 1, 1, -1, 0; k = 10)
        @test q6j(Level(3:8), 1, 1, 1, 1, 1, 1) == q6j(1, 1, 1, 1, 1, 1; k = 3:8)
        @test q6j(Level(5), L) == q6j(L; k = 5)
        @test occursin("Level(5)", string(Level(5)))

        # the eager route still works, deprecated
        @test q6j(1, 1, 1, 1, 1, 1; k = 10, eager = true) ≈ q6j(1, 1, 1, 1, 1, 1; k = 10) atol = 1e-12

        # queries agree with the values they predict
        for kk in (5, 9, 14), l in L[1:60]
            @test iszero_at(kk, l...) == (q6j(l...; k = kk) == 0.0)
            @test !issingular_at(kk, l...)
        end
        zs = iszero_at(5, L)
        @test zs == [iszero_at(5, l...) for l in L]
        @test zs == [q6j(l...; k = 5) == 0.0 for l in L]
        @test iszero_at(8, 5, 5, 5, 5, 5, 5)                      # not admissible at k = 8, so zero
        @test iszero_at(20, 5, 5, 5, 5, 5, 5)                      # a cancellation zero
        @test !iszero_at(15, 5, 5, 5, 5, 5, 5)
        @test iszero_at(q3j, 10, 1, 1, 1, 2, -2, 0)                # m outside its range
        # issingular_at asks whether the labels are singular at the level, admissibility aside
        @test issingular_at(12, 10, 10, 10, 10, 10, 10)                # triangle sums 30 > 2k
        @test !issingular_at(30, 10, 10, 10, 10, 10, 10)               # admissible from k = 30
        @test !issingular_at(4, 1, 1, 5, 1, 1, 5)                      # no triangle: nothing to be singular
        @test q6j(10, 10, 10, 10, 10, 10; k = 12) == 0.0           # the value still follows the TV convention
        @test issingular_at(q6j, 12, 10, 10, 10, 10, 10, 10)
        @test any(k -> issingular_at(k, 5, 5, 5, 5, 5, 5), 8:14) && !issingular_at(20, 5, 5, 5, 5, 5, 5)
        @test iszero_at(q6j, 5, L) == zs

        # the level spectrum matches what the symbol does at each level
        spec = level_spectrum(5, 5, 5, 5, 5, 5; k = 8:24)
        @test spec[1:7] == fill(:inadmissible, 7)                   # admissible only from k = 15
        @test spec[13] === :cancels                                 # k = 20
        for (n, kk) in enumerate(8:24)
            @test (spec[n] in (:inadmissible, :zero, :cancels)) == (q6j(5, 5, 5, 5, 5, 5; k = kk) == 0.0)
        end
        @test level_spectrum(q6j, 1, 1, 1, 1, 1, 1; k = 3:6, cancellation = false) ==
              level_spectrum(q6j, 1, 1, 1, 1, 1, 1; k = 3:6, threads = 1, cancellation = false)
    end

    @testset "Multiplicity bounds: finiteness and segment count" begin
        QR = QRecoupling
        rng = MersenneTwister(23)
        nrules = 0
        while nrules < 400
            k = rand(rng, 1:150); h = k + 2
            J = Tuple(rand(rng, 0:min(2k, 100), 6))
            QR._qδtet(J..., k) || continue
            s = QR.sixj_sum(J...)
            QR.is_empty_sum(s) && continue
            nrules += 1
            # the prefactor never has negative multiplicity for admissible labels (α_i ≤ k)
            @test QR.prefactor_valuation(s, h) >= 0
            for z in s.zlo:s.zhi
                E = QR.term_valuation(s, z, h)
                ρ = sum(f.c < 0 ? mod(QR._arg(f, z), h) : 0 for f in s.fac)
                @test E == fld(1 + ρ, h)          # closed form from Σ_r n_r = z
                @test 0 <= E <= 6                  # carry counting over seven arguments
            end
            # therefore no poles, and the fixed-size classification always applies
            @test !issingular_at(k, (J .// 2)...)
            st, segs = QR.classify_at_level(s, k)
            @test st !== :pole && length(segs) <= 8
            @test QR._classify_small(s, k) !== nothing
            @test s.zhi - s.zlo < h                # range ≤ Σα/12 ≤ k/3
        end
        # the :pole status belongs to labels outside the level's admissible set
        @test QR.classify_at_level(QR.sixj_sum(20, 20, 20, 20, 20, 20), 12)[1] === :pole
        @test_throws DomainError QR.value_at_level(QR.sixj_sum(20, 20, 20, 20, 20, 20), 12, Float64;
                                                   fallback = () -> 0.0)
    end

    @testset "Allocation-free single-symbol path" begin
        QR = QRecoupling
        # rules are bits values, and the 6j and G rules do not depend on the labelling
        @test isbits(QR.sixj_sum(2, 2, 2, 2, 2, 2)) && isbits(QR.threej_sum(2, 2, 2, 0, 0))
        @test QR.is_empty_sum(QR.sixj_sum(1, 1, 1, 1, 1, 1)) && QR.is_empty_sum(QR.threej_sum(2, 2, 3, 0, 0))
        perms(t) = (t, (t[2], t[1], t[3], t[5], t[4], t[6]), (t[3], t[2], t[1], t[6], t[5], t[4]),
                    (t[1], t[3], t[2], t[4], t[6], t[5]), (t[4], t[5], t[3], t[1], t[2], t[6]),
                    (t[1], t[5], t[6], t[4], t[2], t[3]), (t[4], t[2], t[6], t[1], t[5], t[3]))
        for js in ((5, 4, 3, 4, 5, 3), (7//2, 3, 5//2, 3, 7//2, 2), (20, 17, 15, 22, 19, 14))
            for p in perms(js)
                @test q6j(p...; k = 60) === q6j(js...; k = 60)
                @test q6j(p...; q = 1) === q6j(js...; q = 1)
                @test gsymbol(p...; k = 60) === gsymbol(js...; k = 60)
            end
        end
        # the allocation-free classification agrees with the general one, poles included
        rng = MersenneTwister(7)
        nchecked = 0
        while nchecked < 3000
            k = rand(rng, 1:90)
            J = Tuple(rand(rng, 0:min(2k, 100), 6))
            s = isodd(nchecked) ? QR.sixj_sum(J...) :
                QR.threej_sum(J[1], J[2], J[3], rand(rng, -J[1]:2:J[1]), rand(rng, -J[2]:2:J[2]))
            QR.is_empty_sum(s) && continue
            a = QR.classify_at_level(s, k)
            b = QR._classify_small(s, k)
            @test b !== nothing && a[1] === b[1] && (a[1] === :pole || a[2] == collect(b[2]))
            nchecked += 1
        end
        # The level kernel itself allocates nothing for a short sum (longer ones allocate only the lazy
        # ratio buffer). The public wrapper is checked with a bound rather than zero: its return type is a
        # union (value, DCR, exact, complex), which Julia 1.10 boxes at ~16 B per call and 1.11+ elides.
        # (referenced through the module, not the local `QR` alias: a call through a non-const local is
        # dynamically dispatched and would allocate for that reason alone)
        sfix = QRecoupling.sixj_sum(QRecoupling.doubled(1, 1, 1, 1, 1, 1)...)
        tabfix = QRecoupling.qint_tables(Float64, 10)
        segfix = (sfix.zlo:sfix.zhi,)
        QRecoupling._certified_value(sfix, segfix, 10, tabfix)
        QRecoupling.level_pass1(sfix, 10, tabfix)
        # The bounds are what the invariant needs: nothing may be allocated per term, so a regression there
        # would cost hundreds of bytes. A few words of boxing on the older compiler are tolerated.
        slack = VERSION >= v"1.11" ? 0 : 64
        @test (@allocated QRecoupling._certified_value(sfix, segfix, 10, tabfix)) <= slack
        @test (@allocated QRecoupling.level_pass1(sfix, 10, tabfix)) <= slack
        f() = q6j(1, 1, 1, 1, 1, 1; k = 10) + q6j(5//2, 2, 3//2, 2, 5//2, 3; k = 20) + q6j(1, 2, 3, 3, 2, 1; k = 9)
        f()
        @test (@allocated f()) <= 3 * slack
    end

    @testset "Families by the three-term recurrence" begin
        QR = QRecoupling
        rel(a, b) = iszero(b) ? abs(a) : Float64(abs(a - b) / abs(b))
        rng = MersenneTwister(19)
        work = QR.ColumnWork()
        # level columns, entry by entry against 256-bit single-symbol values
        ncol = 0
        while ncol < 25
            k = rand(rng, 3:80)
            J2, J3, L1, L2, L3 = rand(rng, 0:min(k, 60), 5)
            X = QR.sixj_column_range(J2, J3, L1, L2, L3, k)
            length(X) >= 2 || continue
            out = zeros(length(X))
            QR.sixj_column!(out, J2, J3, L1, L2, L3, k, QR.qint_tables(Float64, k), work)
            for (i, x2) in enumerate(X)
                ref = setprecision(BigFloat, 256) do
                    q6j(x2 // 2, J2 // 2, J3 // 2, L1 // 2, L2 // 2, L3 // 2; k = k, T = BigFloat)
                end
                @test iszero(ref) ? abs(out[i]) < 1e-25 : rel(out[i], ref) < 4e-16
            end
            ncol += 1
        end
        # classical columns, including spins where single symbols cancel heavily
        for (J2, J3, L1, L2, L3) in ((40, 36, 30, 44, 38), (300, 280, 260, 310, 290), (7, 5, 5, 6, 4))
            X = QR.sixj_column_range(J2, J3, L1, L2, L3, nothing)
            out = zeros(length(X))
            QR.sixj_column!(out, J2, J3, L1, L2, L3, QR.ClassicalQ(), work)
            for i in unique([1, 2, length(X) ÷ 2, length(X)])
                x2 = X[i]
                ref = setprecision(BigFloat, 1024) do
                    q6j(x2 // 2, J2 // 2, J3 // 2, L1 // 2, L2 // 2, L3 // 2; q = 1, T = BigFloat)
                end
                @test rel(out[i], ref) < 4e-16
            end
        end
        # routing: batches with runs agree with single symbols, zeros stay exact zeros
        L = all_6j(k = 9)
        b = q6j(L; k = 9, threads = 1)
        s1 = [q6j(l...; k = 9) for l in L]
        @test all(i -> (b[i] == 0) == (s1[i] == 0) && rel(b[i], s1[i]) < 1e-13, eachindex(L))
        @test q6j(L; k = 9, threads = 4) == b
        F = vec([(3, 5//2, x, 3//2, 2, y) for x in 1//2:1:8, y in 0:8 if QR._qδtet(Int.(2 .* (3, 5//2, x, 3//2, 2, y))..., 20)])
        @test length(F) >= 16
        for f in (q6j, fsymbol, gsymbol)
            bf = f(F; k = 20)
            @test all(i -> rel(bf[i], f(F[i]...; k = 20)) < 1e-13, eachindex(F))
            bc = f(F; q = 1)
            @test all(i -> rel(bc[i], f(F[i]...; q = 1)) < 1e-13, eachindex(F))
        end
        # doubled labels are the same labels, and give the same values
        D = all_6j(k = 9, doubled = true)
        @test D isa QR.DoubledLabels && collect(D) == all_6j(k = 9)
        @test q6j(D; k = 9) == q6j(all_6j(k = 9); k = 9)
        @test fsymbol(D; k = 9) == fsymbol(all_6j(k = 9); k = 9)
        @test q6j(D; q = 1) == q6j(all_6j(k = 9); q = 1)
        @test iszero_at(9, D) == iszero_at(9, all_6j(k = 9))

        # a column where a single symbol needs extended precision: every entry to full precision
        col = [(x, 60, 55, 50, 65, 58) for x in 0:130 if QR._qδtet(Int.(2 .* (x, 60, 55, 50, 65, 58))..., 300)]
        bc = q6j(col; k = 300)
        for i in (1, length(col) ÷ 3, length(col) ÷ 2, length(col))
            ref = setprecision(BigFloat, 512) do
                q6j(col[i]...; k = 300, T = BigFloat)
            end
            @test rel(bc[i], ref) < 4e-16
        end
    end

    @testset "Hard single symbols from their column" begin
        QR = QRecoupling
        rel(a, b) = Float64(abs(a - b) / abs(b))
        # symbols whose Racah sum has lost its digits: the recurrence tier must be used and be accurate
        for (js, k) in (((50, 50, 50, 50, 50, 50), 400), ((60, 55, 50, 45, 65, 70), 300),
                        ((100, 100, 100, 100, 100, 100), 500), ((80, 70, 60, 50, 90, 100), 400))
            Jd = QR.doubled(js...)
            s = QR.sixj_sum(Jd...)
            _, st, segs = QR.level_pass1(s, k, QR.qint_tables(Float64, k))
            v = QR.sixj_entry(Jd, k)
            @test v !== nothing
            r = setprecision(BigFloat, 512) do
                q6j(js...; k = k, T = BigFloat)
            end
            @test rel(v, r) < 1e-15
            @test rel(q6j(js...; k = k), r) < 1e-15          # through the public path
            if st !== :done                                   # and the tier is what answered
                @test rel(QR.level_escalate(s, segs, k, Float64, QR.level_zero_table(k); labels = Jd), r) < 1e-15
            end
        end
        # any of the six positions may run the column, and the choice does not change the value
        Jd = QR.doubled(40, 36, 30, 44, 38, 34)
        vals = [QR.sixj_entry((Jd[σ[1]], Jd[σ[2]], Jd[σ[3]], Jd[σ[4]], Jd[σ[5]], Jd[σ[6]]), 200)
                for σ in QR._TO_FIRST]
        @test all(v -> v !== nothing && rel(v, vals[1]) < 1e-14, vals)
        c, steps = QR._best_column(Jd, QR.LevelQ(QR.qint_tables(Float64, 200), 200))
        @test steps >= 0 && QR.sixj_entry(c, 200) !== nothing
        # an exact zero must never be answered by the recurrence (its error estimate rejects it)
        @test q6j(5, 5, 5, 5, 5, 5; k = 20) == 0.0
        @test QR.sixj_entry(QR.doubled(5, 5, 5, 5, 5, 5), 20) === nothing
        # classical hard symbols too
        for js in ((100, 100, 100, 100, 100, 100), (120, 110, 100, 90, 130, 140))
            r = setprecision(BigFloat, 1024) do
                q6j(js...; q = 1, T = BigFloat)
            end
            @test rel(q6j(js...; q = 1), r) < 1e-15
        end
    end

    # include("factorial_rules.jl")
    # include("api_v04.jl")

    # @testset "Large spins stay finite" begin
    #     val = q6j(50.0, 50.0, 50.0, 50.0, 50.0, 50.0, k=2000)
    #     @test isfinite(val)
    #     @test abs(val) < 1.0
    # end
end
