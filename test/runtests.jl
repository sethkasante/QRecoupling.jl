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
        m = qdim(1/2)
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
        @test q6j(1, 1, 3, 1, 1, 1).base.sign == 0
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
        @test QR.rmatrix_mono(2, 2, 2) == rmatrix(1, 1, 1)
        d = q6j(1, 1, 1, 1, 1, 1)
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

    @testset "Large spins stay finite" begin
        val = q6j(50.0, 50.0, 50.0, 50.0, 50.0, 50.0, k=2000)
        @test isfinite(val)
        @test abs(val) < 1.0
    end
end
