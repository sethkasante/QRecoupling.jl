# Run from the repository root:
#   julia --project=. benchmark/qcg_verification.jl
#
# A small exact-arithmetic reference, not a replacement for the numerical kernels.
# It demonstrates exact level values for the weighted qcg/q3j rules and checks
# their defining coproduct identities. All algebraic checks use ==, without a
# floating-point tolerance. Only comparisons with the public numerical API use ≈.
# The internal rule constructors used here may change between package versions.

module QCGVerification

using QRecoupling, Nemo, Test
const QR = QRecoupling

"Evaluate a weighted factorial rule exactly at q = exp(iπ/(k+2))."
function weighted_exact_level(s::FactorialSum, w::Int, e2::Int, k::Integer)
    k >= 0 || throw(DomainError(k, "level must be nonnegative"))
    QR.is_empty_sum(s) && return QQBar(0)
    h = Int(k) + 2
    N = QR.max_argument(s)
    N < h || throw(ArgumentError("this reference requires all factorial arguments below k+2"))
    r = root_of_unity(QQBar, 4h)            # exactly q^(1/2), with the chosen embedding
    q = r^2
    facts = [QQBar(1) for _ in 0:N]
    den = q - inv(q)
    for n in 1:N
        facts[n+1] = facts[n] * (q^n - q^(-n)) / den
    end
    pre = QQBar(1)
    for (n, c) in s.pre
        pre *= facts[Int(n)+1]^Int(c)
    end
    # Every factorial here is positive in the unitary embedding, so this root
    # agrees with the public level kernel. This is not a generic-complex branch rule.
    pre = s.sqrt_pre ? sqrt(pre) : pre
    total = QQBar(0)
    for z in s.zlo:s.zhi
        term = q^(w*z)
        for f in s.fac
            term *= facts[Int(f.a)*z+Int(f.b)+1]^Int(f.c)
        end
        total += s.alternating && isodd(z) ? -term : term
    end
    return Int(s.sign0) * r^e2 * pre * total
end

"Exact algebraic level value of qcg; standalone reference, not a public package method."
function qcg_exact_level(k::Integer, j1, m1, j2, m2, j, m=m1+m2)
    k >= 0 || throw(DomainError(k, "level must be nonnegative"))
    J1, M1, J2, M2, J, M = QR.doubled(j1, m1, j2, m2, j, m)
    QR._qδ(J1, J2, J, Int(k)) || return QQBar(0)
    s, (w, e2) = QR.qcg_rule(J1, M1, J2, M2, J, M)
    return weighted_exact_level(s, w, e2, k)
end

"Exact algebraic level value of q3j in the package's CG-derived convention."
function q3j_exact_level(k::Integer, j1, j2, j3, m1, m2, m3=-m1-m2)
    k >= 0 || throw(DomainError(k, "level must be nonnegative"))
    J1, J2, J3, M1, M2, M3 = QR.doubled(j1, j2, j3, m1, m2, m3)
    QR._qδ(J1, J2, J3, Int(k)) || return QQBar(0)
    s, (w, e2) = QR.q3j_rule(J1, J2, J3, M1, M2, M3)
    return weighted_exact_level(s, w, e2, k)
end

# These generators are built directly from the representation, independently of
# the factorial rule. At a unitary level all their nonzero q-numbers are positive.
qn(n, q) = (q^n - q^(-n)) / (q - inv(q))
raise(j, m, q) = abs(m) > j || m == j ? QQBar(0) : sqrt(qn(Int(j-m), q) * qn(Int(j+m+1), q))
lower(j, m, q) = abs(m) > j || m == -j ? QQBar(0) : sqrt(qn(Int(j+m), q) * qn(Int(j-m+1), q))

function run_checks()
    @testset "symbolic spin-half singlet" begin
        # q = r^2. Remove the common normalisation sqrt([2]); the two
        # coefficients are r and -1/r, so the coproduct is rational in r.
        R, t = polynomial_ring(QQ, "r")
        F = fraction_field(R); r = F(t)
        a, b = r, -inv(r)
        @test a / r + b * r == 0          # Δ(E)|0,0> = Δ(F)|0,0> = 0
        @test a^2 + b^2 == r^2 + r^(-2)  # divide by [2] for transpose norm one
    end

    worst = 0.0
    @testset "exact weighted values and coproduct" begin
        # k=2, j1=j2=1 retains only j=0: also test a fusion-truncated sector.
        for k in (2, 3, 4), (j1, j2) in ((1//2, 1//2), (1, 1//2), (1, 1))
            r = root_of_unity(QQBar, 4(k+2)); q = r^2
            # Memoise only within this small representation pair.
            values = Dict{Tuple,typeof(QQBar(0))}()
            cg(x, y, j, m) = get!(values, (x, y, j, m)) do
                qcg_exact_level(k, j1, x, j2, y, j, m)
            end
            for j in abs(j1-j2):min(j1+j2, k-j1-j2), m in -j:j
                norm = QQBar(0)
                for m1 in -j1:j1, m2 in -j2:j2
                    c = cg(m1, m2, j, m)
                    norm += c^2
                    @test r^Int(2(m1+m2)) * c == r^Int(2m) * c  # Δ(K) C = C K
                    # Coefficient of |m1,m2> after raising/lowering a coupled state.
                    e = raise(j1, m1-1, q) * r^Int(2m2) * cg(m1-1, m2, j, m) +
                        r^Int(-2m1) * raise(j2, m2-1, q) * cg(m1, m2-1, j, m)
                    f = lower(j1, m1+1, q) * r^Int(2m2) * cg(m1+1, m2, j, m) +
                        r^Int(-2m1) * lower(j2, m2+1, q) * cg(m1, m2+1, j, m)
                    @test e == raise(j, m, q) * cg(m1, m2, j, m+1)
                    @test f == lower(j, m, q) * cg(m1, m2, j, m-1)
                    a = qcg(Level(k), j1, m1, j2, m2, j, m)
                    worst = max(worst, abs(a - ComplexF64(c)) / max(1.0, abs(ComplexF64(c))))
                    @test a ≈ ComplexF64(c) atol = 2e-14 rtol = 2e-14
                    # q3j enforces m1+m2=-m3 through its own rule.
                    v = q3j_exact_level(k, j1, j2, j, m1, m2, -m)
                    @test v * sqrt(qn(Int(2j+1), q)) == (-1)^Int(j1-j2+m) * c
                    @test q3j(Level(k), j1, j2, j, m1, m2, -m) ≈ ComplexF64(v) atol = 2e-14 rtol = 2e-14
                end
                @test norm == 1  # bilinear norm: no conjugation at complex q
            end
        end
        @test iszero(qcg_exact_level(2, 1, 1, 1, -1, 1)) # forbidden by fusion
        @test_throws DomainError qcg_exact_level(-1, 1, 1, 1, -1, 0)
    end
    println("Worst numerical/exact scaled difference: ", worst)
    println("All coproduct and norm checks above are exact algebraic equalities.")
    println("Example at k=3: ", qcg_exact_level(3, 1//2, 1//2, 1//2, -1//2, 0))
    return worst
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    QCGVerification.run_checks()
end
