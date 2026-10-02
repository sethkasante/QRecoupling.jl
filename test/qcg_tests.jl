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
            @test q3j(j1, j2, j3, m1, m2, m3) === QR.q3j_factorial(j1, j2, j3, m1, m2, m3)
            @test q3j(Exact(), j1, j2, j3, m1, m2, m3) == QR.q3j_factorial(Exact(), j1, j2, j3, m1, m2, m3)
            cg = (-1)^Int(j1 - j2 - m3) * sqrt(2j3 + 1) * QR.q3j_factorial(j1, j2, j3, m1, m2, m3)
            @test qcg(j1, m1, j2, m2, j3) ≈ cg atol = 1e-15
        end
        @test qcg(Exact(), 1//2, 1//2, 1//2, -1//2, 0)^2 == 1//2
        @test q3j(1, 1, 1, 1, -1, 0; q = 1) === q3j(1, 1, 1, 1, -1, 0)
    end

    @testset "orthogonality: C Cᵀ = Cᵀ C = 1, and the same through q3j" begin
        # residuals scaled by Σ|C_a||C_b|: at complex q and at a level the entries are not bounded by 1
        function defect(f; q3 = nothing, admissible = j -> true)
            worst = 0.0
            for (j1, j2) in pairs, m in -(j1 + j2):(j1 + j2)
                js = [j for j in abs(j1 - j2):(j1 + j2) if abs(m) <= j && admissible(j1 + j2 + j)]
                isempty(js) && continue
                m1s = [m1 for m1 in -j1:j1 if abs(m - m1) <= j2]
                C = q3 === nothing ? [f(j1, m1, j2, m - m1, j) for j in js, m1 in m1s] :
                    # [2j+1]^{1/2} (−1)^{j1−j2+m} (j1 j2 j; m1 m2 −m)_q is the same coefficient
                    [sqrt(q3(Int(2j + 1))) * (-1)^Int(j1 - j2 + m) * f(j1, j2, j, m1, m - m1, -m) for j in js, m1 in m1s]
                A = abs.(C)
                worst = max(worst, maximum(abs.(C * transpose(C) - _eye(length(js))) ./ (A * transpose(A))))
                # columns: only when every j of the block is present
                length(js) == length([j for j in abs(j1 - j2):(j1 + j2) if abs(m) <= j]) &&
                    (worst = max(worst, maximum(abs.(transpose(C) * C - _eye(length(m1s))) ./ (transpose(A) * A))))
            end
            return worst
        end
        for q in (0.5, 0.8, 0.999, 1.3, 2.0, 0.7 + 0.4im, cispi(0.1373))
            @test defect((a...) -> qcg(a...; q = q)) < 1e-13
            @test defect((a...) -> q3j(a...; q = q); q3 = n -> qint(n; q = q)) < 1e-13
        end
        # at a level: whole blocks inside the fusion rule (2(j₁ + j₂) ≤ k for these pairs) are orthogonal;
        # where the rule cuts a block, the admissible rows stay orthonormal (the columns cannot be complete)
        @test defect((a...) -> qcg(Level(11), a...)) < 1e-13
        @test defect((a...) -> q3j(Level(11), a...); q3 = n -> qint(n; k = 11)) < 1e-13
        @test defect((a...) -> qcg(Level(5), a...); admissible = s -> s <= 5) < 1e-13
        # the substituted formula is not orthogonal: this is why it is no longer `q3j`
        old = (j1, m1, j2, m2, j) -> (-1)^Int(j1 - j2 + m1 + m2) * sqrt(qint(Int(2j + 1); q = 0.8)) *
                                     QR.q3j_factorial(j1, j2, j, m1, m2, -m1 - m2; q = 0.8)
        @test defect(old) > 0.5
    end

    @testset "coassociativity: four coefficients give the F-symbol" begin
        # Σ C(j1m1 j2m2|j12) C(j12 j3m3|JM) C(j2m2 j3m3|j23) C(j1m1 j23|JM) = (−1)^{j1+j2+j3+J} √([2j12+1][2j23+1]) {6j}
        # — an independent path (q6j's Racah sum), and a check that the square-root branches agree
        for (t, ev, fs) in (((; q = 0.8), (a...) -> qcg(a...; q = 0.8), (a...) -> fsymbol(a...; q = 0.8)),
                            ((; q = -0.8), (a...) -> qcg(a...; q = -0.8), (a...) -> fsymbol(a...; q = -0.8)),
                            ((; q = -1.3), (a...) -> qcg(a...; q = -1.3), (a...) -> fsymbol(a...; q = -1.3)),
                            ((; q = 0.7 + 0.4im), (a...) -> qcg(a...; q = 0.7 + 0.4im), (a...) -> fsymbol(a...; q = 0.7 + 0.4im)),
                            ((; k = 6), (a...) -> qcg(Level(6), a...), (a...) -> fsymbol(Level(6), a...)))
            worst = 0.0
            ok(a, b, c) = !haskey(t, :k) || a + b + c <= t.k
            for j1 in 0:1//2:3//2, j2 in 0:1//2:3//2, j3 in 0:1//2:3//2,
                j12 in abs(j1 - j2):(j1 + j2), j23 in abs(j2 - j3):(j2 + j3)
                ok(j1, j2, j12) && ok(j2, j3, j23) || continue
                for J in max(abs(j12 - j3), abs(j1 - j23)):min(j12 + j3, j1 + j23)
                    ok(j12, j3, J) && ok(j1, j23, J) || continue
                    lhs = zero(ComplexF64); mass = 0.0
                    for m1 in -j1:j1, m2 in -j2:j2
                        m3 = J - m1 - m2
                        abs(m3) <= j3 && abs(m1 + m2) <= j12 && abs(m2 + m3) <= j23 || continue
                        x = ev(j1, m1, j2, m2, j12) * ev(j12, m1 + m2, j3, m3, J) *
                            ev(j2, m2, j3, m3, j23) * ev(j1, m1, j23, m2 + m3, J)
                        lhs += x; mass += abs(x)
                    end
                    worst = max(worst, abs(lhs - fs(j1, j2, j12, j3, J, j23)) / mass)
                end
            end
            @test worst < 1e-13
        end
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

    @testset "qcg_matrix" begin
        for kw in ((;), (; q = 0.8), (; q = 1.3), (; q = 0.999), (; k = 12), (; q = 0.7 + 0.4im), (; q = -0.8))
            for (j1, j2) in ((1, 1//2), (2, 3//2), (3, 2))
                C, prod, coup = qcg_matrix(j1, j2; kw...)
                @test size(C) == (length(prod), length(coup))
                @test eltype(C) == (haskey(kw, :k) || !(get(kw, :q, 1.0) isa Real) || get(kw, :q, 1.0) < 0 ? ComplexF64 : Float64)
                worst = 0.0
                for (b, (j, m)) in enumerate(coup), (a, (m1, m2)) in enumerate(prod)
                    r = m1 + m2 == m ? qcg(j1, m1, j2, m2, j, m; kw...) : 0.0
                    worst = max(worst, abs(C[a, b] - r) / max(abs(r), 1.0))
                end
                @test worst < 1e-13
                A = abs.(C)
                @test maximum(abs.(transpose(C) * C - _eye(size(C, 2))) ./ max.(transpose(A) * A, 1.0)) < 1e-13
            end
        end
        # kron order: a column is |j m⟩ in the basis kron(e_{m1}, e_{m2})
        C, prod, coup = qcg_matrix(1, 1//2; q = 0.8)
        @test prod[1] == (1, 1//2) && prod[2] == (1, -1//2) && coup[1] == (1//2, 1//2)
        # sectors, and a level that cuts the block
        S, m1, js = qcg_matrix(3, 2, 1; q = 0.9)
        @test m1 == [3, 2, 1, 0, -1] && js == [1, 2, 3, 4, 5]
        @test S[2, 3] ≈ qcg(3, 2, 2, -1, 3, 1; q = 0.9) atol = 1e-15
        @test size(qcg_matrix(3, 2, 1; k = 5)[1]) == (5, 0)
        @test size(qcg_matrix(1, 1; k = 0)[1]) == (9, 0)
        @test eltype(qcg_matrix(1, 1; q = 0.8, T = BigFloat)[1]) == BigFloat
        @test_throws ArgumentError qcg_matrix(1, 1, 1//2)
        @test_throws ArgumentError qcg_matrix(1, 1; q = 0.8, k = 3)
        # larger blocks against arbitrary precision, sampled
        for kw in ((; q = 0.97), (; k = 60))
            ref(a) = haskey(kw, :k) ? qcg(Level(kw.k; T = BigFloat), a...) :
                     setprecision(() -> qcg(a...; q = big(kw.q), T = BigFloat), BigFloat, 256)
            worst = 0.0; c = 0
            for m in -25:25
                S, m1s, js = qcg_matrix(15, 12, m; kw...)
                for (b, j) in enumerate(js), (a, m1) in enumerate(m1s)
                    (c += 1) % 37 == 0 || continue
                    r = ref((15, m1, 12, m - m1, j))
                    abs(r) > 1e-280 && (worst = max(worst, Float64(abs(S[a, b] - r) / abs(r))))
                end
            end
            @test worst < 1e-14
        end
    end

    @testset "qcg_row: the dual recurrence in j" begin
        for kw in ((;), (; q = 0.8), (; q = 1.3), (; q = 0.999), (; q = 2.0), (; k = 10), (; k = 13), (; q = 0.7 + 0.4im), (; q = -0.8))
            worst = 0.0; nrm = 0.0
            for J1 in 0:6, J2 in 0:6, M1 in -J1:2:J1, M2 in -J2:2:J2
                c, js = qcg_row(J1//2, M1//2, J2//2, M2//2; kw...)
                for (i, j) in enumerate(js)
                    r = qcg(J1//2, M1//2, J2//2, M2//2, j; kw...)
                    worst = max(worst, abs(c[i] - r) / max(abs(r), 1.0))
                end
                (!haskey(kw, :k) && get(kw, :q, 1.0) isa Real && get(kw, :q, 1.0) > 0) && (nrm = max(nrm, abs(sum(abs2, c) - 1)))
            end
            @test worst < 1e-13
            @test nrm < 1e-14
        end
        # large rows against arbitrary precision (the recurrence is where these come from)
        for (l, kw) in (((30, 7, 25, -11), (; q = 0.999)), ((40, -13, 35, 20), (; q = 0.8)), ((30, 7, 25, -11), (; k = 120)))
            c, js = qcg_row(l...; kw...)
            for i in 1:5:length(js)
                a = (l..., js[i])
                r = haskey(kw, :k) ? qcg(Level(kw.k; T = BigFloat), a...) :
                    setprecision(() -> qcg(a...; q = big(kw.q), T = BigFloat), BigFloat, 256)
                abs(r) > 1e-280 && @test abs(c[i] - r) / abs(r) < 1e-14
            end
        end
        # the formulas: s_j and e_j are the entries of Cᵀ diag(q^{2m1}) C, divided by q − q⁻¹
        q = 0.8; C, m1s, js = qcg_matrix(3, 2, 1; q = q)
        Tm = transpose(C) * [i == k ? q^(2m1s[i]) : 0.0 for i in eachindex(m1s), k in eachindex(m1s)] * C
        @test maximum(abs(Tm[a, b]) for a in axes(Tm, 1), b in axes(Tm, 2) if abs(a - b) > 1) < 1e-14
        @test_throws ArgumentError qcg_row(1, 2, 1, 0)
        @test qcg_row(1, 1, 1, 0; k = 1)[1] == ComplexF64[]
    end

    @testset "family precision follows the scalar API" begin
        for kw in ((; k = 10, T = BigFloat), (; q = big"-0.8"),
                   (; q = big"0.8", T = Float64), (; q = 0.8 + 0.2im, T = BigFloat))
            row, js = qcg_row(2, 1, 2, -1; kw...)
            C, ms, cols = qcg_matrix(2, 2, 0; kw...)
            expected = haskey(kw, :q) && kw.q isa Real && kw.q > 0 ? BigFloat : Complex{BigFloat}
            @test eltype(row) === eltype(C) === expected
            refs = [qcg(2, 1, 2, -1, j; kw...) for j in js]
            @test row == refs
            @test C[findfirst(==(1), ms), :] == refs
        end
        @test eltype(qcg_matrix(1, 1; k = 5, T = BigFloat)[1]) === Complex{BigFloat}
    end

    @testset "strong cancellation is not an exact zero" begin
        # Previously both initial widths were labelled noise and this returned 0;
        # its magnitude is about 8.68e33, with log2(condition) about 267.
        args = (600, 1, 600, -2, 600)
        v = setprecision(() -> qcg(args...; k = 12000, T = BigFloat), BigFloat, 128)
        r = setprecision(() -> qcg(args...; k = 12000, T = BigFloat), BigFloat, 512)
        @test !iszero(v)
        @test abs(v - r) / abs(r) < big(2.0)^(-120)
        # A complete geometric sum really is zero at a primitive 14th root.
        # The common inverse factorial also exercises exact denominator handling.
        s = FactorialSum(0:13; factors = ((0, 3, -1),))
        @test QR._weighted_level_zero(s, 1, 5)
        @test iszero(QR._weighted_level_value(s, 1, 0, 5, Float64))
        @test !QR._weighted_level_zero(FactorialSum(0:12), 1, 5)
    end

    @testset "unit circle: accurate tables, same values" begin
        for q in (cispi(0.1373), cispi(0.41), 1.02 * cispi(0.3))
            for (f, l) in ((q6j, (7, 6, 5, 4, 6, 7)), (QR.q3j_factorial, (9, 7, 5, 3, -2, -1)), (q3j, (9, 7, 5, 3, -2, -1)))
                r = setprecision(() -> f(l...; q = Complex{BigFloat}(q), T = BigFloat), BigFloat, 320)
                @test abs(f(l...; q = q) - r) / abs(r) < 1e-14
            end
            @test QR._analytic_table(q, 50).accurate
        end
        @test !QR._analytic_table(0.7 + 0.4im, 50).accurate        # off the circle the recurrence is kept
    end

    @testset "near-edge tier and real-q accuracy at large spin" begin
        # labels where the direct sum cancels: the tier answers from the column
        for (l, kw) in (((109//2, -9//2, 229//2, -17//2, 143), (; q = 0.999)), ((56, -13, 52, 0, 68), (; k = 400)))
            J = Int.(2 .* l)
            c = QR._cg_entry(get(kw, :q, nothing), get(kw, :k, nothing), J...)
            r = haskey(kw, :k) ? qcg(Level(kw.k; T = BigFloat), l...) :
                setprecision(() -> qcg(l...; q = big(kw.q), T = BigFloat), BigFloat, 256)
            @test c !== nothing && abs(c - r) / abs(r) < 1e-15
            @test abs(qcg(l...; kw...) - r) / abs(r) < 2.0^-40
        end
        # the real-q table is accurate to half an ulp per entry (it drifted by 10⁴ u near q = 1 before)
        for q in (0.8, 0.999, 1.3, -0.8)
            tab = QR._analytic_table(q, 400)
            worst = setprecision(BigFloat, 256) do
                qb = big(q); fb = big(1.0); w = 0.0
                for n in 1:400
                    fb *= n == 1 ? big(1.0) : (qb^n - qb^(-n)) / (qb - inv(qb))
                    a = tab.facts[n+1]
                    w = max(w, Float64(abs(big(a.m) * big(2.0)^a.e - fb) / abs(fb)))
                end
                w
            end
            @test worst <= eps()
        end
        for l in ((249//2, -43//2, 127, 110, 467//2), (115, -65, 231//2, 219//2, 397//2))
            r = setprecision(() -> qcg(l...; q = big(0.999), T = BigFloat), BigFloat, 256)
            @test abs(qcg(l...; q = 0.999) - r) / abs(r) < 2.0^-40
        end
    end

    @testset "CG workspaces reuse and invalidate sector coefficients" begin
        work = QR.EvaluationWorkspace()
        cases = (((109//2,-9//2,229//2,-17//2,143),(;q=0.999)),
                 ((109//2,-7//2,229//2,-19//2,142),(;q=0.999)), # same sector, different entry
                 ((56,-13,52,0,68),(;k=400)),                 # real to complex scratch
                 ((56,-12,52,0,68),(;k=400)),                 # new magnetic sector
                 ((56,-13,52,0,68),(;k=401)),                 # new target
                 ((20,3,15,-2,25),(;q=0.8)),                  # back to real, smaller arrays
                 ((109//2,-9//2,229//2,-17//2,143),(;q=0.999))) # grow again
        for (args,kw) in cases
            @test qcg(args...;kw...,workspace=work) == qcg(args...;kw...)
            j1,m1,j2,m2,j=args
            @test q3j(j1,j2,j,m1,m2,-m1-m2;kw...,workspace=work) == q3j(j1,j2,j,m1,m2,-m1-m2;kw...)
            # Exercise cache invalidation even if this label's scalar sum is easy.
            doubled=Int.(2 .* args)
            q=get(kw,:q,nothing); k=get(kw,:k,nothing)
            @test isequal(QR._cg_entry(q,k,doubled...;workspace=work),QR._cg_entry(q,k,doubled...))
        end
        # Clearing global tables must not leave stale coefficients in caller scratch.
        args,kw=first(cases)
        QR.clear_cg_caches!()
        @test qcg(args...;kw...,workspace=work) == qcg(args...;kw...)
        labels=fill(args,40)
        @test qcg(labels;kw...,threads=2) == [qcg(a...;kw...) for a in labels]
    end

    @testset "conventions and refusals" begin
        @test qcg(1, 1, 1, 0, 1, 0) == 0                   # m ≠ m₁ + m₂
        @test qcg(1, 2, 1, 0, 1) == 0                      # |m₁| > j₁
        @test q3j(1, 1, 3, 0, 0, 0; q = 0.8) == 0          # triangle
        @test qcg(1//2, 1//2, 1, 0, 1//2; q = -0.8) isa Complex      # negative q: the complex branch, for every label
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
