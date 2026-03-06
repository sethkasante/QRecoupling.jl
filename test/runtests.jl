

using Test
using QRacahSymbols
using Nemo

# Set global precision high enough to rigorously test BigFloat vs Exact Field projections
setprecision(BigFloat, 256)

@testset "QRacahSymbols.jl Master Validation Suite" begin

    @testset "1. Admissibility & Gatekeepers" begin
        # Classical Triangle
        @test QRacahSymbols.δ(1, 1, 1) == true
        @test QRacahSymbols.δ(1, 1, 3) == false  # Fails classical triangle
        
        # Quantum Triangle (k limit)
        @test QRacahSymbols.qδ(1, 1, 1, 3) == true
        @test QRacahSymbols.qδ(1, 1, 1, 2) == false # Sum exceeds k
        
        # 6j Symbol Admissibility
        # Valid configuration
        @test QRacahSymbols.qδtet(1, 1, 1, 1, 1, 1, 4) == true
        
        # Invalid configuration (returns 0 structurally)
        res_numeric = q6j(1, 1, 1, 1, 1, 1, 2; mode=:numeric)
        res_exact = q6j(1, 1, 1, 1, 1, 1, 2; mode=:exact)
        res_generic = q6j(1, 1, 1, 1, 1, 1; mode=:generic) # k-independent
        
        @test res_numeric == 0.0
        @test res_exact.pref_sq == 0 && res_exact.sum_cf == 0
        @test res_generic.pref_sq.sign == 0 
    end

    @testset "2. Tetrahedral Symmetries (S_4)" begin
        # The 6j symbol is invariant under 24 permutations of the tetrahedron
        j1, j2, j3, j4, j5, j6 = 1, 1.5, 2, 2.5, 1, 1.5
        k = 10
        
        base_val = q6j(j1, j2, j3, j4, j5, j6, k; mode=:numeric, T=Float64)
        
        # Swap columns 1 and 2
        @test q6j(j2, j1, j3, j5, j4, j6, k; mode=:numeric, T=Float64) ≈ base_val
        # Swap columns 2 and 3
        @test q6j(j1, j3, j2, j4, j6, j5, k; mode=:numeric, T=Float64) ≈ base_val
        # Flip columns 1 and 2
        @test q6j(j4, j5, j3, j1, j2, j6, k; mode=:numeric, T=Float64) ≈ base_val
    end

    @testset "3. The Golden Convergence (Exact vs Numeric)" begin
        # We test that the exact algebraic number field projects perfectly 
        # to the numeric log-sum-exp floating point evaluator.
        
        j_val = 1.0
        k_val = 4
        
        # 1. 6j Convergence
        res_num_6j = q6j(j_val, j_val, j_val, j_val, j_val, j_val, k_val; mode=:numeric, T=BigFloat)
        res_exact_6j = q6j(j_val, j_val, j_val, j_val, j_val, j_val, k_val; mode=:exact)
        res_proj_6j = evaluate_exact(res_exact_6j)
        
        # We expect precision matching up to BigFloat eps (~1e-70 for 256-bit)
        @test isapprox(res_num_6j, res_proj_6j, atol=1e-60)
        
        # 2. 3j Convergence
        m_val = 0.0
        res_num_3j = q3j(j_val, j_val, j_val, m_val, m_val, m_val, k_val; mode=:numeric, T=BigFloat)
        res_exact_3j = q3j(j_val, j_val, j_val, m_val, m_val, m_val, k_val; mode=:exact)
        res_proj_3j = evaluate_exact(res_exact_3j)
        
        @test isapprox(res_num_3j, res_proj_3j, atol=1e-60)
    end

    @testset "4. Symbolic Algebra Engine" begin
        # Verify the in-place array manipulation math is structurally sound
        a = CycloMonomial(1, 2, [0, 1, 2]) # 1 * z^2 * Φ_1^1 * Φ_2^2
        b = CycloMonomial(-1, 3, [0, 0, 1, 4]) # -1 * z^3 * Φ_2^1 * Φ_3^4
        
        # Multiplication
        c = a * b
        @test c.sign == -1
        @test c.z_pow == 5
        @test c.exps == [0, 1, 3, 4]
        
        # Division
        d = c / a
        @test d.sign == -1
        @test d.z_pow == 3
        @test d.exps == [0, 0, 1, 4]
    end

    @testset "5. TQFT Category Data Validation" begin
        j = 1.5
        k = 5
        
        # 1. Quantum Dimensions
        num_dim = qdim(j, k; mode=:numeric, T=BigFloat)
        exact_dim = qdim(j, k; mode=:exact)
        
        # Project Exact dimension to BigFloat for comparison
        target_z = cispi(big"1.0" / (k + 2))
        proj_dim = real(QRacahSymbols.horner_eval(exact_dim, target_z))
        @test isapprox(num_dim, proj_dim, atol=1e-60)
        
        # 2. R-Matrix (Braiding) Phase consistency
        num_rmat = rmatrix(1, 1, 1, k; mode=:numeric, T=BigFloat)
        exact_rmat = rmatrix(1, 1, 1, k; mode=:exact)
        proj_rmat = QRacahSymbols.horner_eval(exact_rmat, target_z)
        @test isapprox(real(num_rmat), real(proj_rmat), atol=1e-60)
        @test isapprox(imag(num_rmat), imag(proj_rmat), atol=1e-60)
    end

    @testset "6. Classical Limit (q -> 1)" begin
        # Test that the Ponzano-Regge classical evaluator runs without errors
        # and correctly parses a generic result into a Float64 limit.
        class_6j = q6j(1, 1, 1, 1, 1, 1; mode=:classical)
        @test class_6j isa Float64
        
        class_3j = q3j(1, 1, 1, 0, 0, 0; mode=:classical)
        @test class_3j isa Float64
    end

end