
#add tests 

using Test
using QRacahSymbols
# using Nemo


@testset "QRacahSymbols.jl Core Tests" begin
    
    @testset "Quantum Integers" begin
       
        k = 2
        # [2] at k=2 is sin(2π/4)/sin(π/4) = 1/0.707... = sqrt(2)
        q2_num = qint(2, k; mode=:numeric)
        @test q2_num ≈ sqrt(2.0)

        # exact algebraic [2]
        exact_q2 = qint(2, k; mode=:exact)
        @test evaluate_exact(exact_q2, Float64) ≈ q2_num atol=1e-14

        # Dimension of spin-1/2 should be [2(1/2)+1]_q = [2]_q
        dim2_num = qdim(1/2, k; mode=:numeric)
        dim2_ext = qdim(1/2, k; mode=:exact)
        @test dim2_num ≈ sqrt(2.0) atol=1e-14
        @test evaluate_exact(dim2_ext, Float64) ≈ dim2_num atol=1e-14
        
        # Symbolic Mode
        symb_q2 = qint(2; mode=:generic)
        @test symb_q2.exps[2] == 1 # [2]_q = Φ_2(q)
        
    end

    @testset "3j and 6j Consistency" begin
        j = 1
        k = 10
        
        # Compute the same 6j symbol across all 3 modes
        val_num = q6j(j, j, j, j, j, j, k; mode=:numeric)
        res_ext = q6j(j, j, j, j, j, j, k; mode=:exact)
        res_gen = q6j(j, j, j, j, j, j; mode=:generic)
        
        # Evaluate Exact and Generic to Float64
        val_ext_float = evaluate_exact(res_ext, Float64)
        val_gen_float = evaluate_generic(res_gen, k,Float64)
        
        # All three must match perfectly
        @test isapprox(val_num, val_ext_float, atol=1e-12)
        @test isapprox(val_num, val_gen_float, atol=1e-12)
    end

    @testset "Admissibility" begin
        # Triangle inequality failure
        # Triangle inequality violation: |1 - 1| <= 3 <= 1 + 1 is FALSE
        @test q6j(1, 1, 3, 1, 1, 1; mode=:generic).pref_sq.sign == 0
        @test q6j(1, 1, 3, 1, 1, 1, 10; mode=:numeric) == 0.0
        
        # Level-k cutoff violation: j1+j2+j3 > k (1+1+1 = 3 > 2)
        @test q6j(1, 1, 1, 1, 1, 1, 2; mode=:numeric) == 0.0
        @test q6j(1, 1, 1, 1, 1, 1, 2; mode=:exact).pref_sq == 0
    end

    @testset "6j Tetrahedral Symmetries" begin
        # Swapping any two columns of a 6j symbol leaves it invariant
        j1, j2, j3, j4, j5, j6 = 1, 1/2, 1, 1/2, 1, 1/2
        k = 5
        
        val1 = q6j(j1, j2, j3, j4, j5, j6, k; mode=:numeric)
        val2 = q6j(j2, j1, j3, j5, j4, j6, k; mode=:numeric) # Swap col 1 & 2
        val3 = q6j(j1, j3, j2, j4, j6, j5, k; mode=:numeric) # Swap col 2 & 3
        
        @test val1 ≈ val2 
        @test val1 ≈ val3 
    end

    @testset "5. TQFT Category Data" begin
        j = 1
        k = 10
        
        # F-symbol normalization check
        # F = (-1)^{...} * sqrt([2j3+1][2j6+1]) * {6j}
        f_val = fsymbol(j, j, j, j, j, j, k; mode=:numeric)
        q6j_val = q6j(j, j, j, j, j, j, k; mode=:numeric)
        dim_j = qdim(j, k; mode=:numeric)
        
        expected_f = ((-1)^(4j)) * sqrt(dim_j * dim_j) * q6j_val
        @test f_val ≈ expected_f atol=1e-12
        
        # R-Matrix Phase
        r_val = rmatrix(1, 1, 1, k; mode=:numeric)
        @test abs(r_val) ≈ 1.0 atol=1e-14 # R-matrix is unitary
    end

    @testset "High-Spin Stable Symbolic Evaluation" begin
        # Test that evaluate_generic doesn't crash or drift on larger spins
        jj = 20
        kk = 100
        
        qs6 = q6j(jj, jj, jj, jj, jj, jj; mode=:generic)
        val_symb = evaluate_generic(qs6, kk, Float64; prec=256)
        val_num  = q6j(jj, jj, jj, jj, jj, jj, kk; mode=:numeric)
        
        @test isapprox(val_symb, val_num, atol=1e-10)
    end

    @testset "Cache Management" begin
        # Ensure the cache clearing function runs without error
        @test clear_caches!() === nothing
        # Ensure doing it twice doesn't break anything
        @test clear_caches!() === nothing
    end
end

# Set global precision high enough to rigorously test BigFloat vs Exact Field projections
# setprecision(BigFloat, 256)

# @testset "QRacahSymbols.jl Master Validation Suite" begin

#     @testset "1. Admissibility & Gatekeepers" begin
#         # Classical Triangle
#         @test QRacahSymbols.δ(1, 1, 1) == true
#         @test QRacahSymbols.δ(1, 1, 3) == false  # Fails classical triangle
        
#         # Quantum Triangle (k limit)
#         @test QRacahSymbols.qδ(1, 1, 1, 3) == true
#         @test QRacahSymbols.qδ(1, 1, 1, 2) == false # Sum exceeds k
        
#         # 6j Symbol Admissibility
#         # Valid configuration
#         @test QRacahSymbols.qδtet(1, 1, 1, 1, 1, 1, 4) == true
        
#         # Invalid configuration (returns 0 structurally)
#         res_numeric = q6j(1, 1, 1, 1, 1, 1, 2; mode=:numeric)
#         res_exact = q6j(1, 1, 1, 1, 1, 1, 2; mode=:exact)
#         res_generic = q6j(1, 1, 1, 1, 1, 1; mode=:generic) # k-independent
        
#         @test res_numeric == 0.0
#         @test res_exact.pref_sq == 0 && res_exact.sum_cf == 0
#         @test res_generic.pref_sq.sign == 0 
#     end

#     @testset "2. Tetrahedral Symmetries (S_4)" begin
#         # The 6j symbol is invariant under 24 permutations of the tetrahedron
#         j1, j2, j3, j4, j5, j6 = 1, 1.5, 2, 2.5, 1, 1.5
#         k = 10
        
#         base_val = q6j(j1, j2, j3, j4, j5, j6, k; mode=:numeric, T=Float64)
        
#         # Swap columns 1 and 2
#         @test q6j(j2, j1, j3, j5, j4, j6, k; mode=:numeric, T=Float64) ≈ base_val
#         # Swap columns 2 and 3
#         @test q6j(j1, j3, j2, j4, j6, j5, k; mode=:numeric, T=Float64) ≈ base_val
#         # Flip columns 1 and 2
#         @test q6j(j4, j5, j3, j1, j2, j6, k; mode=:numeric, T=Float64) ≈ base_val
#     end

#     @testset "3. The Golden Convergence (Exact vs Numeric)" begin
#         # We test that the exact algebraic number field projects perfectly 
#         # to the numeric log-sum-exp floating point evaluator.
        
#         j_val = 1.0
#         k_val = 4
        
#         # 1. 6j Convergence
#         res_num_6j = q6j(j_val, j_val, j_val, j_val, j_val, j_val, k_val; mode=:numeric, T=BigFloat)
#         res_exact_6j = q6j(j_val, j_val, j_val, j_val, j_val, j_val, k_val; mode=:exact)
#         res_proj_6j = evaluate_exact(res_exact_6j)
        
#         # We expect precision matching up to BigFloat eps (~1e-70 for 256-bit)
#         @test isapprox(res_num_6j, res_proj_6j, atol=1e-60)
        
#         # 2. 3j Convergence
#         m_val = 0.0
#         res_num_3j = q3j(j_val, j_val, j_val, m_val, m_val, m_val, k_val; mode=:numeric, T=BigFloat)
#         res_exact_3j = q3j(j_val, j_val, j_val, m_val, m_val, m_val, k_val; mode=:exact)
#         res_proj_3j = evaluate_exact(res_exact_3j)
        
#         @test isapprox(res_num_3j, res_proj_3j, atol=1e-60)
#     end

#     @testset "4. Symbolic Algebra Engine" begin
#         # Verify the in-place array manipulation math is structurally sound
#         a = CycloMonomial(1, 2, [0, 1, 2]) # 1 * z^2 * Φ_1^1 * Φ_2^2
#         b = CycloMonomial(-1, 3, [0, 0, 1, 4]) # -1 * z^3 * Φ_2^1 * Φ_3^4
        
#         # Multiplication
#         c = a * b
#         @test c.sign == -1
#         @test c.z_pow == 5
#         @test c.exps == [0, 1, 3, 4]
        
#         # Division
#         d = c / a
#         @test d.sign == -1
#         @test d.z_pow == 3
#         @test d.exps == [0, 0, 1, 4]
#     end

#     @testset "5. TQFT Category Data Validation" begin
#         j = 1.5
#         k = 5
        
#         # 1. Quantum Dimensions
#         num_dim = qdim(j, k; mode=:numeric, T=BigFloat)
#         exact_dim = qdim(j, k; mode=:exact)
        
#         # Project Exact dimension to BigFloat for comparison
#         target_z = cispi(big"1.0" / (k + 2))
#         proj_dim = real(QRacahSymbols.horner_eval(exact_dim, target_z))
#         @test isapprox(num_dim, proj_dim, atol=1e-60)
        
#         # 2. R-Matrix (Braiding) Phase consistency
#         num_rmat = rmatrix(1, 1, 1, k; mode=:numeric, T=BigFloat)
#         exact_rmat = rmatrix(1, 1, 1, k; mode=:exact)
#         proj_rmat = QRacahSymbols.horner_eval(exact_rmat, target_z)
#         @test isapprox(real(num_rmat), real(proj_rmat), atol=1e-60)
#         @test isapprox(imag(num_rmat), imag(proj_rmat), atol=1e-60)
#     end

#     @testset "6. Classical Limit (q -> 1)" begin
#         # Test that the Ponzano-Regge classical evaluator runs without errors
#         # and correctly parses a generic result into a Float64 limit.
#         class_6j = q6j(1, 1, 1, 1, 1, 1; mode=:classical)
#         @test class_6j isa Float64
        
#         class_3j = q3j(1, 1, 1, 0, 0, 0; mode=:classical)
#         @test class_3j isa Float64
#     end

# end