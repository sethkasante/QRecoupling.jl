using Test
using QRacahSymbols
using Nemo

@testset "QRacahSymbols.jl Master Test Suite" begin
    
    @testset "1. Quantum Dimensions & R-Matrices" begin
        k = 2 # At k=2 (Ising), [2(1/2)+1]_q = [2]_q = sqrt(2)
        
        dim_num = qdim(1/2, k; mode=:numeric)
        @test dim_num ≈ sqrt(2.0) atol=1e-14
        
        # Exact Dimension Check
        dim_ext = qdim(1/2, k; mode=:exact)
        @test isempty(dim_ext.pref_rad.exps) # Should have no radical remainder
        
        # R-Matrix Unitarity Check (|R| = 1)
        r_val = rmatrix(1, 1, 1, 10; mode=:numeric)
        @test abs(r_val) ≈ 1.0 atol=1e-14
    end

    @testset "2. Admissibility & Topological Zeros" begin
        # Triangle inequality violation: |1 - 1| <= 3 <= 1 + 1 is FALSE
        res_cyclo = q6j(1, 1, 3, 1, 1, 1; mode=:cyclo)
        @test res_cyclo.m_min.sign == 0
        @test q6j(1, 1, 3, 1, 1, 1, 10; mode=:numeric) == 0.0
        
        # Level-k cutoff violation: j1+j2+j3 > k (1+1+1 = 3 > 2)
        @test q6j(1, 1, 1, 1, 1, 1, 2; mode=:numeric) == 0.0
        
        # Ensure exact engine correctly catches the zero
        ext_zero = q6j(1, 1, 1, 1, 1, 1, 2; mode=:exact)
        @test iszero(ext_zero.sum_part)
    end

    @testset "3. Engine Consistency (Cyclo vs Numeric)" begin
        j, k = 1, 10
        
        # Compute using the direct Log-Sum-Exp engine
        val_direct = q6j(j, j, j, j, j, j, k; mode=:numeric)
        
        # Compute by building the CycloResult first, then mapping to Floats
        res_cyclo = q6j(j, j, j, j, j, j; mode=:cyclo)
        val_mapped = cyclo_to_numeric(res_cyclo, Float64; k=k)
        
        # They must match perfectly
        @test isapprox(val_direct, val_mapped, atol=1e-12)
        
        # Check Unit Circle projection vs Analytic continuation at theta = pi/4
        theta_val = π/4
        q_val = exp(im * theta_val)
        
        val_unit = cyclo_to_numeric(res_cyclo, ComplexF64; theta=theta_val)
        val_analytic = cyclo_to_numeric(res_cyclo, ComplexF64; q=q_val)
        
        @test isapprox(val_unit, val_analytic, atol=1e-12)
    end

    @testset "4. Classical Limits (q -> 1)" begin
        j1, j2, j3, j4, j5, j6 = 2, 2, 2, 2, 2, 2
        
        # Fast Floats
        val_class_f64 = q6j(j1, j2, j3, j4, j5, j6; mode=:classical)
        
        # Zero-Allocation Exact Rational
        res_class_exact = q6j(j1, j2, j3, j4, j5, j6; mode=:classical_exact)
        val_class_exact_f64 = Float64(res_class_exact)
        
        @test isapprox(val_class_f64, val_class_exact_f64, atol=1e-12)
    end

    @testset "5. Tetrahedral Symmetries" begin
        # Swapping any two columns of a 6j symbol leaves it invariant
        j1, j2, j3, j4, j5, j6 = 1, 1/2, 1, 1/2, 1, 1/2
        k = 5
        
        val_base = q6j(j1, j2, j3, j4, j5, j6, k; mode=:numeric)
        val_swap1 = q6j(j2, j1, j3, j5, j4, j6, k; mode=:numeric) # Swap col 1 & 2
        val_swap2 = q6j(j1, j3, j2, j4, j6, j5, k; mode=:numeric) # Swap col 2 & 3
        
        @test val_base ≈ val_swap1 atol=1e-14
        @test val_base ≈ val_swap2 atol=1e-14
    end

    @testset "6. Exact Orthogonality (The Bubble Move)" begin
        j1, j2, j3, j4, j5 = 1, 1, 1, 1, 1
        k = 6
        h = k + 2
        K, z = cyclotomic_field(2 * h, "ζ")
        
        # We will compute: Sum_x [2x+1]_q * {j1 j2 j3; j4 j5 x}^2
        sum_val = QRacahSymbols.ExactResult(k, QRacahSymbols.EMPTY_MONOMIAL, K(0))
        
        x_min = max(abs(j1 - j5), abs(j2 - j4))
        x_max = min(j1 + j5, j2 + j4, k - max(j1+j5, j2+j4))
        
        for x in x_min:1:x_max
            res_A = q6j(j1, j2, j3, j4, j5, x, k; mode=:exact)
            res_B = q6j(j1, j2, j3, j4, j5, x, k; mode=:exact)
            
            # { } * { } squares the remainder, making it EMPTY_MONOMIAL!
            squared_symbol = res_A * res_B 
            dim_x = qdim(x, k; mode=:exact)
            
            term = squared_symbol * dim_x
            sum_val = sum_val + term
        end
        
        # The exact right-hand side of Orthogonality: 1 / [2j3+1]_q
        dim_j3 = qdim(j3, k; mode=:exact)
        expected_rhs = inv(dim_j3.sum_part)
        
        @test sum_val.sum_part == expected_rhs
    end

    @testset "7. Biedenharn-Elliott Pentagon Identity" begin
        # Sum_x (-1)^{...} [2x+1]_q {}{}{} = {}{}
        j1, j2, j3 = 1, 1, 1
        l1, l2, l3 = 1, 1, 1
        j23, j12 = 1, 1
        k = 8
        
        lhs_sum = 0.0
        x_min = 0
        x_max = k
        
        for x in x_min:x_max
            # Check bounds implicitly via the numeric engine (returns 0 if invalid)
            w1 = q6j(j1, j2, j12, l1, l2, x, k; mode=:numeric)
            w2 = q6j(j1, x, l2, l3, j3, j23, k; mode=:numeric)
            w3 = q6j(l1, j2, x, j3, l3, l2, k; mode=:numeric)
            
            if w1 != 0.0 && w2 != 0.0 && w3 != 0.0
                dim_x = qdim(x, k; mode=:numeric)
                phase_val = iseven(round(Int, j12 + j23 + x + l2)) ? 1.0 : -1.0
                lhs_sum += phase_val * dim_x * w1 * w2 * w3
            end
        end
        
        rhs_w1 = q6j(j12, j3, j23, l3, l1, l2, k; mode=:numeric)
        rhs_w2 = q6j(j1, j2, j12, j3, j23, l1, k; mode=:numeric)
        rhs_val = rhs_w1 * rhs_w2
        
        @test lhs_sum ≈ rhs_val atol=1e-12
    end

    @testset "8. Cache Management" begin
        # Populate caches
        q6j(1, 1, 1, 1, 1, 1, 10; mode=:numeric)
        q6j(1, 1, 1, 1, 1, 1, 10; mode=:exact)
        
        # Clear them
        clear_caches!()
        
        # Verify sizes (assuming they are empty)
        @test length(QRacahSymbols.LOGQFACT_CACHE) == 0
        @test length(QRacahSymbols.EXACT_PHI_CACHE) == 0
    end
end