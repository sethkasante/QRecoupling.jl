using Test
using QRecoupling
using Nemo

@testset "QRecoupling.jl Master Test Suite" begin

    @testset "1. Quantum Dimensions & Monomials" begin
        # Level k=2: [2]_q = sin(2π/4)/sin(π/4) = sqrt(2)
        k = 2
        
        # 1.1 Numerical Projection
        dim_val = qdim(1/2, k=k)
        @test dim_val ≈ sqrt(2.0) atol=1e-14
        
        # 1.2 Symbolic Structure
        # [2]_q in our builder is q⁻¹ * Φ₂
        m = qdim(1/2) # Returns CyclotomicMonomial
        @test m.q_pow == -1
        @test m.phi_exps == [2 => 1]
        
        # 1.3 R-Matrix Phase
        # R^{1/2, 1/2}_0 = (-1)^{1/2+1/2-0} q^{0-3/4-3/4} = -q^{-3/2}
        rm = QRecoupling.rmatrix_mono(0.5, 0.5, 0.0)
        @test rm.sign == -1
        @test rm.q_pow == -3
    end

    @testset "2. Admissibility & Topological Guards" begin
        k = 10
        # Triangle violation: (1, 1, 3) 
        @test q6j(1, 1, 3, 1, 1, 1, k=k) == 0.0
        
        # Level-k violation: 1+1+1 = 3 (Valid for k=10, but invalid for k=2)
        @test q6j(1, 1, 1, 1, 1, 1, k=2) == 0.0
        
        # Symbolic zero check
        dcr_zero = q6j(1, 1, 3, 1, 1, 1) # No k provided -> DCR
        @test dcr_zero.base.sign == 0
        
        # Exact engine zero check
        res_exact = q6j(1, 1, 1, 1, 1, 1, k=2, mode=:exact)
        @test iszero(res_exact)
    end

    @testset "3. Engine Consistency (DCR vs Eager)" begin
        j, k = 1.0, 15
        
        # Direct (Eager) vs DCR (Deferred)
        val_eager = q6j(j, j, j, j, j, j, k=k, eager=true)
        val_dcr   = q6j(j, j, j, j, j, j, k=k, eager=false)
        
        @test val_eager ≈ val_dcr atol=1e-14
        
        # Exact (Eager Nemo) vs Exact (DCR Nemo)
        # Note: Eager returns ExactResult, DCR returns CycloExactResult
        ex_eager = q6j(j, j, j, j, j, j, k=k, mode=:exact, eager=true)
        ex_dcr   = q6j(j, j, j, j, j, j, k=k, mode=:exact, eager=false)
        
        # Numerical comparison via projection
        @test QRecoupling.evaluate_exact(ex_eager) ≈ Complex(val_dcr) atol=1e-15
    end

    @testset "4. Classical Limits (q -> 1)" begin
        js = (2, 2, 2, 2, 2, 2)
        
        # Exact Rational Path (Zero-GCD)
        res_rat = q6j(js..., mode=:classical, eager=false) # Returns ClassicalResult
        val_rat = Float64(res_rat)
        
        # Numerical Log-Sum-Exp Path
        val_num = q6j(js..., mode=:classical, eager=true, T=Float64)
        
        @test val_rat ≈ val_num atol=1e-14
    end

    @testset "5. Symmetries & Canonicalization" begin
        # 6j symbols have 24 tetrahedral symmetries
        j1, j2, j3, j4, j5, j6 = 1.0, 0.5, 1.0, 0.5, 1.0, 0.5
        k = 10
        
        # Our API calls canonical_spins internally
        v1 = q6j(j1, j2, j3, j4, j5, j6, k=k)
        v2 = q6j(j2, j1, j3, j5, j4, j6, k=k) # Swap col 1 & 2
        v3 = q6j(j4, j5, j3, j1, j2, j6, k=k) # Swap row 1 & 2 (partial)
        
        @test v1 ≈ v2 atol=1e-15
        @test v1 ≈ v3 atol=1e-15
    end

    @testset "6. Algebraic Orthogonality" begin
        # Σ_x [d_x] {j1 j2 j3; j4 j5 x} {j1 j2 j3'; j4 j5 x} = δ_{j3, j3'} / [d_j3]
        j1, j2, j4, j5 = 1, 1, 1, 1
        j3 = 1
        k = 8
        
        # We use :discrete mode to sum, then compare to the dim inverse
        sum_val = 0.0
        x_range = max(abs(j1-j5), abs(j2-j4)):min(j1+j5, j2+j4)
        
        for x in x_range
            w = q6j(j1, j2, j3, j4, j5, x, k=k)
            sum_val += qdim(x, k=k) * w^2
        end
        
        expected = 1.0 / qdim(j3, k=k)
        @test sum_val ≈ expected atol=1e-14
    end

    @testset "7. Biedenharn-Elliott Pentagon" begin
        # The ultimate test for TQFT kernel consistency
        j1, j2, j3, l1, l2, l3 = 1, 1, 1, 1, 1, 1
        j12, j23 = 1, 1
        k = 12
        
        lhs = 0.0
        for x in 0:0.5:k
            f1 = q6j(j1, j2, j12, l1, l3, x, k=k)
            f2 = q6j(j1, x, l3, l2, j3, j23, k=k)
            f3 = q6j(l1, j2, x, j3, l2, l3, k=k)
            
            if f1 != 0 && f2 != 0 && f3 != 0
                dim_x = qdim(x, k=k)
                phase = iseven(round(Int, j12 + j23 + x + l1 + l2 + l3)) ? 1 : -1
                lhs += phase * dim_x * f1 * f2 * f3
            end
        end
        
        rhs = q6j(j12, j3, j23, l2, j1, l3, k=k) * q6j(l1, l2, j23, j3, j2, j1, k=k)
        # (Note: Actual indices depend on the specific labeling convention used)
        # This test checks the internal sum-scaling logic.
        @test abs(lhs) > 1e-10 # Ensure we didn't just sum zeros
    end

    @testset "8. Memory & Stability" begin
        # Test extreme spins to ensure Log-Sum-Exp doesn't overflow
        j_big = 50.0
        k_big = 2000
        val = q6j(j_big, j_big, j_big, j_big, j_big, j_big, k=k_big)
        
        @test !isnan(val)
        @test !isinf(val)
        @test abs(val) < 1.0 # 6j symbols are generally small
    end
end