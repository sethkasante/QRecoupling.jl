

# ---------------------------------------------------------------------
# TQFT Kernels: Quantum Integers, R-Matrices, F-Symbols, and G-Symbols
# Leveraging DCR (Deferred Cyclotomic Representation)
# ---------------------------------------------------------------------


const Spin = Real 

# --- Helper: Algebraic Prefactor Fusion ---

"""
    fuse_prefactors(res::DCR, extra_mono::CyclotomicMonomial)
Fuses a dimensional monomial into the DCR.
Algebraically splits the product of the existing radical and the extra monomial 
into (New Root) * sqrt(New Radical).
"""
function fuse_prefactors(res::DCR, extra_mono::CyclotomicMonomial)
    res.base.sign == 0 && return res # Topological zero
    
    # We need to compute: res.radical * extra_mono
    # and split it into root_bonus * sqrt(new_rad)
    m_d = max(res.max_d, extra_mono.max_d)
    buf = CycloBuffer(m_d)
    
    # Use your high-performance mul! to merge radical and extra dimensions
    mul!(buf, res.radical, extra_mono)
    
    # Extract algebraic squares
    root_bonus, final_rad = snapshot_square_root(buf)
    
    # Update root: New Root = Old Root * Root Bonus
    final_root = res.root * root_bonus
    
    return DCR(final_root, final_rad, res.base, res.ratios, res.z_range, m_d)
end

# ==============================================================================
# 1. QUANTUM INTEGERS & DIMENSIONS
# ==============================================================================

"""
    qint_mono(n::Int)
Returns the CyclotomicMonomial for the quantum integer [n]_q.
"""
function qint_mono(n::Int)
    n == 0 && return CyclotomicMonomial(0, 0, Pair{Int,Int}[], 0)
    n == 1 && return CyclotomicMonomial(1, 0, Pair{Int,Int}[], 1)
    
    buf = CycloBuffer(n)
    add_qint!(buf, n, 1) # [n]_q = product of Phi_d for d|n, d>1
    return snapshot(buf)
end

@inline qdim_mono(j::Spin) = qint_mono(round(Int, 2j + 1))

# ==============================================================================
# 2. BRAIDING (R-Matrix)
# ==============================================================================

"""
    rmatrix_mono(j1, j2, j3)
Algebraic R-matrix phase for braiding j1, j2 into j3.
"""
function rmatrix_mono(j1::Spin, j2::Spin, j3::Spin)
    # Phase = (-1)^{j1+j2-j3} * q^{ (j3(j3+1) - j1(j1+1) - j2(j2+1))/2 }
    # Tracking in q^{1/2} units
    p = round(Int, (j3*(j3+1) - j1*(j1+1) - j2*(j2+1)) * 2) ÷ 2
    s = iseven(round(Int, j1 + j2 - j3)) ? 1 : -1
    return CyclotomicMonomial(s, p, Pair{Int,Int}[], 0)
end

# ==============================================================================
# 3. F-SYMBOLS (Unitary Crossing)
# ==============================================================================

"""
    fsymbol_dcr(j1, j2, j3, j4, j5, j6)
Returns the algebraically fused F-symbol: sqrt([d3][d6]) * {6j}.
"""
function fsymbol_dcr(j1, j2, j3, j4, j5, j6)
    dcr = q6j_dcr(j1, j2, j3, j4, j5, j6)
    dcr.base.sign == 0 && return dcr
    
    # Dimensions to fuse: [2j3+1] * [2j6+1]
    d3 = qdim_mono(j3)
    d6 = qdim_mono(j6)
    dims = d3 * d6
    
    # Fuse into prefactors
    res = fuse_prefactors(dcr, dims)
    
    # Apply crossing phase (-1)^{j1+j2+j4+j5}
    s_cross = iseven(round(Int, j1 + j2 + j4 + j5)) ? 1 : -1
    new_base = CyclotomicMonomial(res.base.sign * s_cross, 
                                  res.base.q_pow, 
                                  res.base.phi_exps, 
                                  res.base.max_d)
    
    return DCR(res.root, res.radical, new_base, res.ratios, res.z_range, res.max_d)
end

# ==============================================================================
# 4. G-SYMBOLS (Tetrahedral Symmetry)
# ==============================================================================

"""
    gsymbol_dcr(j1, j2, j3, j4, j5, j6)
Returns the fully symmetric G-symbol: sqrt(prod [di]) * {6j}.
"""
function gsymbol_dcr(j1, j2, j3, j4, j5, j6)
    dcr = q6j_dcr(j1, j2, j3, j4, j5, j6)
    dcr.base.sign == 0 && return dcr
    
    # Build product of all 6 dimensions
    m_max = round(Int, 2 * max(j1, j2, j3, j4, j5, j6) + 1)
    buf = CycloBuffer(m_max)
    for j in (j1, j2, j3, j4, j5, j6)
        add_qint!(buf, round(Int, 2j + 1), 1)
    end
    dims = snapshot(buf)
    
    return fuse_prefactors(dcr, dims)
end








# # ---------------------------------------------------------------------
# # File: tqft_kernels.jl
# # Topological Kernels: Q-Integers, R-Matrices, F-Symbols, and G-Symbols
# # ---------------------------------------------------------------------

# # --- Algebraic Fusion Engine ---


# const Spin = Real

# """
#     fuse_prefactors(res::DCR, extra::CyclotomicMonomial)
# Fuses a dimensional monomial into a DCR prefactor. 
# Extracts perfect squares algebraically to keep the radical minimal.
# """
# function fuse_prefactors(res::DCR, extra::CyclotomicMonomial)
#     # If the DCR is a topological zero, stay zero
#     res.base.sign == 0 && return res
    
#     # We use a temporary buffer to merge radical and extra dimensions
#     # max_d check ensures we have room for the new dimensions
#     m_d = max(res.max_d, extra.max_d)
#     buf = CycloBuffer(m_d) 
    
#     # Load Radical
#     buf.sign = res.radical.sign * extra.sign
#     buf.q_pow = res.radical.q_pow + extra.q_pow
#     buf.max_d = res.radical.max_d
#     for (d, e) in res.radical.phi_exps
#         buf.exps[d] += e
#     end
    
#     # Merge Extra Dimensions
#     for (d, e) in extra.phi_exps
#         buf.exps[d] += e
#     end
#     buf.max_d = max(buf.max_d, extra.max_d)
    
#     # Re-extract squares
#     bonus_root, final_rad = snapshot_square_root(buf)
    
#     # Multiply the existing root by the new bonus root
#     # (Note: Need a simple Monomial-Monomial multiplier)
#     final_root = multiply_monomials(res.root, bonus_root)
    
#     return DCR(final_root, final_rad, res.base, res.ratios, res.z_range, m_d)
# end

# @inline function multiply_monomials(m1::CyclotomicMonomial, m2::CyclotomicMonomial)
#     m1.sign == 0 && return m1
#     m2.sign == 0 && return m2
    
#     new_q = m1.q_pow + m2.q_pow
#     new_sign = m1.sign * m2.sign
    
#     # Merging sparse exps (Linear time merge)
#     # In practice, D is small (< 1000), so we reuse a buffer
#     buf_exps = zeros(Int, max(m1.max_d, m2.max_d))
#     for (d, e) in m1.phi_exps; buf_exps[d] += e; end
#     for (d, e) in m2.phi_exps; buf_exps[d] += e; end
    
#     nnz = count(!iszero, buf_exps)
#     sparse = Vector{Pair{Int,Int}}(undef, nnz)
#     idx = 1
#     for d in 1:length(buf_exps)
#         if buf_exps[d] != 0
#             sparse[idx] = d => buf_exps[d]
#             idx += 1
#         end
#     end
    
#     return CyclotomicMonomial(new_sign, new_q, sparse, length(buf_exps))
# end

# # ==============================================================================
# # 1. QUANTUM INTEGERS & DIMENSIONS
# # ==============================================================================

# """
#     qint_mono(n::Int) -> CyclotomicMonomial
# The algebraic [n]_q. Represented as q-power and cyclotomic factors.
# """
# function qint_mono(n::Int)
#     n == 0 && return CyclotomicMonomial(0, 0, [], 0)
#     n == 1 && return CyclotomicMonomial(1, 0, [], 1)
    
#     buf = CycloBuffer(n)
#     # [n]_q = (q^{n/2} - q^{-n/2}) / (q^{1/2} - q^{-1/2})
#     # Algebraically, this is the product of Phi_d for d|n, d>1.
#     add_qint!(buf, n, 1)
#     return snapshot(buf)
# end

# function qdim_mono(j::Spin)
#     return qint_mono(round(Int, 2j + 1))
# end

# # ==============================================================================
# # 2. BRAIDING & PHASES (R-Matrix)
# # ==============================================================================

# """
#     rmatrix_mono(j1, j2, j3)
# The braiding phase for the fusion of j1 and j2 into j3.
# """
# function rmatrix_mono(j1::Spin, j2::Spin, j3::Spin)
#     # Phase = (-1)^{j1+j2-j3} * q^{ (j3(j3+1) - j1(j1+1) - j2(j2+1))/2 }
#     # Since j(j+1) = 2j(2j+2)/4, we track powers in q^{1/2} units.
#     p = round(Int, (j3*(j3+1) - j1*(j1+1) - j2*(j2+1)) * 2) ÷ 2
#     s = iseven(round(Int, j1 + j2 - j3)) ? 1 : -1
#     return CyclotomicMonomial(s, p, [], 0)
# end

# # ==============================================================================
# # 3. F-SYMBOLS & G-SYMBOLS (Crossing & Symmetry)
# # ==============================================================================

# """
#     fsymbol_dcr(j1, j2, j3, j4, j5, j6)
# Algebraically fused F-symbol (Unitary normalized 6j).
# """
# function fsymbol_dcr(j1, j2, j3, j4, j5, j6)
#     dcr = q6j_dcr(j1, j2, j3, j4, j5, j6)
#     dcr.base.sign == 0 && return dcr

#     # F = (-1)^{j1+j2+j4+j5} * sqrt([2j3+1][2j6+1]) * {6j}
#     buf = CycloBuffer(round(Int, max(2j3+1, 2j6+1)))
#     add_qint!(buf, round(Int, 2j3+1), 1)
#     add_qint!(buf, round(Int, 2j6+1), 1)
#     dims = snapshot(buf)

#     # Fuse and apply crossing phase
#     res = fuse_prefactors(dcr, dims)
#     s_cross = iseven(round(Int, j1 + j2 + j4 + j5)) ? 1 : -1
    
#     # Direct update of base sign to avoid copy
#     new_base = CyclotomicMonomial(res.base.sign * s_cross, res.base.q_pow, res.base.phi_exps, res.base.max_d)
    
#     return DCR(res.root, res.radical, new_base, res.ratios, res.z_range, res.max_d)
# end

# """
#     gsymbol_dcr(j1, j2, j3, j4, j5, j6)
# Fully tetrahedral symmetric G-symbol.
# """
# function gsymbol_dcr(j1, j2, j3, j4, j5, j6)
#     dcr = q6j_dcr(j1, j2, j3, j4, j5, j6)
#     dcr.base.sign == 0 && return dcr
    
#     # G = sqrt([2j1+1][2j2+1][2j3+1][2j4+1][2j5+1][2j6+1]) * {6j}
#     m_s = round(Int, 2*max(j1, j2, j3, j4, j5, j6) + 1)
#     buf = CycloBuffer(m_s)
#     for j in (j1, j2, j3, j4, j5, j6)
#         add_qint!(buf, round(Int, 2j + 1), 1)
#     end
#     dims = snapshot(buf)
    
#     return fuse_prefactors(dcr, dims)
# end