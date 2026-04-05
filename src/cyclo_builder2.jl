# ---------------------------------------------------------------------
# File: builders.jl
# High-level builders for DCR Symbols and TQFT Kernels.
# ---------------------------------------------------------------------

# --- Constants ---
const ZERO_DCR = DCR(ZERO_MONOMIAL, ONE_MONOMIAL, ZERO_MONOMIAL, CyclotomicMonomial[], 0:-1, 0)

# ==============================================================================
# 1. ALGEBRAIC UTILITIES
# ==============================================================================

"""
    fuse_prefactors(res::DCR, extra_mono::CyclotomicMonomial)
Fuses an additional dimensional monomial into a DCR.
Algebraically extracts perfect squares from the product of the existing 
radical and the extra monomial, moving them into the 'root'.
"""
function fuse_prefactors(res::DCR, extra_mono::CyclotomicMonomial)
    res.base.sign == 0 && return res
    
    m_d = max(res.max_d, extra_mono.max_d)
    buf = CycloBuffer(m_d)
    
    # Merge existing radical with the new dimensions
    mul!(buf, res.radical, extra_mono)
    
    # Extract squares: radical * extra = (bonus_root)^2 * final_rad
    root_bonus, final_rad = snapshot_square_root(buf)
    
    # Update prefactors
    final_root = res.root * root_bonus
    
    return DCR(final_root, final_rad, res.base, res.ratios, res.z_range, m_d)
end

# ==============================================================================
# 2. MONOMIAL BUILDERS (Dimensions & Phases)
# ==============================================================================

"""
    qint_mono(n::Int) -> CyclotomicMonomial
Returns the algebraic representation of the quantum integer [n]q.
"""
function qint_mono(n::Int)
    n == 0 && return ZERO_MONOMIAL
    n == 1 && return ONE_MONOMIAL
    buf = CycloBuffer(n)
    add_qint!(buf, n, 1)
    return snapshot(buf)
end

@inline qdim_mono(j::Spin) = qint_mono(round(Int, 2j + 1))

"""
    rmatrix_mono(j1, j2, j3) -> CyclotomicMonomial
Algebraic R-matrix phase for braiding j1, j2 into j3.
Uses the q^{1/2} power convention.
"""
function rmatrix_mono(j1::Spin, j2::Spin, j3::Spin)
    # Phase = (-1)^{j1+j2-j3} * q^{ (j3(j3+1) - j1(j1+1) - j2(j2+1))/2 }
    # Tracking in q^{1/2} units: P = j3(j3+1) - j1(j1+1) - j2(j2+1)
    p = round(Int, (j3*(j3+1) - j1*(j1+1) - j2*(j2+1)) * 2) ÷ 2
    s = iseven(round(Int, j1 + j2 - j3)) ? 1 : -1
    return CyclotomicMonomial(s, p, Pair{Int,Int}[], 0)
end

# ==============================================================================
# 3. DCR ASSEMBLY ENGINE
# ==============================================================================

"""
    build_dcr!(buf, pre_func, base_func, ratio_func, z_min, z_max)
Universal skeletal assembler for q-hypergeometric series.
"""
function build_dcr!(buf, pre_func, base_func, ratio_func, z_min::Int, z_max::Int)
    reset!(buf)
    pre_func(buf)
    m_root, m_rad = snapshot_square_root(buf)
    g_max_d = buf.max_d

    # Topological Zero or Empty Range
    (z_min > z_max || buf.sign == 0) && return DCR(m_root, m_rad, ZERO_MONOMIAL, CyclotomicMonomial[], z_min:z_max, g_max_d)

    # Base Term at z_min
    reset!(buf, iseven(z_min) ? 1 : -1)
    base_func(buf, z_min)
    m_base = snapshot(buf)
    g_max_d = max(g_max_d, buf.max_d)

    # Ratio Sequence
    ratios = Vector{CyclotomicMonomial}(undef, z_max - z_min)
    for (i, z) in enumerate(z_min:z_max-1)
        reset!(buf, -1) # Ratios in these series usually carry a -1 flip
        ratio_func(buf, z)
        ratios[i] = snapshot(buf)
        g_max_d = max(g_max_d, buf.max_d)
    end

    return DCR(m_root, m_rad, m_base, ratios, z_min:z_max, g_max_d)
end

# ==============================================================================
# 4. SYMBOL BUILDERS (3j & 6j)
# ==============================================================================

@inline function qtriangle!(buf, j1, j2, j3)
    add_qfact!(buf, round(Int, j1 + j2 - j3))
    add_qfact!(buf, round(Int, j1 - j2 + j3))
    add_qfact!(buf, round(Int, -j1 + j2 + j3))
    add_qfact!(buf, round(Int, j1 + j2 + j3 + 1), -1)
end

@inline function qtetrahedron!(buf, j1, j2, j3, j4, j5, j6)
    qtriangle!(buf, j1, j2, j3); qtriangle!(buf, j1, j5, j6)
    qtriangle!(buf, j2, j4, j6); qtriangle!(buf, j3, j4, j5)
end

function q3j_dcr(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    # 1. Topological & Selection Rule Admissibility
    (!δ(j1, j2, j3) || m1 + m2 + m3 != 0) && return ZERO_DCR
    
    # Standard summation bounds for Wigner 3j
    α = (round(Int, j3 - j2 + m1), round(Int, j3 - j1 - m2))
    β = (round(Int, j1 + j2 - j3), round(Int, j1 - m1), round(Int, j2 + m2))
    z_min, z_max = max(0, -α[1], -α[2]), min(β[1], β[2], β[3])
    
    # Initialize buffer with sufficient capacity
    buf = CycloBuffer(max(z_max + 2, round(Int, j1 + j2 + j3 + 1)))

    return build_dcr!(buf,
        # Prefactor: Triangle factorials * (j ± m) factorials
        b -> begin
            qtriangle!(b, j1, j2, j3)
            # Add (j ± m)! terms explicitly for clarity and speed
            add_qfact!(b, round(Int, j1 + m1)); add_qfact!(b, round(Int, j1 - m1))
            add_qfact!(b, round(Int, j2 + m2)); add_qfact!(b, round(Int, j2 - m2))
            add_qfact!(b, round(Int, j3 + m3)); add_qfact!(b, round(Int, j3 - m3))
        end,
        # Base Term (at z_min)
        (b, z) -> begin
            add_qfact!(b, z, -1)
            for a in α; add_qfact!(b, a + z, -1); end
            for bv in β; add_qfact!(b, bv - z, -1); end
        end,
        # Ratio R_z = Term(z+1)/Term(z)
        (b, z) -> begin
            for bv in β; add_qint!(b, bv - z); end
            add_qint!(b, z + 1, -1)
            for a in α; add_qint!(b, a + z + 1, -1); end
        end,
        z_min, z_max
    )
end

function q6j_dcr(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    !δtet(j1, j2, j3, j4, j5, j6) && return ZERO_DCR
    α = (round(Int, j1+j2+j3), round(Int, j1+j5+j6), round(Int, j2+j4+j6), round(Int, j3+j4+j5))
    β = (round(Int, j1+j2+j4+j5), round(Int, j1+j3+j4+j6), round(Int, j2+j3+j5+j6))
    z_min, z_max = max(α...), min(β...)
    buf = CycloBuffer(max(z_max + 2, β...))

    return build_dcr!(buf,
        b -> qtetrahedron!(b, j1, j2, j3, j4, j5, j6),
        (b, z) -> (add_qfact!(b, z+1); for a in α add_qfact!(b, z-a, -1) end; for bv in β add_qfact!(b, bv-z, -1) end),
        (b, z) -> (add_qint!(b, z+2); for bv in β add_qint!(b, bv-z) end; for a in α add_qint!(b, z+1-a, -1) end),
        z_min, z_max
    )
end

# ==============================================================================
# 5. TQFT KERNELS (F & G Symbols)
# ==============================================================================

"""
    fsymbol_dcr(j1, j2, j3, j4, j5, j6)
Algebraically fused F-symbol: sqrt([2j3+1][2j6+1]) * {6j}.
Matches the unitary crossing matrix in Spin Networks.
"""
function fsymbol_dcr(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    dcr = q6j_dcr(j1, j2, j3, j4, j5, j6)
    dcr.base.sign == 0 && return dcr
    
    dims = qdim_mono(j3) * qdim_mono(j6)
    res = fuse_prefactors(dcr, dims)
    
    phase = iseven(round(Int, j1 + j2 + j4 + j5)) ? 1 : -1
    # Inline sign update to avoid full DCR copy
    new_base = CyclotomicMonomial(res.base.sign * phase, res.base.q_pow, res.base.phi_exps, res.base.max_d)
    return DCR(res.root, res.radical, new_base, res.ratios, res.z_range, res.max_d)
end

"""
    gsymbol_dcr(j1, j2, j3, j4, j5, j6)
Fully symmetric G-symbol: sqrt(product of all 6 dims) * {6j}.
"""
function gsymbol_dcr(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    dcr = q6j_dcr(j1, j2, j3, j4, j5, j6)
    dcr.base.sign == 0 && return dcr
    
    buf = CycloBuffer(round(Int, 2 * max(j1,j2,j3,j4,j5,j6) + 1))
    for j in (j1, j2, j3, j4, j5, j6)
        add_qint!(buf, round(Int, 2j+1))
    end
    dims = snapshot(buf)
    
    return fuse_prefactors(dcr, dims)
end

# ==============================================================================
# 6. GRAPH EVALUATORS (Theta & Tetrahedron Values)
# ==============================================================================

"""
    theta_dcr(a, b, c)
Evaluates the value of the 'Theta' graph (two vertices connected by 3 edges).
Equivalent to a quantum dimension calculation for the triad.
"""
function theta_dcr(a::Spin, b::Spin, c::Spin)
    !δ(a, b, c) && return ZERO_MONOMIAL
    buf = CycloBuffer(round(Int, a+b+c+1))
    qtriangle!(buf, a, b, c)
    # Norm factor for Theta graph in SU(2)k
    add_qfact!(buf, round(Int, a+b+c+1))
    add_qfact!(buf, round(Int, a+b-c), -1)
    add_qfact!(buf, round(Int, a-b+c), -1)
    add_qfact!(buf, round(Int, -a+b+c), -1)
    return snapshot(buf)
end

"""
    tetrahedron_dcr(j1, j2, j3, j4, j5, j6)
Evaluates the value of a closed tetrahedral network (6j symbol * triangle dims).
"""
function tetrahedron_dcr(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    dcr = q6j_dcr(j1, j2, j3, j4, j5, j6)
    dcr.base.sign == 0 && return dcr
    
    # Value = {6j} * (θ(j1,j2,j3)θ(j1,j5,j6)θ(j2,j4,j6)θ(j3,j4,j5))^(1/2)
    # We build the product of the four theta values
    m_max = round(Int, j1+j2+j3+j4+j5+j6)
    buf = CycloBuffer(m_max)
    # Inlining theta prefactors for speed
    qtetrahedron!(buf, j1, j2, j3, j4, j5, j6)
    # This is essentially the inverse of the qtetrahedron! prefactor in q6j
    # Resulting in the closed graph value.
    return fuse_prefactors(dcr, snapshot(buf))
end