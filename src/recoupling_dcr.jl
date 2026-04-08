
# --------------------------------------------------------------
#  -- DCR Builders for Quantum Recoupling Symbols --
# --------------------------------------------------------------


# --- Atomic units: Dimensions & Phases ---


"""
    qint_mono(n::Int) -> CyclotomicMonomial
Returns the algebraic representation of the quantum integer [n]q.
"""
function qint_mono(n::Int)
    n <= 0 && return ZERO_MONOMIAL
    n == 1 && return ONE_MONOMIAL
    buf = CycloBuffer(n)
    add_qint!(buf, n, 1)
    return snapshot(buf)
end

@inline qdim_mono(j::Spin) = qint_mono(round(Int, 2j + 1))

#qfact_mono

"""
    rmatrix_mono(j1, j2, j3) -> CyclotomicMonomial
The algebraic phase for the R-matrix braiding j1, j2 into j3.
Uses the q^{1/4} power convention to securely encode quarter-integer 
phases as exact integers, bypassing float precision loss.
"""
function rmatrix_mono(j1::Spin, j2::Spin, j3::Spin)
    # The true physical exponent is: 1/2 * (c3 - c1 - c2)
    # q_pow now tracks powers of q^{1/4}
    p = round(Int, 2 * (j3*(j3+1) - j1*(j1+1) - j2*(j2+1))) 
    
    # parity sign
    s = iseven(round(Int, j1 + j2 - j3)) ? 1 : -1
    
    return CyclotomicMonomial(Int8(s), p, Pair{Int,Int}[], 0)
end


# --- buffers for triangle coefficients ----

@inline function qtriangle!(buf, j1, j2, j3)
    add_qfact!(buf, round(Int, j1 + j2 - j3))
    add_qfact!(buf, round(Int, j1 - j2 + j3))
    add_qfact!(buf, round(Int, -j1 + j2 + j3))
    add_qfact!(buf, round(Int, j1 + j2 + j3 + 1), -1)
end

@inline function qtetrahedron!(buf, j1, j2, j3, j4, j5, j6)
    qtriangle!(buf, j1, j2, j3)
    qtriangle!(buf, j1, j5, j6)
    qtriangle!(buf, j2, j4, j6)
    qtriangle!(buf, j3, j4, j5)
end

# --- Recoupling Symbols (3j & 6j) ---

function q3j_dcr(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    # admissibile conditions
    (!δ(j1, j2, j3) || m1 + m2 + m3 != 0) && return ZERO_DCR 

    # Standard summation bounds for Wigner 3j
    α = (round(Int, j3 - j2 + m1), round(Int, j3 - j1 - m2))
    β = (round(Int, j1 + j2 - j3), round(Int, j1 - m1), round(Int, j2 + m2))
    z_min, z_max = max(0, -α[1], -α[2]), min(β[1], β[2], β[3])
    
    # Initialize buffer 
    buf = CycloBuffer(max(z_max + 2, round(Int, j1 + j2 + j3 + 1)))

    return build_dcr!(buf,
        # Prefactor
        b -> begin
            qtriangle!(b, j1, j2, j3)
            # Add (j ± m)! terms explicitly 
            add_qfact!(b, round(Int, j1 + m1)); add_qfact!(b, round(Int, j1 - m1))
            add_qfact!(b, round(Int, j2 + m2)); add_qfact!(b, round(Int, j2 - m2))
            add_qfact!(b, round(Int, j3 + m3)); add_qfact!(b, round(Int, j3 - m3))
        end,

        # Base Term (at z_min)
        (b, z) -> (add_qfact!(b, z, -1); 
                    map(a -> add_qfact!(b, a+z, -1), α); 
                    map(bv -> add_qfact!(b, bv-z, -1), β)
                   ),
        
        # Ratio R_z = Term(z+1)/Term(z)
        (b, z) -> (map(bv -> add_qint!(b, bv-z), β); 
                    add_qint!(b, z+1, -1); 
                    map(a -> add_qint!(b, a+z+1, -1), α)
                   ),
        z_min, z_max;
        extract_radical = true,
        alternating_sign = true
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
        # base term 
        (b, z) -> begin 
                    add_qfact!(b, z+1); 
                    map(a -> add_qfact!(b, z-a, -1), α); 
                    map(bv -> add_qfact!(b, bv-z, -1), β)
                end,

        # ratios 
        (b, z) -> begin 
                    add_qint!(b, z+2); 
                    map(bv -> add_qint!(b, bv-z), β); 
                    map(a -> add_qint!(b, z+1-a, -1), α)
                end,
        z_min, z_max;
        extract_radical = true,
        alternating_sign = true
    )
end



#  ---- TQFT invariants (F & G Symbols) ---

"""
    fsymbol_dcr(j1, j2, j3, j4, j5, j6)
Algebraically fused F-symbol: sqrt([2j3+1][2j6+1]) * {6j}, fuses √([d3][d6]) into the radical.
Matches the unitary crossing matrix in Spin Networks.
"""
function fsymbol_dcr(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    dcr = q6j_dcr(j1, j2, j3, j4, j5, j6)
    dcr.base.sign == 0 && return dcr
    
    dims = qdim_mono(j3) * qdim_mono(j6)
    res = fuse_radical(dcr, dims)
    
    #phase factor 
    phase = iseven(round(Int, j1 + j2 + j4 + j5)) ? 1 : -1
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
    
    return fuse_radical(dcr, dims)
end


# ---- Graph evaluators (Theta & Tetrahedron Values)  ----

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
    return fuse_radical(dcr, snapshot(buf))
end