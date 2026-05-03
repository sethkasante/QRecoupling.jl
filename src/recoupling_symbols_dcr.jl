
# --------------------------------------------------------------
#  -- DCR Builders for Quantum Recoupling Symbols --
# --------------------------------------------------------------


# --- Atomic units: Dimensions & Phases ---


# quantum dimension [2j+1] -> [J+1]
"""
    qdim_mono(J::Int) -> CyclotomicMonomial
Returns the algebraic quantum dimension [J+1]_q.
Inputs use twice spins (J = 2j).
"""
@inline qdim_mono(J::Int) = qint(J + 1)


#qfact_mono

"""
    rmatrix_mono(J1, J2, J3) -> CyclotomicMonomial
The algebraic phase for the R-matrix braiding J1, J2 into J3.
Uses the q^{1/4} power convention to securely encode quarter-integer 
Inputs are 2*Spins (J = 2j ∈ ℤ).
"""
function rmatrix_mono(J1::Int, J2::Int, J3::Int)
    # The phase is (J3(J3+2) - J1(J1+2) - J2(J2+2)) / 2
    # By SU(2) admissibility, J1+J2-J3 is even, guaranteeing this numerator is even.
    p = (J3*(J3+2) - J1*(J1+2) - J2*(J2+2)) ÷ 2
    
    # parity sign
    s = iseven((J1 + J2 - J3) ÷ 2) ? 1 : -1
    
    return iseven(p) ? CyclotomicMonomial(Int8(s), p ÷ 2, Pair{Int,Int}[], 0) : "R-matrix phase for these spins requires a half-integer power: $(s) * q^($(p)/2). Not stored as CyclotomicMonomial"
end


# --- buffers for triangle coefficients ----

@inline function qtriangle!(buf, J1::Int, J2::Int, J3::Int)
    add_qfact!(buf, (J1 + J2 - J3) ÷ 2)
    add_qfact!(buf, (J1 - J2 + J3) ÷ 2)
    add_qfact!(buf, (-J1 + J2 + J3) ÷ 2)
    add_qfact!(buf, (J1 + J2 + J3) ÷ 2 + 1, -1)
end

@inline function qtetrahedron!(buf, J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    qtriangle!(buf, J1, J2, J3)
    qtriangle!(buf, J1, J5, J6)
    qtriangle!(buf, J2, J4, J6)
    qtriangle!(buf, J3, J4, J5)
end

# --- Recoupling Symbols (3j & 6j) ---

function q3j_dcr(J1::Int, J2::Int, J3::Int, M1::Int, M2::Int, M3::Int = -M1-M2)
    # admissibile conditions
    (!_δ(J1, J2, J3) || M1 + M2 + M3 != 0) && return ZERO_DCR

    # Standard summation bounds for Wigner 3j
    α1 = (J3 - J2 + M1) ÷ 2; α2 = (J3 - J1 - M2) ÷ 2
    β1 = (J1 + J2 - J3) ÷ 2; β2 = (J1 - M1) ÷ 2; β3 = (J2 + M2) ÷ 2
    
    z_min = max(0, -α1, -α2)
    z_max = min(β1, β2, β3)
    
    # Initialize buffer 
    buf = CycloBuffer(max(z_max + 2, (J1 + J2 + J3) ÷ 2 + 1))

    return build_dcr!(buf,
        # Prefactor
        b -> begin
            qtriangle!(b, J1, J2, J3)
            # Add (j ± m)! terms explicitly 
            add_qfact!(b, (J1 + M1) ÷ 2); add_qfact!(b, (J1 - M1) ÷ 2)
            add_qfact!(b, (J2 + M2) ÷ 2); add_qfact!(b, (J2 - M2) ÷ 2)
            add_qfact!(b, (J3 + M3) ÷ 2); add_qfact!(b, (J3 - M3) ÷ 2)
        end,

        # Base Term (at z_min)
        (b, z) -> begin
            add_qfact!(b, z, -1)
            for a in (α1, α2); add_qfact!(b, a+z, -1); end
            for bv in (β1, β2, β3); add_qfact!(b, bv-z, -1); end
        end,
        
        # Ratio R_z = Term(z+1)/Term(z)
        (b, z) -> begin
            for bv in (β1, β2, β3); add_qint!(b, bv-z); end
            add_qint!(b, z+1, -1)
            for a in (α1, α2); add_qint!(b, a+z+1, -1); end
        end,
        z_min, z_max;
        extract_radical = true,
        alternating_sign = true
    )
end

function q6j_dcr(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    !_δtet(J1, J2, J3, J4, J5, J6) && return ZERO_DCR

    α1 = (J1+J2+J3) ÷ 2; α2 = (J1+J5+J6) ÷ 2; α3 = (J2+J4+J6) ÷ 2; α4 = (J3+J4+J5) ÷ 2
    β1 = (J1+J2+J4+J5) ÷ 2; β2 = (J1+J3+J4+J6) ÷ 2; β3 = (J2+J3+J5+J6) ÷ 2
    
    α = (α1, α2, α3, α4)
    β = (β1, β2, β3)
    
    z_min = max(α1, α2, α3, α4)
    z_max = min(β1, β2, β3)
    buf = CycloBuffer(max(z_max + 2, β1, β2, β3))

    return build_dcr!(buf,
        b -> qtetrahedron!(b, J1, J2, J3, J4, J5, J6),
        
        # base term 
        (b, z) -> begin 
            add_qfact!(b, z+1)
            for a in α; add_qfact!(b, z-a, -1); end
            for bv in β; add_qfact!(b, bv-z, -1); end
        end,

        # ratios 
        (b, z) -> begin 
            add_qint!(b, z+2)
            for bv in β; add_qint!(b, bv-z); end
            for a in α; add_qint!(b, z+1-a, -1); end
        end,
        z_min, z_max;
        extract_radical = true,
        alternating_sign = true
    )
end



#  ---- TQFT invariants (F & G Symbols) ---

"""
    fsymbol_dcr(J1, J2, J3, J4, J5, J6)
Algebraically fused F-symbol: sqrt([2j3+1][2j6+1]) * {6j}, fuses √([d3][d6]) into the radical.
Matches the unitary crossing matrix in Spin Networks.
Inputs are twice spins (J = 2j ∈ ℤ).
"""
function fsymbol_dcr(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    dcr = q6j_dcr(J1, J2, J3, J4, J5, J6)
    dcr.base.sign == 0 && return dcr
    
    dims = qdim_mono(J3) * qdim_mono(J6)
    res = fuse_radical(dcr, dims)
    
    # phase factor 
    phase = iseven((J1 + J2 + J4 + J5) ÷ 2) ? 1 : -1
    new_base = phase * res.base
    
    return DCR(res.root, res.radical, new_base, res.ratios, res.z_range, res.max_d)
end


"""
    gsymbol_dcr(J1, J2, J3, J4, J5, J6)
Fully symmetric G-symbol: sqrt(product of all 6 dims) * {6j}.
"""
function gsymbol_dcr(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    dcr = q6j_dcr(J1, J2, J3, J4, J5, J6)
    dcr.base.sign == 0 && return dcr
    
    buf = CycloBuffer(max(J1, J2, J3, J4, J5, J6) + 1)
    for J in (J1, J2, J3, J4, J5, J6)
        add_qint!(buf, J + 1)
    end
    
    return fuse_radical(dcr, snapshot(buf))
end


# ---- Graph evaluators (Theta & Tetrahedron Values)  ----

"""
    theta_mono(A, B, C)
Evaluates the value of the 'Theta' graph (two vertices connected by 3 edges).
Equivalent to a quantum dimension calculation for the triad.
"""
function theta_mono(A::Int, B::Int, C::Int)
    !_δ(A, B, C) && return ZERO_MONOMIAL
    buf = CycloBuffer((A + B + C) ÷ 2 + 1)
    qtriangle!(buf, A, B, C)
    # Norm factor for Theta graph in SU(2)k
    add_qfact!(buf, (A + B + C) ÷ 2 + 1)
    add_qfact!(buf, (A + B - C) ÷ 2, -1)
    add_qfact!(buf, (A - B + C) ÷ 2, -1)
    add_qfact!(buf, (-A + B + C) ÷ 2, -1)
    return snapshot(buf)
end

"""
    tetrahedron_dcr(J1, J2, J3, J4, J5, J6)
Evaluates the value of a closed tetrahedral network (6j symbol * triangle dims).

"""
function tetrahedron_dcr(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    dcr = q6j_dcr(J1, J2, J3, J4, J5, J6)
    dcr.base.sign == 0 && return dcr
    
    # Value = {6j} * (θ(j1,j2,j3)θ(j1,j5,j6)θ(j2,j4,j6)θ(j3,j4,j5))^(1/2)
    # We build the product of the four theta values
    m_max = (J1 + J2 + J3 + J4 + J5 + J6) ÷ 2
    
    # inline theta prefactors 
    buf = CycloBuffer(m_max)
    qtetrahedron!(buf, J1, J2, J3, J4, J5, J6)

    # This is essentially the inverse of the qtetrahedron! prefactor in q6j
    # Resulting in the closed graph value.
    return fuse_radical(dcr, snapshot(buf))
end
