
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
    rmatrix_mono(J1, J2, J3) -> QPhase
The braiding phase (-1)^{j1+j2-j3} q^{j3(j3+1) - j1(j1+1) - j2(j2+1)} as an exact `QPhase`;
the exponent can be a half-integer. Inputs are doubled spins (J = 2j).
"""
function rmatrix_mono(J1::Int, J2::Int, J3::Int)
    p = (J3*(J3+2) - J1*(J1+2) - J2*(J2+2)) ÷ 2
    s = iseven((J1 + J2 - J3) ÷ 2) ? Int8(1) : Int8(-1)
    return QPhase(s, p // 2)
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

# Symbol formulas live in factorial_rule.jl; these internal methods use doubled labels.
q3j_dcr(J1::Int,J2::Int,J3::Int,M1::Int,M2::Int,M3::Int=-M1-M2) =
    _factorial_dcr(threej_sum(J1,J2,J3,M1,M2,M3))
q6j_dcr(J1::Int,J2::Int,J3::Int,J4::Int,J5::Int,J6::Int) =
    _factorial_dcr(sixj_sum(J1,J2,J3,J4,J5,J6))
fsymbol_dcr(J1::Int,J2::Int,J3::Int,J4::Int,J5::Int,J6::Int) =
    _factorial_dcr(fsymbol_sum(J1,J2,J3,J4,J5,J6))
gsymbol_dcr(J1::Int,J2::Int,J3::Int,J4::Int,J5::Int,J6::Int) =
    _factorial_dcr(gsymbol_sum(J1,J2,J3,J4,J5,J6))

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
tetrahedron_dcr(J1::Int,J2::Int,J3::Int,J4::Int,J5::Int,J6::Int) =
    _factorial_dcr(tetrahedron_sum(J1,J2,J3,J4,J5,J6))
