
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
@inline qdim_mono(J::Int) = qint_mono(J + 1)


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

The theta net — two vertices joined by three edges coloured `A`, `B`, `C` (doubled labels) — as a
cyclotomic monomial:

    θ(A,B,C) = [T+1]! [m₁]! [m₂]! [m₃]! / ([A]! [B]! [C]!),
    m₁ = (A+B−C)/2,  m₂ = (A−B+C)/2,  m₃ = (−A+B+C)/2,  T = m₁+m₂+m₃ = (A+B+C)/2.

**Correction, replacing a value of 1.** Until this was fixed the "norm factor" below the triangle
coefficient divided by `[m₁]![m₂]![m₃]!` instead of by `[A]![B]![C]!`, which is exactly what
`qtriangle!` had just multiplied in: the two cancelled and `theta_value` returned `1.0` for every
admissible triad, on every evaluation path. Nothing in the package consumed the value, which is why the
suite never saw it.

No convention is being chosen here. Closing a theta net on a single edge gives the quantum dimension,
`θ(A, A, 0) = [A+1]`, and that pins the formula completely once the package's unsigned `qdim` is taken as
the convention — the only remaining freedom is the overall `(−1)^T` of Kauffman–Lins, which is the same
sign their `Δ_a = (−1)^a[a+1]` carries and which `qdim` already drops. The signed variant is `(−1)^T`
times this.
"""
function theta_mono(A::Int, B::Int, C::Int)
    !_δ(A, B, C) && return ZERO_MONOMIAL
    m1 = (A + B - C) ÷ 2; m2 = (A - B + C) ÷ 2; m3 = (-A + B + C) ÷ 2
    t = (A + B + C) ÷ 2 + 1                      # the largest argument, so the buffer needs no more
    buf = CycloBuffer(t)
    add_qfact!(buf, t)
    add_qfact!(buf, m1); add_qfact!(buf, m2); add_qfact!(buf, m3)
    add_qfact!(buf, A, -1); add_qfact!(buf, B, -1); add_qfact!(buf, C, -1)
    return snapshot(buf)
end

"""
    tetrahedron_dcr(J1, J2, J3, J4, J5, J6)
Evaluates the value of a closed tetrahedral network (6j symbol * triangle dims).

"""
tetrahedron_dcr(J1::Int,J2::Int,J3::Int,J4::Int,J5::Int,J6::Int) =
    _factorial_dcr(tetrahedron_sum(J1,J2,J3,J4,J5,J6))
