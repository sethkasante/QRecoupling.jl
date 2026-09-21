
#admissibility tests for all classical and quantum symbols

# Admissibility tests strictly using twice spins (J = 2j ∈ ℤ)

@inline _δ(J1::Int, J2::Int, J3::Int)::Bool =
    iseven(J1 + J2 + J3) && (abs(J1 - J2) <= J3 <= J1 + J2)

@inline _qδ(J1::Int, J2::Int, J3::Int, k::Int)::Bool =
    _δ(J1, J2, J3) && (J1 + J2 + J3) <= 2k

@inline _δtet(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)::Bool =
    _δ(J1, J2, J3) && _δ(J1, J5, J6) && _δ(J2, J4, J6) && _δ(J3, J4, J5)

@inline _qδtet(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int, k::Int)::Bool =
    _qδ(J1, J2, J3, k) && _qδ(J1, J5, J6, k) && _qδ(J2, J4, J6, k) && _qδ(J3, J4, J5, k)

"Magnetic labels for a 3j symbol (doubled): |M| ≤ J with J ± M even, and M1 + M2 + M3 = 0."
@inline _mproj_ok(J1::Int, J2::Int, J3::Int, M1::Int, M2::Int, M3::Int)::Bool =
    M1 + M2 + M3 == 0 && abs(M1) <= J1 && abs(M2) <= J2 && abs(M3) <= J3 &&
    iseven(J1 + M1) && iseven(J2 + M2) && iseven(J3 + M3)


# twice spins (J = 2j); labels that are not multiples of 1/2 are rejected
@inline function doubled(j::Spin)
    J = 2j
    isinteger(J) || throw(ArgumentError("spin labels must be multiples of 1/2, got $j"))
    return Int(J)
end
@inline doubled(js...) = map(doubled, js)


# apis

@inline δ(j1::Spin, j2::Spin, j3::Spin)::Bool =
    _δ(doubled(j1, j2, j3)...)

@inline qδ(j1::Spin, j2::Spin, j3::Spin, k::Int)::Bool =
    _qδ(doubled(j1, j2, j3)..., k)

@inline δtet(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::Bool =
    _δtet(doubled(j1, j2, j3, j4, j5, j6)...)

@inline qδtet(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int)::Bool =
    _qδtet(doubled(j1, j2, j3, j4, j5, j6)..., k)
