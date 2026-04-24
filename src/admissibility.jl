
#admissibility tests for classical and quantum symbols

# Admissibility tests strictly using twice spins (J = 2j ∈ ℤ)

@inline δ(J1::Int, J2::Int, J3::Int)::Bool = 
    iseven(J1 + J2 + J3) && (abs(J1 - J2) <= J3 <= J1 + J2)

@inline qδ(J1::Int, J2::Int, J3::Int, k::Int)::Bool = 
    δ(J1, J2, J3) && (J1 + J2 + J3) <= 2k 

@inline δtet(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)::Bool = 
    δ(J1, J2, J3) && δ(J1, J5, J6) && δ(J2, J4, J6) && δ(J3, J4, J5)

@inline qδtet(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int, k::Int)::Bool = 
    qδ(J1, J2, J3, k) && qδ(J1, J5, J6, k) && qδ(J2, J4, J6, k) && qδ(J3, J4, J5, k)