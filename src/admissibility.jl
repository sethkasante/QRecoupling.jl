
#admissibility tests for classical and quantum symbols

@inline ishalfInt(j::Spin)::Bool = isinteger(2*j) && j ≥ 0

@inline δ(j1::Spin, j2::Spin, j3::Spin)::Bool = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)

@inline qδ(j1::Spin, j2::Spin, j3::Spin, k::Int)::Bool = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 

@inline δtet(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::Bool = 
    ishalfInt(j1) && ishalfInt(j2) && ishalfInt(j3) && 
    ishalfInt(j4) && ishalfInt(j5) && ishalfInt(j6) && 
    δ(j1, j2, j3) && δ(j1, j5, j6) && δ(j2, j4, j6) && δ(j3, j4, j5)

@inline qδtet(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int)::Bool = 
    ishalfInt(j1) && ishalfInt(j2) && ishalfInt(j3) && 
    ishalfInt(j4) && ishalfInt(j5) && ishalfInt(j6) && 
    qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)

