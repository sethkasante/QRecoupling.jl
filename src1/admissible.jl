
# src/admissible.jl

"""
    ishalfInt(j) -> Bool
"""
@inline ishalfInt(j)::Bool = isinteger(2*j) && j ≥ 0

"""
    δ(j1, j2, j3) -> Bool

Classical angular momentum triangle inequality.
"""
@inline δ(j1, j2, j3)::Bool = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)

"""
    qδ(j1, j2, j3, k) -> Bool

Quantum triangle inequality: satisfies classical δ and sum j1+j2+j3 ≤ k.
"""
@inline qδ(j1, j2, j3, k)::Bool = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 

"""
    δtet(j1, j2, j3, j4, j5, j6) -> Bool

Admissibility condition for a classical tetrahedral (6j symbol) configuration.
"""
@inline δtet(j1, j2, j3, j4, j5, j6)::Bool = 
    ishalfInt(j1) && ishalfInt(j2) && ishalfInt(j3) && 
    ishalfInt(j4) && ishalfInt(j5) && ishalfInt(j6) && 
    δ(j1, j2, j3) && δ(j1, j5, j6) && δ(j2, j4, j6) && δ(j3, j4, j5)

"""
    qδtet(j1, j2, j3, j4, j5, j6, k) -> Bool

Admissibility condition for a quantum tetrahedral configuration at level k.
"""
@inline qδtet(j1, j2, j3, j4, j5, j6, k)::Bool = 
    ishalfInt(j1) && ishalfInt(j2) && ishalfInt(j3) && 
    ishalfInt(j4) && ishalfInt(j5) && ishalfInt(j6) && 
    qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)

