# src/Symmetry.jl

"""
    canonical_spins(j1, j2, j3, j4, j5, j6) -> NTuple{6, Float64}

Exploits the 24-fold tetrahedral symmetry (S_4) of the 6j symbol to return a 
unique, lexicographically maximal tuple of spins. 
"""
@inline function canonical_spins(j1, j2, j3, j4, j5, j6)
    # Fast-path for the regular tetrahedron
    if all(x -> x == j1, (j2, j3, j4, j5, j6))
        return (Float64(j1), Float64(j2), Float64(j3), Float64(j4), Float64(j5), Float64(j6))
    end
    
    # Ensure uniform types for the tuple to avoid type-instability in caching
    t = (Float64(j1), Float64(j2), Float64(j3), Float64(j4), Float64(j5), Float64(j6))
    
    # 1. Base Column Permutations (S_3)
    p1 = t
    p2 = (t[2], t[1], t[3], t[5], t[4], t[6])
    p3 = (t[3], t[2], t[1], t[6], t[5], t[4])
    p4 = (t[1], t[3], t[2], t[4], t[6], t[5])
    p5 = (t[2], t[3], t[1], t[5], t[6], t[4])
    p6 = (t[3], t[1], t[2], t[6], t[4], t[5])

    # 2. Helper closure for row-flips (Z_2 x Z_2)
    @inline flips(x) = (
        x,
        (x[4], x[5], x[3], x[1], x[2], x[6]), # Flip cols 1 & 2
        (x[1], x[5], x[6], x[4], x[2], x[3]), # Flip cols 2 & 3
        (x[4], x[2], x[6], x[1], x[5], x[3])  # Flip cols 1 & 3
    )

    # 3. Generate all 24 configurations and return the max
    all_perms = (
        flips(p1)..., flips(p2)..., flips(p3)..., 
        flips(p4)..., flips(p5)..., flips(p6)...
    )

    return reduce(max, all_perms)
end