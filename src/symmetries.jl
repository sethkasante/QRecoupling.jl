
#  canonical symmetries

# @inline function canonical_spins(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
#     if all(x -> x == j1, (j2, j3, j4, j5, j6))
#         return (Float64(j1), Float64(j2), Float64(j3), Float64(j4), Float64(j5), Float64(j6))
#     end
    
#     t = (Float64(j1), Float64(j2), Float64(j3), Float64(j4), Float64(j5), Float64(j6))
#     p1 = t; p2 = (t[2], t[1], t[3], t[5], t[4], t[6]); p3 = (t[3], t[2], t[1], t[6], t[5], t[4])
#     p4 = (t[1], t[3], t[2], t[4], t[6], t[5]); p5 = (t[2], t[3], t[1], t[5], t[6], t[4]); p6 = (t[3], t[1], t[2], t[6], t[4], t[5])

#     @inline flips(x) = (
#         x,
#         (x[4], x[5], x[3], x[1], x[2], x[6]), 
#         (x[1], x[5], x[6], x[4], x[2], x[3]), 
#         (x[4], x[2], x[6], x[1], x[5], x[3])  
#     )

#     all_perms = (flips(p1)..., flips(p2)..., flips(p3)..., flips(p4)..., flips(p5)..., flips(p6)...)
#     return reduce(max, all_perms)
# end 

@inline function canonical_spins(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    # convert to doubled spins (J = 2j) 
    t = doubled(j1, j2, j3, j4, j5, j6)
    # t = (J1, J2, J3, J4, J5, J6)

    if allequal(t)
        return t
    end
    
    p1 = t
    p2 = (t[2], t[1], t[3], t[5], t[4], t[6])
    p3 = (t[3], t[2], t[1], t[6], t[5], t[4])
    p4 = (t[1], t[3], t[2], t[4], t[6], t[5])
    p5 = (t[2], t[3], t[1], t[5], t[6], t[4])
    p6 = (t[3], t[1], t[2], t[6], t[4], t[5])

    @inline flips(x) = (
        x,
        (x[4], x[5], x[3], x[1], x[2], x[6]), 
        (x[1], x[5], x[6], x[4], x[2], x[3]), 
        (x[4], x[2], x[6], x[1], x[5], x[3])  
    )

    all_perms = (flips(p1)..., flips(p2)..., flips(p3)..., flips(p4)..., flips(p5)..., flips(p6)...)
    
    # Returns an NTuple{6, Int} of doubled spins
    return reduce(max, all_perms) 
end



"""
    regge_canonical(J1, J2, J3, J4, J5, J6) -> NTuple{6,Int}

Doubled labels of one fixed representative of the class of a 6j symbol under all 144 symmetries (the 24
tetrahedral relabellings times Regge's 6). The value depends only on the sorted triangle sums α and
quadrilateral sums β (`racah_sums`), and the 144 symmetries are exactly the independent permutations of α and
of β. So placing the sorted sums in fixed positions, `j = (α_i + α_k − β_m)/2`, gives a class invariant.
For admissible labels the result is admissible, at every level where the input is. The G-symbol is not
Regge invariant, because its dimension factors are not; use `canonical_spins` for it.
"""
@inline function regge_canonical(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    α, β = racah_sums(J1, J2, J3, J4, J5, J6)
    return (α[1] + α[2] - β[3], α[1] + α[3] - β[2], α[1] + α[4] - β[1],
            α[3] + α[4] - β[3], α[2] + α[4] - β[2], α[2] + α[3] - β[1])
end


#TODO: implement the symmetries of the 3j symbol (Regge's 72)?



