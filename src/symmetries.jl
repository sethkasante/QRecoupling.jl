
#  canonical symmetries

@inline function canonical_spins(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    if all(x -> x == j1, (j2, j3, j4, j5, j6))
        return (Float64(j1), Float64(j2), Float64(j3), Float64(j4), Float64(j5), Float64(j6))
    end
    
    t = (Float64(j1), Float64(j2), Float64(j3), Float64(j4), Float64(j5), Float64(j6))
    p1 = t; p2 = (t[2], t[1], t[3], t[5], t[4], t[6]); p3 = (t[3], t[2], t[1], t[6], t[5], t[4])
    p4 = (t[1], t[3], t[2], t[4], t[6], t[5]); p5 = (t[2], t[3], t[1], t[5], t[6], t[4]); p6 = (t[3], t[1], t[2], t[6], t[4], t[5])

    @inline flips(x) = (
        x,
        (x[4], x[5], x[3], x[1], x[2], x[6]), 
        (x[1], x[5], x[6], x[4], x[2], x[3]), 
        (x[4], x[2], x[6], x[1], x[5], x[3])  
    )

    all_perms = (flips(p1)..., flips(p2)..., flips(p3)..., flips(p4)..., flips(p5)..., flips(p6)...)
    return reduce(max, all_perms)
end 

#TODO: implement symmetries of 3j and other Regge symmetries?