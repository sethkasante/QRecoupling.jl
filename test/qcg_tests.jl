@testset "quantum coupling" begin
    # Closed-form singlet: both coproduct raising/lowering operators annihilate it.
    for kw in ((; q=0.8), (; q=0.7+0.4im), (; k=8))
        q = haskey(kw,:k) ? cispi(1/(kw.k+2)) : complex(kw.q)
        a = qcg(1//2,1//2,1//2,-1//2,0; kw...)
        b = qcg(1//2,-1//2,1//2,1//2,0; kw...)
        @test a ≈ sqrt(q)/sqrt(q+inv(q)) atol=1e-14
        @test abs(a/sqrt(q)+b*sqrt(q)) < 1e-14
        @test q3j(1//2,1//2,0,1//2,-1//2,0; kw...) ≈ a atol=1e-14
    end
    @test qcg(Exact(), 1//2,1//2,1//2,-1//2,0)^2 == 1//2
    @test qcg((1,1),(1,-1),(1,0)) ≈ 1/sqrt(2)

    # One complete sector and one fusion-truncated sector.
    for kw in ((; q=0.8), (; k=2))
        C, ms, js = qcg_matrix(1,1,0; kw...)
        @test C ≈ [qcg(1,m,1,-m,j; kw...) for m in ms, j in js] atol=1e-13
        @test transpose(C)*C ≈ [i==j for i in eachindex(js), j in eachindex(js)] atol=1e-13
    end
    row, js = qcg_row(1,1,1,-1; q=0.8)
    @test row ≈ [qcg(1,1,1,-1,j; q=0.8) for j in js] atol=1e-13
    @test qcg(1,1,1,0,1,0) == 0 # magnetic labels do not sum to m
    @test_throws ArgumentError q3j(Exact(5), 1,1,1,1,-1,0)
    @test_throws ArgumentError qcg(Symbolic(), 1,1,1,-1,1)
    for f in (qcg,q3j)
        @test_throws ArgumentError f(NTuple{6,Int}[]; k=3,q=0.8)
    end
    # Previously inaccurate large-spin level coupling, including the recurrence fallback.
    l = (709//2, 650//2, 257//2, -425//2, 476//2, -51//2)
    ref = setprecision(() -> q3j(Level(1601; T=BigFloat), l...), BigFloat, 256)
    @test q3j(Level(1601), l...) ≈ ref rtol=1e-13
    cg = setprecision(() -> qcg(Level(1601; T=BigFloat), l[1], l[4], l[2], l[5], l[3]), BigFloat, 256)
    @test QR._cg_entry(nothing, 1601, 709, -425, 650, 476, 257) ≈ cg rtol=1e-13
end
