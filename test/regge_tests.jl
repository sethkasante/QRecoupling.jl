# The 144 symmetries of the 6j symbol (tetrahedral relabellings times Regge's) and `regge_canonical`.

"Doubled labels with triangle sums α and quadrilateral sums β in the given positions (the inverse of `racah_sums`)."
_labels_from_sums(α, β) = (α[1] + α[2] - β[3], α[1] + α[3] - β[2], α[1] + α[4] - β[1],
                           α[3] + α[4] - β[3], α[2] + α[4] - β[2], α[2] + α[3] - β[1])
_perms(v) = length(v) <= 1 ? [collect(v)] : [vcat(v[i], p) for i in eachindex(v) for p in _perms(deleteat!(collect(v), i))]

@testset "Regge symmetries" begin
    samples = all_6j(k = 6)[1:41:end]
    orbit_sizes = Int[]
    for js in samples
        J = QR.doubled(js...)
        α, β = QR.racah_sums(J...)
        c = QR.regge_canonical(J...)
        @test QR._δtet(c...) && QR.racah_sums(c...) == (α, β)
        @test QR.regge_canonical(c...) == c   v# idempotent
        k = maximum(α) + 1     # the smallest admissible level
        v_q, v_k, v_e = q6j(js...; q = 0.8), q6j(js...; k = k), q6j(Exact(), js...)
        images = Set(_labels_from_sums(α[pa], β[pb]) for pa in _perms(1:4) for pb in _perms(1:3))
        push!(orbit_sizes, length(images))
        @test all(I -> QR._δtet(I...) && QR.regge_canonical(I...) == c, images)
        @test all(I -> rel(q6j((I .// 2)...; q = 0.8), v_q) < 1e-13, images)
        @test all(I -> rel(q6j((I .// 2)...; k = k), v_k) < 1e-13, images)
        @test all(I -> q6j(Exact(), (I .// 2)...) == v_e, images)
    end
    @test maximum(orbit_sizes) > 24     # Regge adds to the tetrahedral orbit
    # a canonical sweep keeps one label set per class, and every label set has its class represented
    k = 4
    canon(l) = QR.regge_canonical(QR.doubled(l...)...)
    classes = Set(canon(l) for l in all_6j(k = k))
    reps = all_6j(k = k, canonical = true)
    @test length(reps) == length(classes)
    @test Set(canon(l) for l in reps) == classes
end
