@testset "Regge symmetry and enumeration" begin
    # A Regge transform preserves the sorted Racah sums and the coefficient.
    labels = (2,3//2,5//2,2,5//2,3//2)
    J = QR.doubled(labels...)
    c = QR.regge_canonical(J...)
    @test QR.regge_canonical(c...) == c
    @test QR.racah_sums(c...) == QR.racah_sums(J...)
    @test q6j(Exact(), labels...) == q6j(Exact(), (c .// 2)...)
    canon(l) = QR.regge_canonical(QR.doubled(l...)...)
    classes = Set(canon(l) for l in all_6j(k=3))
    reps = all_6j(k=3,canonical=true)
    @test length(reps) == length(classes)
    @test Set(canon(l) for l in reps) == classes
end
