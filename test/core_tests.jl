# Run in an environment containing QRecoupling but no Nemo dependency.
using Test, QRecoupling
const QR = QRecoupling

@testset "Nemo-free core" begin
    @test Base.get_extension(QR, :QRecouplingNemoExt) === nothing
    @test !any(id -> id.name == "Nemo", keys(Base.loaded_modules))
    @test q6j(1,1,1,1,1,1) ≈ 1/6
    @test q6j(Exact(),1,1,1,1,1,1) == 1//6
    @test numeric_value(q6j(Exact(),1,1,1,1,1,1)) ≈ big(1)/6
    s = q6j(Symbolic(),1,1,1,1,1,1)
    @test occursin("SymbolicValue", string(s))
    @test_throws ArgumentError x_form(s)
    @test_throws ArgumentError q6j(Exact(10),1,1,1,1,1,1)
    @test_throws ArgumentError q6j(Exact(1),1,1,1,1,1,1)
    @test_throws ArgumentError qdim(Exact(10),1)
    for q in (1+1e-8, 0.8, 1.3)
        @test qform(q6j,2,2,2,2,2,2)(q) ≈ q6j(At(q),2,2,2,2,2,2) rtol=1e-13
        # `f(q)` evaluates the rule, so test the form as written separately
        written = QR._qform_value(qform(q6j,2,2,2,2,2,2), n -> (big(q)^n - big(q)^(-n)) / (big(q) - 1/big(q)))
        @test Float64(written) ≈ q6j(At(q),2,2,2,2,2,2) rtol=1e-12
    end
    @test iszero(q6j(Level(10),1,3//2,3//2,2,3//2,3//2))
    @test iszero(qform(q6j,1,5//2,5//2,5//2,3,3;k=10))
    for k in (3,10,12,30)
        v = QR.level_value(QR._rule_for(q6j,1,1,1,1,1,1),k)
        @test float(v) ≈ q6j(Level(k),1,1,1,1,1,1)
        @test v * inv(v) == 1
        @test !isempty(string(v))
    end
    @test isfinite(abs(qcg(1,0,1,0,0,0;q=0.8)))
    F, es, fs = fmatrix(1,1,1,1;k=6)
    @test transpose(F)*F ≈ [i==j for i in eachindex(fs), j in eachindex(fs)] atol=1e-13
    # Distinct quadratic characters certify a nonzero radical sum; coincident characters do not.
    a = QR.LevelValue(3, [big(1),big(0)], big(1), big(2), Int[])
    b = QR.LevelValue(3, [big(1),big(0)], big(1), big(3), Int[])
    @test QR._distinct_root_classes(QR.LevelValueSum(a) + QR.LevelValueSum(b))
    c = QR.LevelValue(3, [big(1),big(0)], big(1), big(5), Int[])
    g = QR.LevelValue(3, QR._quadratic_roots(QR._level_field(5))[1][2], big(1), big(1), Int[])
    dependent = QR.LevelValueSum(c) - QR.LevelValueSum(g)
    @test !QR._distinct_root_classes(dependent)
    @test iszero(dependent)
    # Modular interpolation tables use batch inversion, with the same exact residues.
    field = QR._level_field(14); t = QR._places(field,1,3)
    for j in 2:field.d, i in j:field.d
        @test QR.mont_mul(t.m,t.nd[j-1][i-j+1],QR.mont_sub(t.m,t.X[i],t.X[i-j+1])) == t.m.one
    end
    # The cache limit controls retained tables, including growth after insertion.
    limit = QR.LEVEL_FIELD_CACHE_LIMIT[]
    try
        QR.LEVEL_FIELD_CACHE_LIMIT[] = 1
        QR._places(field,2,4)
        @test isempty(QR.LEVEL_FIELDS)
    finally
        QR.LEVEL_FIELD_CACHE_LIMIT[] = limit
    end
    @test Base.get_extension(QR, :QRecouplingNemoExt) === nothing
end
