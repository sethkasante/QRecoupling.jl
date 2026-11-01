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
    @test s isa QForm && string(s) == "([5] − 1)/[4]!" && s(1) == 1//6
    @test_throws ArgumentError x_form(s)
    @test occursin("x_form", sprint(showerror, try x_form(s) catch e; e end))
    # what needs Nemo still explains itself without it
    VERSION >= v"1.11" && @test all(n -> Docs.hasdoc(QRecoupling, n), filter(!=(:QRecoupling), names(QRecoupling)))
    # exact values are forms and need no backend: only radicals and polynomials do
    v = q6j(Exact(10),1,1,1,1,1,1)
    @test v isa QForm && float(v) ≈ q6j(Level(10),1,1,1,1,1,1)
    @test v * inv(v) == 1 && iszero(v - v)
    @test iszero(q6j(Exact(1),1,1,1,1,1,1))
    @test float(qdim(Exact(10),1)) ≈ qdim(Level(10),1)
    @test q6j(Symbolic(),1,1,1,1,1,1) != 1//6 && q6j(Symbolic(),1,1,1,1,1,1) != q6j(Symbolic(),2,2,2,2,2,2)
    @test float(qcg(Exact(9),3//2,1//2,1,0,3//2,1//2)) ≈ qcg(Level(9),3//2,1//2,1,0,3//2,1//2)
    @test_throws ArgumentError radical(v)
    for q in (1+1e-8, 0.8, 1.3)
        @test QR.qform(q6j,2,2,2,2,2,2)(q) ≈ q6j(At(q),2,2,2,2,2,2) rtol=1e-13
        # `f(q)` evaluates the rule, so test the form as written separately
        written = QR._qform_value(QR.qform(q6j,2,2,2,2,2,2), n -> (big(q)^n - big(q)^(-n)) / (big(q) - 1/big(q)))
        @test Float64(written) ≈ q6j(At(q),2,2,2,2,2,2) rtol=1e-12
    end
    @test QR.qform(qcg,1,0,1,0,1,0)(1.3) ≈ qcg(1,0,1,0,1,0;q=1.3)
    @test float(QR.qform(qcg,3//2,1//2,1,0,3//2,1//2;k=9)) ≈ qcg(3//2,1//2,1,0,3//2,1//2;k=9)
    @test iszero(q6j(Level(10),1,3//2,3//2,2,3//2,3//2))
    @test iszero(QR.qform(q6j,1,5//2,5//2,5//2,3,3;k=10))
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
