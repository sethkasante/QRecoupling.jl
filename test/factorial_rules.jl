# References here build factorials directly, without using rule ratios or DCR projections.
function direct_factorial_rule(s, q)
    qn(n) = q == 1 ? one(q)*n : sum(q^(n-1-2i) for i in 0:n-1; init=zero(q))
    fact(n) = prod(qn(i) for i in 1:n; init=one(q))
    pre = prod(fact(Int(n))^Int(c) for (n,c) in s.pre; init=one(q))
    s.sqrt_pre && (pre=sqrt(pre))
    return s.sign0*pre*sum((s.alternating && isodd(z) ? -1 : 1)*
                          prod(fact(Int(f.a)*z+Int(f.b))^Int(f.c) for f in s.fac; init=one(q))
                          for z in s.zlo:s.zhi; init=zero(q))
end

@testset "General factorial rules" begin
    QR=QRecoupling
    s=qseries([(2,0,1)],1:3)
    @test s isa FactorialSum
    @test qeval(s;q=1) ≈ 746
    @test Float64(qeval(s;q=1,exact=true)) == 746
    @test qeval(Classical(),s) ≈ 746
    @test qeval(qseries([(-2,8,1)],1:3);q=1) ≈ 746
    @test qeval(qseries([(0,3,2)],-2:2);q=1) ≈ 180
    @test qeval(qseries([(2,0,-1)],1:3);q=1) ≈ 391/720
    @test qeval(FactorialSum(-2:2;alternating=true,sign=-1);q=1) == -1
    @test qeval(FactorialSum(2:1);q=1) == 0
    @test qeval(FactorialSum(1:3;sign=0);k=5) == 0
    cancelled=FactorialSum(1:4;factors=[(2,0,3),(2,0,-3),(0,3,0)],prefactor=[3=>2,3=>-2])
    @test isempty(cancelled.fac) && isempty(cancelled.pre)
    @test qeval(cancelled;q=1) == 4
    @test_throws DomainError FactorialSum(0:3;factors=[(2,-1,1)])
    @test_throws DomainError FactorialSum(0:3;prefactor=[-1=>1])
    @test_throws ArgumentError FactorialSum(0:3;sign=2)
    @test_throws ArgumentError qeval(s)
    @test_throws ArgumentError qeval(s;k=5,q=1)
    @test_throws DomainError qeval(s;k=-1)

    rng=MersenneTwister(78456)
    work=EvaluationWorkspace()
    old=QR.POLICY[]
    try
        for i in 1:100
            lo,hi=-1,rand(rng,2:5)
            factors=map(1:rand(rng,1:5)) do _
                a=rand(rng,-3:3)
                b=rand(rng,0:3)-min(a*lo,a*hi)
                (a,b,rand(rng,(-2,-1,1,2)))
            end
            r=FactorialSum(lo:hi;factors=factors,prefactor=[3=>rand(rng,(-1,1))],
                           sqrt_prefactor=isodd(i),alternating=isodd(i÷2),sign=isodd(i÷3) ? -1 : 1)
            k=2QR.max_argument(r)+10
            for mode in (:default,:lazy,:compensated_only)
                QR.POLICY[]=mode
                ref=setprecision(BigFloat,512) do
                    direct_factorial_rule(r,BigFloat(1))
                end
                @test isapprox(qeval(r;q=1,workspace=work),ref;rtol=2e-12,atol=1e-25)
                refk=setprecision(BigFloat,512) do
                    direct_factorial_rule(r,cispi(big(1)/(k+2)))
                end
                @test isapprox(qeval(r;k=k,workspace=work),real(refk);rtol=2e-12,atol=1e-25)
            end
            if i<=20
                q=0.8+0.3im
                @test isapprox(qeval(r;q=q),direct_factorial_rule(r,q);rtol=1e-10,atol=1e-20)
                @test Float64(qeval(r;q=1,exact=true)) ≈ Float64(direct_factorial_rule(r,BigFloat(1))) rtol=1e-12 atol=1e-25
            end
        end
    finally
        QR.POLICY[]=old
    end

    # Classifier versus literal valuation at each z, including slopes larger than h.
    for i in 1:500
        lo,hi=0,rand(rng,1:30)
        a=rand(rng,-12:12)
        b=rand(rng,0:12)-min(a*lo,a*hi)
        r=FactorialSum(lo:hi;factors=[(a,b,rand(rng,(-2,-1,1,2))),(0,3,1)])
        k=rand(rng,0:12); h=k+2
        vs=[QR.term_valuation(r,z,h) for z in lo:hi]
        status,segs=QR.classify_at_level(r,k)
        expected=any(<(0),vs) ? :pole : all(>(0),vs) ? :zero : :finite
        @test status === expected
        status === :finite && @test collect(Iterators.flatten(segs)) == [z for z in lo:hi if vs[z+1]==0]
        small=QR._classify_small(r,k)
        if small!==nothing
            @test small[1] === status
            status === :finite && @test collect(small[2]) == segs
        end
    end

    # Leading specialized zero must not discard the later nonzero term.
    d=qseries(z->z==1 ? qint(5) : qint(1),1:2)
    @test QR.evaluate_exact(qeval(d;k=3,exact=true)) ≈ 1
    r=qseries([(-2,7,1)],1:3) # [5]! + [3]! + [1]! at k=3
    @test qeval(r;k=3) ≈ real(QR.evaluate_exact(qeval(r;k=3,exact=true)))
    r=qseries([(0,6,1),(0,5,-1)],1:3) # cancels [5]! before specialization
    @test qeval(r;k=3) ≈ -3
    @test QR.evaluate_exact(qeval(r;k=3,exact=true)) ≈ -3
    @test_throws DomainError qeval(qseries([(0,5,-1)],1:2);k=3)
    @test_throws DomainError qeval(qseries([(0,5,-1)],1:2);k=3,exact=true)
end

@testset "Validated intervals and owned scratch" begin
    QR=QRecoupling
    work=EvaluationWorkspace()
    for k in 0:6, js in all_6j(k=k)
        J=QR.doubled(js...)
        for (f,rule,family) in ((q6j,QR.sixj_sum,Val(:sixj)),(fsymbol,QR.fsymbol_sum,Val(:f)),(gsymbol,QR.gsymbol_sum,Val(:g)))
            s=rule(J...)
            status,segs=QR.classify_at_level(s,k)
            @test status===:finite && segs==[s.zlo:min(s.zhi,k)]
            @test f(js...;k=k,workspace=work) == f(js...;k=k)
        end
    end
    for js in ((1,1,1,1,1,1),(30,30,30,30,30,30),(100,100,100,100,100,100)), k in (300,500)
        @test q6j(js...;k=k,workspace=work) == q6j(js...;k=k)
        @test q6j(js...;q=1,workspace=work) == q6j(js...;q=1)
    end
    for f in (q6j,fsymbol,gsymbol), T in (Float64,BigFloat)
        L=[(2,2,2,2,2,2),(0,0,0,0,0,0),(1,1,1,1,1,1),(0,0,1,0,0,1)]
        @test f(L;k=3,T=T,threads=4) == [f(l...;k=3,T=T) for l in L]
    end
    L3=[(1,1,1,0,0,0),(2,2,2,0,0,0),(1,1,1,2,-2,0)]
    @test q3j(L3;k=3,threads=4)==[q3j(l...;k=3) for l in L3]
    @test q3j(1,1,1,1,-1,0;k=10,workspace=work)==q3j(1,1,1,1,-1,0;k=10)
    L=repeat([(30,30,30,30,30,30),(50,50,50,50,50,50),(100,100,100,100,100,100)],20)
    @test q6j(L;k=500,threads=4)==q6j(L;k=500,threads=1)
    @test_throws ArgumentError q6j(1,1,1,1,1,1;k=3:5,workspace=work)
end
