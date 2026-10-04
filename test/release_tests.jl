# Representative regressions for previously fixed bugs; broad sweeps stay outside CI.
@testset "release regressions" begin
    @testset "modular candidates need exact confirmation" begin
        # Deliberately collide in small valid fields, including a ratio beyond machine Int.
        s = FactorialSum(0:1; factors=((1,1,120),), alternating=true)
        tabs = (QR.ClassicalModTable(UInt64(5),2), QR.ClassicalModTable(UInt64(13),2))
        @test QR._sum_vanishes_mod(s,tabs)
        @test qeval(Exact(),s) == 1-big(2)^120
        @test !QR.is_classical_zero(s,tabs)
        s = FactorialSum(0:1; factors=((1,1,6),), alternating=true)
        tabs = (QR._generic_mod_table(UInt64(7),2,UInt64(2)),
                QR._generic_mod_table(UInt64(13),2,UInt64(2)))
        @test QR._sum_vanishes_mod(s,tabs) && !QR.is_generic_zero(s,tabs)

        # k=4: 1-[2]!^6=-26, but both conjugates vanish modulo 13.
        p = 13; m = QR.Montgomery(p)
        root = Int(QR.root_of_unity(UInt64(p),12))
        facts = Matrix{UInt64}(undef,6,2); invfacts = similar(facts)
        for (c,a) in enumerate((1,5))
            q = powermod(root,a,p); qi = invmod(q,p); f = 1
            for n in 0:5
                n > 0 && (f = mod(f*(powermod(q,n,p)-powermod(qi,n,p))*invmod(q-qi,p),p))
                facts[n+1,c] = QR.to_mont(m,f)
                invfacts[n+1,c] = QR.to_mont(m,invmod(f,p))
            end
        end
        tab = QR.LevelZeroTable(m,facts,invfacts)
        @test QR.is_cancellation_zero(s,(0:1,),4,tab) === true
        @test qeval(Exact(4),s) == -26
        @test !QR.is_zero_at_level(s,4,tab)
        @test QR.level_escalate(s,(0:1,),4,Float64,tab) ≈ -26 rtol=1e-14
    end

    @testset "zeros at the supplied q, with and without weights" begin
        s = FactorialSum(0:2; factors=((1,0,-1),(-1,2,-1)), alternating=true)
        for q in (0.8,0.7+0.4im)
            @test !QR._analytic_sum_iszero(s,q,0) # 2/[2]-1 is nonzero here
            @test iszero(QR.analytic_value(s,q;weight=(1,0))) # (1+q²)/[2]-q = 0
        end
        @test QR._analytic_sum_iszero(s,1,0) # an isolated zero, not a zero function
        q = 1+big(1)//big(2)^500
        @test !QR._analytic_sum_iszero(s,q,0)
        @test qeval(s;q) ≈ BigFloat(2/(q+inv(q))-1) rtol=big"1e-60"
        pair = FactorialSum(0:1; alternating=true)
        @test QR._analytic_sum_iszero(pair,0.8,0) && !QR._analytic_sum_iszero(pair,0.8,1)

        # A collision, a bad rational denominator, a zero complex image, and a zero q-integer.
        s = FactorialSum(0:1; factors=((1,1,6),), alternating=true)
        m = QR.Montgomery(13); ii = QR.to_mont(m,5)
        for q in (2,1//13,-5+im,5)
            r = QR._analytic_rational_q(q)
            @test !QR._analytic_sum_nonzero_mod(s,r,0,m,ii)
            @test !QR._analytic_sum_iszero(s,q,0,m,ii)
        end
    end

    @testset "parallel zero flags and genuine cancellations" begin
        cancelled = (1,3//2,3//2,2,3//2,3//2)
        @test iszero(q6j(Exact(10),cancelled...)) && iszero_at(10,cancelled...)
        # Cross the packed Boolean word boundary once, rather than repeating stress sweeps.
        cases = [cancelled,(1,1,1,1,1,1),(1,1,1,1,1,4)]
        flags = [iszero(q6j(Exact(10),l...)) for l in cases]
        labels = [cases[mod1(i,3)] for i in 1:65]
        for input in (empty(labels),labels), nw in (1,4)
            actual = iszero_at(10,input; threads=nw)
            @test actual isa BitVector
            @test actual == [flags[mod1(i,3)] for i in eachindex(input)]
        end
    end

    @testset "level_spectrum: screened candidates, and proofs on request" begin
        for l in ((5,5,5,5,5,5), (1,1,1,5//2,5//2,5//2), (2,3//2,3//2,9//2,4,5)), nw in (1,4)
            screen = level_spectrum(l...; k=2:40, threads=nw)
            proved = level_spectrum(l...; k=2:40, prove=true, threads=nw)
            # proving changes only :cancels entries, and a proved :cancels is an exact zero
            @test all(screen[i] === proved[i] || screen[i] === :cancels for i in eachindex(screen))
            for (i, kk) in enumerate(2:40)
                proved[i] === :cancels && @test iszero(q6j(Exact(kk), l...))
                proved[i] === :finite && @test !iszero_at(kk, l...)
            end
        end
    end

    @testset "exact equality survives small rational scaling" begin
        a = q6j(Exact(4),0,0,0,1//2,1//2,1//2)
        b = q6j(Exact(4),0,0,0,3//2,3//2,3//2)
        @test a.sqclass != b.sqclass && a == b
        tiny = big(1)//big(10)^100
        @test iszero(tiny*(a-b)) && !iszero(tiny*(a+b))
        @test (a-b)+tiny == tiny
    end

    @testset "tables retain numeric type and precision" begin
        smatrix(11) # Float64 first at 11; BigFloat first at 13
        for bits in (64,128)
            setprecision(BigFloat,bits) do
                for k in (11,13)
                    S, = smatrix(k; T=BigFloat)
                    @test eltype(S) === BigFloat && all(x->precision(x)==bits,S)
                    @test S[1,1] ≈ sqrt(BigFloat(2)/(k+2))*sin(big(pi)/(k+2)) rtol=eps(BigFloat)*4
                end
            end
        end
        @test eltype(first(smatrix(13))) === Float64
    end

    @testset "F matrices preserve scalar values and input precision" begin
        setprecision(BigFloat,256) do
            for kw in ((;k=6,T=BigFloat), (;q=big"0.8",T=Float64),
                       (;q=big"-0.8",T=Float64), (;q=0.8+0.3im,T=BigFloat))
                E = haskey(kw,:q) && !(kw.q isa Real && kw.q>0) ? Complex{BigFloat} : BigFloat
                F, es, fs = fmatrix(1,1,1,1; kw...)
                ref = [fsymbol(1,1,e,1,1,f; kw...,T=BigFloat) for e in es,f in fs]
                identity = [i==j ? one(E) : zero(E) for i in eachindex(fs),j in eachindex(fs)]
                @test eltype(F) === E
                @test maximum(abs,F-ref) < big"1e-65"
                @test maximum(abs,transpose(F)*F-identity) < big"1e-65"
            end
        end
    end

    @testset "negative-q square-root convention" begin
        l = (1//2,1//2,0,1//2,1//2,1); q = -0.8
        for f in (q6j,fsymbol,gsymbol)
            v = f(l...;q)
            r = setprecision(() -> f(l...;q=Complex{BigFloat}(q)),BigFloat,256)
            @test v == f(l...;q=complex(q))
            @test v ≈ r rtol=1e-13
            @test qeval(f(Symbolic(),l...);q) == v
        end
        @test real(fsymbol(l...;q)) < 0
        @test project_analytic(fsymbol(Symbolic(),l...).dcr,q) ≈ fsymbol(l...;q) rtol=1e-13
        @test_logs q6j(1,1,1,1,1,1;q=-nextfloat(1.0))
    end

    @testset "invalid error estimates cannot certify a value" begin
        for (v,b) in ((Inf,0.0),(NaN,0.0),(1.0,-1.0),(1.0,Inf),(0.0,0.0))
            @test !QR._certifies(v,b,1e-12)
            @test QR._surviving_digits(v,b,1) == -Inf
        end
        @test QR._certifies(1.0,1e-14,1e-12)
    end
end
