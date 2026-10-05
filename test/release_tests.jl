@testset "correctness regressions" begin
    # A modular collision must not turn a nonzero exact sum into zero.
    s = FactorialSum(0:1; factors=((1,1,120),), alternating=true)
    tabs = (QR.ClassicalModTable(UInt64(5),2), QR.ClassicalModTable(UInt64(13),2))
    @test QR._sum_vanishes_mod(s,tabs) && !QR.is_classical_zero(s,tabs)
    @test qeval(Exact(), s) == 1-big(2)^120

    # Weights change zero identities; tiny nonzero results must survive cancellation.
    s = FactorialSum(0:2; factors=((1,0,-1),(-1,2,-1)), alternating=true)
    @test !QR._analytic_sum_iszero(s,0.8,0)
    @test iszero(QR.analytic_value(s,0.8; weight=(1,0)))
    q = 1+big(1)//big(2)^500
    @test qeval(s; q) ≈ BigFloat(2/(q+inv(q))-1) rtol=big"1e-60"

    # An actual level cancellation, including writes across a packed-bit boundary.
    cancelled = (1,3//2,3//2,2,3//2,3//2)
    @test iszero(q6j(Exact(10), cancelled...))
    labels = [isodd(i) ? cancelled : ONES for i in 1:65]
    @test iszero_at(10,labels; threads=4) == isodd.(1:65)
    @test level_spectrum(cancelled...; k=[3,10],prove=true) == [:inadmissible,:cancels]
    @test level_spectrum(ONES...; k=Int[],threads=4) == Symbol[]
    @test_throws DomainError iszero_at(-1,ONES...)

    # Dependent radical classes need exact equality, even after tiny rescaling.
    a = q6j(Exact(4), 0,0,0,1//2,1//2,1//2)
    b = q6j(Exact(4), 0,0,0,3//2,3//2,3//2)
    tiny = big(1)//big(10)^100
    @test a == b && iszero(tiny*(a-b)) && !iszero(tiny*(a+b))

    # Retain the sign on the negative-real branch and precision in matrix output.
    l = (1//2,1//2,0,1//2,1//2,1)
    @test fsymbol(l...; q=-0.8) == fsymbol(l...; q=complex(-0.8))
    @test real(fsymbol(l...; q=-0.8)) < 0
    smatrix(7) # populate a Float64 table first
    setprecision(BigFloat,128) do
        F, es, fs = fmatrix(1,1,1,1; q=big"0.8",T=Float64)
        @test eltype(F) === BigFloat
        @test F ≈ [fsymbol(1,1,e,1,1,f; q=big"0.8") for e in es, f in fs] rtol=big"1e-30"
        S, = smatrix(7; T=BigFloat)
        @test eltype(S) === BigFloat && precision(S[1,1]) == 128
        @test S[1,1] ≈ sqrt(big(2)/9)*sinpi(big(1)/9) rtol=eps(BigFloat)*4
    end

    # Scaled intermediates must handle both ends of the Float64 exponent range.
    for q in (nextfloat(0.0),floatmax(Float64))
        expected = Float64(-inv(BigFloat(q)+inv(BigFloat(q))))
        @test q6j(1//2,1//2,0,1//2,1//2,0; q) ≈ expected rtol=1e-13 atol=nextfloat(0.0)
    end
end
