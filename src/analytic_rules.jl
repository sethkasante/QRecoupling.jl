# Generic-q factorial rules. Exponents stay separate throughout products and
# compensated summation; only the final result is converted to the output range.
struct AnalyticScaled{T}
    m::T
    e::Int
end

# ---------------------------------------------------------------------------
#  Double-word mantissas for real q
#
#  The machine-precision table is built by two length-N recurrences, so its
#  entries carry O(N u) relative error (measured: 15,000 u on the unit circle)
#  and no error bound over them can certify a Float64 answer. Carrying the
#  mantissa as a double word makes every operation accurate to u², so the same
#  bound certifies and the arbitrary-precision tier is only needed for genuinely
#  ill-conditioned q. Real q needs no transcendental functions here: q is exact,
#  q^{-2} is one division, and 1 − q^{-2} is exact as a double word even at
#  q ≈ 1, where the Float64 build needs `expm1`.
# ---------------------------------------------------------------------------
struct DWNum <: Real
    hi::Float64
    lo::Float64
end
DWNum(x::Real) = DWNum(Float64(x),0.0)
DWNum(x::DWNum) = x
DWNum(x::BigFloat) = (h = Float64(x); DWNum(h,Float64(x-h)))

# `Complex{DWNum}` is a legal Complex because DWNum <: Real, so products,
# quotients and reciprocals come from Base; only these two are missing.
Base.hypot(a::DWNum,b::DWNum) = sqrt(a*a+b*b)

# Written out rather than taken from Base, whose generic complex square root
# reaches for `nextfloat` and other parts of the AbstractFloat interface.
function Base.sqrt(z::Complex{DWNum})
    x,y = real(z), imag(z)
    if iszero(y)
        return signbit(x) ? Complex(zero(DWNum),flipsign(sqrt(-x),y)) :
                            Complex(sqrt(x),zero(DWNum))
    end
    r = hypot(x,y)
    two = DWNum(2.0)
    if !signbit(x)
        a = sqrt((r+x)/two)
        return Complex(a,y/(two*a))
    end
    b = flipsign(sqrt((r-x)/two),y)
    return Complex(y/(two*b),b)
end
_dwnum(q::Real) = DWNum(q)
_dwnum(q::Complex) = Complex{DWNum}(DWNum(real(q)),DWNum(imag(q)))
_narrow(m::DWNum) = Float64(m)
_narrow(m::Complex{DWNum}) = ComplexF64(Float64(real(m)),Float64(imag(m)))

Base.Float64(a::DWNum) = a.hi + a.lo
Base.convert(::Type{DWNum},x::Real) = DWNum(x)
Base.BigFloat(a::DWNum; kw...) = BigFloat(a.hi; kw...) + BigFloat(a.lo; kw...)
Base.promote_rule(::Type{DWNum},::Type{<:Union{Integer,Float16,Float32,Float64}}) = DWNum
# Wider than a double word: promote outwards, or `promote_type` recurses.
Base.promote_rule(::Type{DWNum},::Type{BigFloat}) = BigFloat
Base.zero(::Type{DWNum}) = DWNum(0.0,0.0)
Base.one(::Type{DWNum}) = DWNum(1.0,0.0)
Base.zero(::DWNum) = zero(DWNum)
Base.one(::DWNum) = one(DWNum)
Base.iszero(a::DWNum) = iszero(a.hi)
Base.isfinite(a::DWNum) = isfinite(a.hi) && isfinite(a.lo)
Base.:(-)(a::DWNum) = DWNum(-a.hi,-a.lo)
Base.abs(a::DWNum) = signbit(a.hi) ? -a : a
Base.signbit(a::DWNum) = signbit(a.hi)
Base.:(<)(a::DWNum,b::DWNum) = Float64(a-b) < 0
Base.:(==)(a::DWNum,b::DWNum) = a.hi == b.hi && a.lo == b.lo
Base.isequal(a::DWNum,b::DWNum) = isequal(a.hi,b.hi) && isequal(a.lo,b.lo)
Base.hash(a::DWNum,h::UInt) = hash(a.lo,hash(a.hi,h))
Base.ldexp(a::DWNum,e::Int) = DWNum(ldexp(a.hi,e),ldexp(a.lo,e))
Base.frexp(a::DWNum) = (fr = frexp(a.hi); (DWNum(fr[1],ldexp(a.lo,-fr[2])),fr[2]))
Base.precision(::Type{DWNum}) = 2 * precision(Float64)
Base.eps(::Type{DWNum}) = eps(Float64)^2 / 2
Base.eps(::DWNum) = eps(DWNum)
# Enough of the AbstractFloat interface for Base's generic complex sqrt.
Base.AbstractFloat(a::DWNum) = a
Base.float(a::DWNum) = a
Base.isnan(a::DWNum) = isnan(a.hi) || isnan(a.lo)
Base.isinf(a::DWNum) = isinf(a.hi)
Base.copysign(a::DWNum,b::Real) = signbit(b) ? -abs(a) : abs(a)
Base.flipsign(a::DWNum,b::Real) = signbit(b) ? -a : a
Base.show(io::IO,a::DWNum) = print(io,"DWNum(",Float64(a),")")

Base.:(+)(a::DWNum,b::DWNum) = DWNum(_dw_add(a.hi,a.lo,b.hi,b.lo)...)
Base.:(-)(a::DWNum,b::DWNum) = a + (-b)
Base.:(*)(a::DWNum,b::DWNum) = DWNum(_dw_mul(a.hi,a.lo,b.hi,b.lo)...)
Base.:(/)(a::DWNum,b::DWNum) = DWNum(_dw_div(a.hi,a.lo,b.hi,b.lo)...)
Base.sqrt(a::DWNum) = DWNum(_dw_sqrt(a.hi,a.lo)...)

"(ah + al) + (bh + bl) as a double word."
@inline function _dw_add(ah,al,bh,bl)
    s, e = _two_sum(ah,bh)
    t, f = _two_sum(al,bl)
    e += t
    s, e = _fast_two_sum(s,e)
    return _fast_two_sum(s,e+f)
end

"(ah + al) / (bh + bl) as a double word."
@inline function _dw_div(ah,al,bh,bl)
    q1 = ah / bh
    p, e = _dw_mul(bh,bl,q1,0.0)
    r = ((ah - p) - e) + al
    q2 = r / bh
    return _fast_two_sum(q1,q2)
end

_aldexp(x::Real,e::Int) = iszero(e) ? x : ldexp(x,e)
_aldexp(x::Complex,e::Int) = iszero(e) ? x : complex(ldexp(real(x),e),ldexp(imag(x),e))

# Scale by the larger component rather than the modulus: `abs` of a Complex is a
# hypot, and one of those per multiply dominates the machine-precision kernel.
# The mantissa then lives in [1/2, 2sqrt(2)) in modulus, which is all that the
# exponent split needs.
@inline _ascale_ref(m::Real) = abs(m)
@inline _ascale_ref(m::Complex) = max(abs(real(m)),abs(imag(m)))

function _ascaled(m::T,e::Int=0) where T
    iszero(m) && return AnalyticScaled(zero(m),0)
    _,k = frexp(_ascale_ref(m))
    return AnalyticScaled(_aldexp(m,-k),Base.checked_add(e,k))
end

# MPFR already owns a wide exponent; avoid copying and renormalizing every
# arbitrary-precision product. Machine-float kernels retain explicit scaling.
_ascaled(m::BigFloat,e::Int=0) = AnalyticScaled(_aldexp(m,e),0)
_ascaled(m::Complex{BigFloat},e::Int=0) = AnalyticScaled(_aldexp(m,e),0)
_amul(a::AnalyticScaled,b::AnalyticScaled) = _ascaled(a.m*b.m,Base.checked_add(a.e,b.e))

function _adiv(a::AnalyticScaled,b::AnalyticScaled)
    iszero(b.m) && throw(DomainError(b.m,"zero q-integer in analytic ratio; use a level target for roots of unity"))
    _ascaled(a.m/b.m,Base.checked_sub(a.e,b.e))
end
function _aadd(a::AnalyticScaled,b::AnalyticScaled)
    iszero(a.m) && return b
    iszero(b.m) && return a
    e=max(a.e,b.e)
    _ascaled(_aldexp(a.m,a.e-e)+_aldexp(b.m,b.e-e),e)
end
_aneg(a::AnalyticScaled) = AnalyticScaled(-a.m,a.e)
_asub(a::AnalyticScaled,b::AnalyticScaled) = _aadd(a,_aneg(b))
_aabs(a::AnalyticScaled) = AnalyticScaled(abs(a.m),a.e)
function _apow(a::AnalyticScaled,n::Integer)
    n==1 && return a
    n < 0 && return _apow(_adiv(_ascaled(one(a.m)),a),-Int(n))
    r=_ascaled(one(a.m)); p=a; k=Int(n)
    while k != 0
        isodd(k) && (r=_amul(r,p))
        k >>= 1
        k != 0 && (p=_amul(p,p))
    end
    r
end

@inline function _amul_power(v,a,n)
    n==1 && return _amul(v,a)
    n==-1 && return _adiv(v,a)
    n==0 && return v
    return _amul(v,_apow(a,n))
end

function _asqrt(a::AnalyticScaled)
    e=a.e
    _ascaled(sqrt(isodd(e) ? 2a.m : a.m),fld(e,2))
end
function _aexp(z)
    # real(z) is allowed to exceed the floating-point exponential range.
    l2=log(oftype(real(z),2))
    e=floor(Int,real(z)/l2)
    _ascaled(exp(z-e*l2),e)
end

"q^r for a rational r, as a scaled value."
_aexp_pow(q,r) = _aexp(r*log(q))

# A double word has no logarithm or exponential, and it needs neither: the
# radical's exponent is a half-integer, so q^r is an integer power of q or of
# its principal square root. Same branch as exp(r log q), since raising a fixed
# complex number to an integer power cannot cross the cut.
function _aexp_pow(q::Complex{DWNum},r::Rational)
    n = numerator(r); d = denominator(r)
    d == 1 && return _apow(_ascaled(q),n)
    d == 2 && return _apow(_asqrt(_ascaled(q)),n)
    throw(ArgumentError("radical phase $r is neither integral nor half-integral"))
end
_avalue(a::AnalyticScaled) = _aldexp(a.m,a.e)

# `inverses` and `balanced` are filled on first use: the first is read only by
# unit-slope negative powers, the second only by a square-root prefactor off the
# positive real axis, and each costs a full pass over the table to build.
mutable struct AnalyticRuleTable{T}
    q::T
    bits::Int
    ints::Vector{AnalyticScaled{T}}
    inverses::Vector{AnalyticScaled{T}}
    facts::Vector{AnalyticScaled{T}}  # index n+1 represents [n]!
    balanced::Vector{AnalyticScaled{T}}
end

function _table_inverses(tab::AnalyticRuleTable)
    length(tab.inverses) == length(tab.ints) && return tab.inverses
    unit = _ascaled(one(tab.q))
    tab.inverses = [_adiv(unit,v) for v in tab.ints]
    return tab.inverses
end

# [n] = product_{d|n,d>=2} Psi_d. These balanced factors retain cyclotomic
# square extraction without constructing any DCR ratios.
function _table_balanced(tab::AnalyticRuleTable)
    length(tab.balanced) == length(tab.ints) && return tab.balanced
    N = length(tab.ints)
    b = copy(tab.ints)
    for d in 2:N÷2
        for n in 2d:d:N
            b[n] = _adiv(b[n],b[d])
        end
    end
    tab.balanced = b
    return b
end

function _analytic_table(q::T,N::Int) where T
    R=typeof(real(q)); bits=precision(R)
    ints=Vector{AnalyticScaled{T}}(undef,N)
    facts=Vector{AnalyticScaled{T}}(undef,N+1)
    facts[1]=_ascaled(one(q))
    # Choose a logarithm with nonnegative real part, using reciprocal invariance
    # of q-integers only (never of the complex radical branch).
    t=q isa Real ? log(abs(q)) : log(q)
    real(t)<0 && (t=-t)
    den=-expm1(-2t)
    decay=exp(-2t)
    numerator=den
    growth=_aexp(t)
    power=_ascaled(one(q))
    for n in 1:N
        v = n==1 ? _ascaled(one(q)) :
            _amul(power,_ascaled(numerator/den))
        q isa Real && q<0 && iseven(n) && (v=_aneg(v))
        ints[n]=v
        facts[n+1]=_amul(facts[n],v)
        power=_amul(power,growth)
        numerator=numerator*decay+den
    end
    empty=AnalyticScaled{T}[]
    AnalyticRuleTable(q,bits,ints,empty,facts,copy(empty))
end

"Double-word table, built from exact arithmetic only (no log/exp/expm1)."
function _analytic_table(q::T,N::Int) where {T<:Union{DWNum,Complex{DWNum}}}
    ints=Vector{AnalyticScaled{T}}(undef,N)
    facts=Vector{AnalyticScaled{T}}(undef,N+1)
    unit=one(T)
    facts[1]=_ascaled(unit)
    # Reciprocal invariance of the q-integers, so that q^{-2} cannot overflow.
    # For real q the sign is carried separately, as [n]_{-q} = (-1)^{n+1} [n]_q.
    aq=q isa Real ? abs(q) : q
    growth=_ascaled(abs(aq) < one(DWNum) ? unit/aq : aq)
    decay=_avalue(_adiv(_ascaled(unit),_amul(growth,growth)))
    # |growth| >= 1 keeps |decay| <= 1; the subtraction is exact as a double
    # word even when q is within an ulp of 1.
    den=unit-decay
    numer=den
    power=_ascaled(unit)
    negative=q isa Real && signbit(q)
    for n in 1:N
        v = n==1 ? _ascaled(unit) : _amul(power,_ascaled(numer/den))
        negative && iseven(n) && (v=_aneg(v))
        ints[n]=v
        facts[n+1]=_amul(facts[n],v)
        power=_amul(power,growth)
        numer=numer*decay+den
    end
    empty=AnalyticScaled{T}[]
    AnalyticRuleTable(q,precision(DWNum),ints,empty,facts,copy(empty))
end

function _analytic_table(q::T,N::Int,w::EvaluationWorkspace) where T
    key=(T,precision(typeof(real(q))))
    tab=get(w.analytic,key,nothing)
    if tab isa AnalyticRuleTable{T} && isequal(tab.q,q) && length(tab.ints)>=N
        return tab
    end
    # Bounded caller-owned storage: one q and at most four precision/type tiers.
    if length(w.analytic)>=4 || any(t->!isequal(t.q,q),values(w.analytic))
        empty!(w.analytic)
    end
    N=max(N,maximum(t->length(t.ints),values(w.analytic);init=0))
    tab=_analytic_table(q,N)
    w.analytic[key]=tab
    tab
end
_analytic_table(q,N,::Nothing) = _analytic_table(q,N)

function _analytic_prefactor(s,tab)
    q=tab.q
    if !s.sqrt_pre || (q isa Real && q>0)
        v=_ascaled(one(q))
        for (n,c) in s.pre
            v=_amul_power(v,tab.facts[n+1],Int(c))
        end
        return s.sqrt_pre ? _asqrt(v) : v
    end
    # Preserve snapshot_square_root's extracted sign and balanced radical branch.
    b=CycloBuffer(length(tab.ints))
    for (n,c) in s.pre
        add_qfact!(b,Int(n),Int(c))
    end
    root,rad=snapshot_square_root(b)
    psi=_table_balanced(tab)
    function balanced(m)
        v=_ascaled(oftype(q,m.sign))
        for (d,e) in m.phi_exps
            v=_amul_power(v,psi[d],e)
        end
        v
    end
    pr=balanced_phase(root); pd=balanced_phase(rad)
    r=balanced(root); v=balanced(rad)
    if q isa Real
        r=_amul(r,_apow(_ascaled(q),pr))
        v=_amul(v,_apow(_ascaled(q),pd))
        return _amul(r,_asqrt(v))
    end
    if abs(imag(v.m)) <= sqrt(eps(one(real(q))))*abs(v.m)
        v=AnalyticScaled(oftype(q,real(v.m)),v.e)
    end
    return _amul(_amul(r,_asqrt(v)),_aexp_pow(q,pr+pd//2))
end

"Mantissas that carry guard digits of their own, so the sum needs no compensation."
_wide_arith(q) = real(q) isa BigFloat || real(q) isa DWNum

function _analytic_pass(s::FactorialSum,tab::AnalyticRuleTable)
    q=tab.q
    ints=tab.ints; facts=tab.facts
    # One unit-slope dividing step anywhere in the rule pays for the table of
    # reciprocals; without one the table is never built. `_factor_step` negates
    # the power for a decreasing factor, so the step divides when a*c < 0.
    inverses=any(f->abs(Int(f.a))==1 && Int(f.a)*Int(f.c)<0,s.fac) ? _table_inverses(tab) : ints
    t=_ascaled(one(q))
    for f in s.fac
        t=_amul_power(t,facts[_arg(f,s.zlo)+1],Int(f.c))
    end
    s.alternating && isodd(s.zlo) && (t=_aneg(t))
    acc=t; comp=_ascaled(zero(q)); mass=_aabs(t)
    for z in s.zlo:s.zhi-1
        ratio=_ascaled(s.alternating ? -one(q) : one(q))
        for f in s.fac
            lo,hi,c=_factor_step(f,z)
            c==0 && continue
            if lo==hi
                ratio=_amul_power(ratio,c>0 ? ints[hi] : inverses[hi],abs(c))
                continue
            end
            # Prefix division makes arbitrary slopes independent of interval length.
            block=lo>hi ? _ascaled(one(q)) : _adiv(facts[hi+1],facts[lo])
            ratio=_amul_power(ratio,block,c)
        end
        t=_amul(t,ratio)
        if _wide_arith(q)
            # Double-word and arbitrary-precision passes already carry guard
            # bits; compensation is needed only in the machine-precision pass.
            acc=_aadd(acc,t)
        else
            y=_asub(t,comp); next=_aadd(acc,y)
            comp=_asub(_asub(next,acc),y); acc=next
        end
        mass=_aadd(mass,_aabs(t))
    end
    v=_amul(_analytic_prefactor(s,tab),acc)
    s.sign0<0 && (v=_aneg(v))
    cond=iszero(acc.m) ? oftype(real(q),Inf) : _avalue(_adiv(mass,_aabs(acc)))
    return v,cond
end

function _analytic_close(a,b,tol)
    iszero(a.m) && iszero(b.m) && return true
    iszero(b.m) && return false
    d=_asub(a,b)
    iszero(d.m) || _avalue(_adiv(_aabs(d),_aabs(b))) <= tol
end

function _analytic_q(q::Number)
    qq=float(q*1.0)
    isfinite(qq) && !iszero(qq) || throw(DomainError(q,"analytic q must be finite and nonzero"))
    # Do not perturb exactly representable singular roots through log/exp.
    # Their cancellations belong to the level/symbolic projection machinery.
    qq in (-1,im,-im) && throw(DomainError(q,"use a level target at roots of unity"))
    return qq
end

"Direct, scaled analytic evaluation with compensated sums and precision escalation."
function analytic_value(s::FactorialSum,q::Number;workspace=nothing)
    qq=_analytic_q(q); T=typeof(qq); R=typeof(real(qq))
    is_empty_sum(s) && return zero(T)
    N=max_argument(s)
    # Conservative operation count: table and prefix-product roundings as well as
    # the accumulated ratio operations. Multiplied by the working precision's
    # unit roundoff and by the sum's own cancellation it bounds the result.
    ops=16*(N+1)*(1+sum(f->abs(Int(f.c)),s.fac;init=0)+
                   sum(p->abs(Int(p.second)),s.pre;init=0))*(s.zhi-s.zlo+1)
    tol=R===BigFloat ? eps(one(real(qq)))*32 : R(64)*eps(R)
    if R===Float64
        # Double-word tier: u² per operation, so this bound is met for every
        # well-conditioned q and the arbitrary-precision tiers below are reached
        # only near a singularity of the rule.
        vd,kd=_analytic_pass(s,_analytic_table(_dwnum(qq),N,workspace))
        if isfinite(kd) && ops*eps(DWNum)*kd <= tol
            return T(_aldexp(_narrow(vd.m),vd.e))
        end
        v=AnalyticScaled(T(_narrow(vd.m)),vd.e)
    else
        v,kappa=_analytic_pass(s,_analytic_table(qq,N,workspace))
        if isfinite(kappa) && ops*eps(one(real(qq)))*kappa <= tol
            return T(_avalue(v))
        end
    end
    # Rebuild from the supplied q at each precision, never from rounded table
    # entries. A tier that certifies its own bound is accepted on its own: the
    # machine-precision value is not a reliable witness, so requiring the two to
    # agree only forced a second arbitrary-precision pass. An exactly cancelling
    # sum carries no bound, and there two tiers must agree instead.
    previous=v
    bits=max(128,precision(R)+32)
    for _ in 1:8
        next,condition=setprecision(BigFloat,bits) do
            qb=q isa Real ? BigFloat(q) : Complex{BigFloat}(q)
            _analytic_pass(s,_analytic_table(qb,N,workspace))
        end
        if iszero(next.m)
            _analytic_close(previous,next,tol/4) && return T(_avalue(next))
        elseif isfinite(condition) && ops*ldexp(one(condition),-bits)*condition < tol/4
            return T(_avalue(next))
        end
        previous=next
        bits*=2
    end
    throw(ErrorException("analytic evaluation did not converge; increase input precision or use an exact level target"))
end
