# Generic-q factorial rules. Exponents stay separate throughout products and
# compensated summation; only the final result is converted to the output range.
struct AnalyticScaled{T}
    m::T
    e::Int
end

_aldexp(x::Real,e::Int) = iszero(e) ? x : ldexp(x,e)
_aldexp(x::Complex,e::Int) = iszero(e) ? x : complex(ldexp(real(x),e),ldexp(imag(x),e))

function _ascaled(m::T,e::Int=0) where T
    iszero(m) && return AnalyticScaled(zero(m),0)
    _,k = frexp(abs(m))
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
_avalue(a::AnalyticScaled) = _aldexp(a.m,a.e)

struct AnalyticRuleTable{T}
    q::T
    bits::Int
    ints::Vector{AnalyticScaled{T}}
    inverses::Vector{AnalyticScaled{T}}
    facts::Vector{AnalyticScaled{T}}  # index n+1 represents [n]!
    balanced::Vector{AnalyticScaled{T}}
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
    # [n] = product_{d|n,d>=2} Ψ_d. These balanced factors retain
    # cyclotomic square extraction without constructing any DCR ratios.
    balanced = q isa Real && q>0 ? AnalyticScaled{T}[] : copy(ints)
    if !isempty(balanced)
        for d in 2:N÷2
            for n in 2d:d:N
                balanced[n]=_adiv(balanced[n],balanced[d])
            end
        end
    end
    unit=_ascaled(one(q))
    inverses=[_adiv(unit,v) for v in ints]
    AnalyticRuleTable(q,bits,ints,inverses,facts,balanced)
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
    function balanced(m)
        v=_ascaled(oftype(q,m.sign))
        for (d,e) in m.phi_exps
            v=_amul_power(v,tab.balanced[d],e)
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
    return _amul(_amul(r,_asqrt(v)),_aexp((pr+pd//2)*log(q)))
end

function _analytic_pass(s::FactorialSum,tab::AnalyticRuleTable)
    q=tab.q
    t=_ascaled(one(q))
    for f in s.fac
        t=_amul_power(t,tab.facts[_arg(f,s.zlo)+1],Int(f.c))
    end
    s.alternating && isodd(s.zlo) && (t=_aneg(t))
    acc=t; comp=_ascaled(zero(q)); mass=_aabs(t)
    for z in s.zlo:s.zhi-1
        ratio=_ascaled(s.alternating ? -one(q) : one(q))
        for f in s.fac
            lo,hi,c=_factor_step(f,z)
            c==0 && continue
            if lo==hi
                ratio=_amul_power(ratio,c>0 ? tab.ints[hi] : tab.inverses[hi],abs(c))
                continue
            end
            # Prefix division makes arbitrary slopes independent of interval length.
            block=lo>hi ? _ascaled(one(q)) : _adiv(tab.facts[hi+1],tab.facts[lo])
            ratio=_amul_power(ratio,block,c)
        end
        t=_amul(t,ratio)
        if real(q) isa BigFloat
            # Escalated passes already carry guard bits; compensation is needed
            # only in the machine-precision pass.
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
    tab=_analytic_table(qq,N,workspace)
    v,kappa=_analytic_pass(s,tab)
    # Conservative diagnostic, not an interval certificate. Includes table and
    # prefix-product errors as well as accumulated ratio operations.
    ops=16*(N+1)*(1+sum(f->abs(Int(f.c)),s.fac;init=0)+
                   sum(p->abs(Int(p.second)),s.pre;init=0))*(s.zhi-s.zlo+1)
    tol=R===BigFloat ? eps(one(real(qq)))*32 : R(64)*eps(R)
    if isfinite(kappa) && ops*eps(one(real(qq)))*kappa <= tol
        return T(_avalue(v))
    end
    # Rebuild from the supplied q at each precision, never from rounded table
    # entries. Agreement is checked in scaled form, before output underflow.
    previous=v
    bits=max(128,precision(R)+32)
    for _ in 1:8
        next,condition=setprecision(BigFloat,bits) do
            qb=q isa Real ? BigFloat(q) : Complex{BigFloat}(q)
            _analytic_pass(s,_analytic_table(qb,N,workspace))
        end
        if _analytic_close(previous,next,tol/4) &&
                (iszero(next.m) || isfinite(condition) && ops*ldexp(one(condition),-bits)*condition < tol/4)
            return T(_avalue(next))
        end
        previous=next
        bits*=2
    end
    throw(ErrorException("analytic evaluation did not converge; increase input precision or use an exact level target"))
end
