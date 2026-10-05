# Generic-q factorial rules. Exponents stay separate throughout products and
# compensated summation; only the final result is converted to the output range.
struct AnalyticScaled{T}
    m::T
    e::Int
end

# ---------------------------------------------------------------------------
#  Double-word mantissas for real q
#
#  Machine-precision tables carry O(N u) error, which no bound over them can certify; double-word mantissas
#  are accurate to u², so the same bound certifies and arbitrary precision is left for ill-conditioned q.
#  Real q needs no transcendental functions (1 − q⁻² is exact as a double word even at q ≈ 1).
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
    circle::Bool                      # |q| = 1 to within the input's own rounding; see `_on_unit_circle`
    ints::Vector{AnalyticScaled{T}}
    inverses::Vector{AnalyticScaled{T}}
    facts::Vector{AnalyticScaled{T}}  # index n+1 represents [n]!
    balanced::Vector{AnalyticScaled{T}}
    roots::Vector{AnalyticScaled{T}}  # √Ψ_d on the package's branch; complex q only
    condmax::Vector{Float64}          # prefix maxima of the q-integer condition; see `_prefactor_cond`
    accurate::Bool                    # entries rounded from double words: each within half an ulp
end
AnalyticRuleTable(q,bits,circle,ints,inverses,facts,balanced,roots,condmax) =
    AnalyticRuleTable(q,bits,circle,ints,inverses,facts,balanced,roots,condmax,false)

"""
    _on_unit_circle(q) -> Bool

Whether `q` is on the unit circle up to the rounding of its own Float64 representation (`cispi(t)` gives
`|q| = 1 ± 2⁻⁵³`); the same tolerance at every precision keeps the tiers in agreement. On the circle each
`Ψ_d` is real and a negative one sits on the branch cut, so `_analytic_prefactor` drops the rounding-noise
imaginary part and fixes the branch at `+i√|Ψ_d|` (otherwise 21 of 72 symbols flipped sign under a one-ulp
change of `q`).
"""
_on_unit_circle(::Real) = false
@inline function _on_unit_circle(q::Complex)
    a = Float64(real(q) * real(q) + imag(q) * imag(q))
    return abs(a - 1.0) <= 8 * eps(Float64)
end

"Drop an imaginary part that is known to be rounding noise, so that `sqrt` lands on a definite branch."
@inline _deimag(a::AnalyticScaled) = a
@inline _deimag(a::AnalyticScaled{Complex{R}}) where {R} =
    AnalyticScaled(Complex{R}(real(a.m), zero(R)), a.e)

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

"""
`√Ψ_d` for every `d`, one principal root per balanced factor (real part only on the unit circle). Cached:
it is the only transcendental in a complex-q prefactor and does not depend on the rule.
"""
function _table_roots(tab::AnalyticRuleTable)
    length(tab.roots) == length(tab.ints) && return tab.roots
    psi = _table_balanced(tab)
    r = tab.circle ? [_asqrt(_deimag(p)) for p in psi] : [_asqrt(p) for p in psi]
    tab.roots = r
    return r
end

function _analytic_table(q::T,N::Int) where T
    q isa Complex && _negative_real_axis(q) && return _negative_axis_table(q,N)
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
    AnalyticRuleTable(q,bits,_on_unit_circle(q),ints,empty,facts,copy(empty),copy(empty),Float64[])
end

"Double-word table, built from exact arithmetic only (no log/exp/expm1)."
function _analytic_table(q::T,N::Int) where {T<:Union{DWNum,Complex{DWNum}}}
    q isa Complex && _negative_real_axis(q) && return _negative_axis_table(q,N)
    ints=Vector{AnalyticScaled{T}}(undef,N)
    facts=Vector{AnalyticScaled{T}}(undef,N+1)
    unit=one(T)
    facts[1]=_ascaled(unit)
    # Reciprocal invariance of the q-integers, so that q^{-2} cannot overflow.
    # For real q the sign is carried separately, as [n]_{-q} = (-1)^{n+1} [n]_q.
    aq=q isa Real ? abs(q) : q
    aqs=_ascaled(aq)
    # Base's Float64 hypot scales its operands; the double-word hypot can
    # underflow/overflow before this reciprocal-selection test is made.
    growth=abs(_narrow(aq)) < 1 ? _adiv(_ascaled(unit),aqs) : aqs
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
    AnalyticRuleTable(q,precision(DWNum),_on_unit_circle(q),ints,empty,facts,copy(empty),copy(empty),Float64[])
end

"""
Real Float64 `q`: the double-word table, rounded once per entry, so every entry is within half an ulp as the
plain tier's bound assumes. The Float64 recurrence drifted by up to 88,000 u in `[n]!` near q ≈ 1, and values
up to 2.3e−12 were accepted against the 9.1e−13 promise. Cached per `q`.
"""
function _analytic_table(q::Float64,N::Int)
    # The short recurrence below rescales once per step. Very large or small q
    # can exhaust that headroom; use fully scaled products in that case.
    if !(0x1p-256 <= abs(q) <= 0x1p256)
        d = _analytic_table(DWNum(q),N)
        narrow(a) = AnalyticScaled(Float64(a.m),a.e)
        empty = AnalyticScaled{Float64}[]
        return AnalyticRuleTable(q,precision(Float64),false,map(narrow,d.ints),empty,
                                 map(narrow,d.facts),copy(empty),copy(empty),Float64[],true)
    end
    # [n+1] = x[n] − [n−1], x = g + 1/g with g = max(|q|, 1/|q|), in double words: near q = 1 the recurrence
    # is only weakly dominant and its error grows like n²·u², which at u² is still far below an ulp. Values
    # carry a separate exponent so that |q|^n cannot overflow; each entry is rounded to Float64 once.
    ints=Vector{AnalyticScaled{Float64}}(undef,N)
    facts=Vector{AnalyticScaled{Float64}}(undef,N+1)
    facts[1]=_ascaled(1.0)
    aq=abs(q)
    g=aq<1 ? inv(DWNum(aq)) : DWNum(aq)
    x=g+inv(g)
    prev=zero(DWNum); cur=one(DWNum); e=0          # [n] = cur·2^e
    fm=one(DWNum); fe=0                             # [n]! = fm·2^fe
    negative=signbit(q)
    for n in 1:N
        v=_ascaled(Float64(cur),e)
        neg=negative && iseven(n)                   # [n]_{−q} = (−1)^{n+1} [n]_q
        ints[n]=neg ? _aneg(v) : v
        fm=fm*cur; fe+=e
        neg && (fm=-fm)
        facts[n+1]=_ascaled(Float64(fm),fe)
        if !(0x1p-400<abs(fm.hi)<0x1p400)                 # keep the running product in range
            ex=exponent(fm.hi); fm=ldexp(fm,-ex); fe+=ex
        end
        nxt=x*cur-prev
        prev,cur=cur,nxt
        if abs(cur.hi)>0x1p500
            prev=ldexp(prev,-500); cur=ldexp(cur,-500); e+=500
        end
    end
    empty=AnalyticScaled{Float64}[]
    AnalyticRuleTable(q,precision(Float64),false,ints,empty,facts,copy(empty),copy(empty),Float64[],true)
end

"""
Complex Float64 `q` near the unit circle: the double-word table, rounded once per entry. Accurate entries need
no `pcond` charge, so the plain pass is accepted on the circle (2–4 µs instead of 15–22 µs in double words).
The value is still the one at the `q` passed. Off the circle the recurrence is kept.
"""
function _analytic_table(q::ComplexF64,N::Int)
    _negative_real_axis(q) && return _negative_axis_table(q,N)
    abs(abs(q)-1) < 0.1 || return invoke(_analytic_table,Tuple{Any,Int},q,N)
    d=_analytic_table(_dwnum(q),N)
    r(a)=_ascaled(_narrow(a.m),a.e)
    empty=AnalyticScaled{ComplexF64}[]
    AnalyticRuleTable(q,precision(Float64),_on_unit_circle(q),map(r,d.ints),empty,map(r,d.facts),
                      copy(empty),copy(empty),Float64[],true)
end

"Negative-axis q-integers are real: lift the real table without logarithmic phase noise."
function _negative_axis_table(q::Complex,N::Int)
    t = _analytic_table(real(q),N)
    lift(a) = AnalyticScaled(complex(a.m),a.e)
    empty = AnalyticScaled{typeof(q)}[]
    return AnalyticRuleTable(q,t.bits,false,map(lift,t.ints),empty,map(lift,t.facts),
                             copy(empty),copy(empty),Float64[],t.accurate)
end

@inline _negative_real_axis(q) = real(q) < 0 && iszero(imag(q))

"Multiply a real scaled value by i^n using exact sign changes and component swaps."
@inline function _quarter_turn(a::AnalyticScaled{R},n::Int) where {R<:Real}
    r = mod(n,4); z = zero(R); m = a.m
    v = r == 0 ? complex(m,z) : r == 1 ? complex(z,m) :
        r == 2 ? complex(-m,z) : complex(z,-m)
    return AnalyticScaled(v,a.e)
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
"""
Process-wide analytic tables for callers without a workspace, keyed by `q`, type and precision; a hit reuses
any table large enough. Shared between tasks: the lazy fields are built locally and assigned whole, so a racing
task rebuilds the same content.
"""
const _ANALYTIC_TABLE_CACHE = LRU{UInt,Any}(maxsize = 12)
const _ANALYTIC_TABLE_LOCK = ReentrantLock()

clear_analytic_caches!() = (@lock _ANALYTIC_TABLE_LOCK empty!(_ANALYTIC_TABLE_CACHE); nothing)

function _analytic_table_cached(q::T,N::Int) where T
    # A hashed key, not the `(type, precision, q)` tuple it stands for: that tuple holds a `DataType` and
    # so is heap-allocated, which is the one allocation a real-q symbol had left. A collision cannot give
    # a wrong answer because the hit is verified against `q` and the type before it is used.
    key = hash(q,hash(precision(typeof(real(q))),objectid(T)))
    tab = @lock _ANALYTIC_TABLE_LOCK get(_ANALYTIC_TABLE_CACHE,key,nothing)
    if tab isa AnalyticRuleTable{T} && isequal(tab.q,q) && length(tab.ints) >= N
        return tab
    end
    # Grow rather than shrink: a sweep that asks for a short table after a long one keeps the long one.
    want = (tab isa AnalyticRuleTable{T} && isequal(tab.q,q)) ? max(N,length(tab.ints)) : N
    fresh = _analytic_table(q,max(want,16))
    @lock _ANALYTIC_TABLE_LOCK (_ANALYTIC_TABLE_CACHE[key] = fresh)
    return fresh
end

_analytic_table(q,N,::Nothing) = _analytic_table_cached(q,N)

"""
`Π Ψ_d^{e_d}` for one balanced monomial, from the table's balanced factors (a function, not a closure: the
closure boxed a variable and allocated 35 times per prefactor).
"""
@inline function _balanced_mono(m, psi, q)
    acc=_ascaled(oftype(q,m.sign))
    for (d,e) in m.phi_exps
        acc=_amul_power(acc,psi[d],e)
    end
    return acc
end

"Largest q-factorial argument in a rule's prefactor: how far up the table the prefactor actually reads."
@inline function _pre_max(s::FactorialSum)
    m = 0
    for (n,c) in s.pre
        c == 0 && continue
        nn = Int(n); nn > m && (m = nn)
    end
    return m
end

"""
    _psi_exponents(s) -> (E, Dmax, qpow)

The Ψ-basis exponents of a rule's prefactor `Π [n]!^c`, and its residual q-power: `[n]! = Π_{d≥2} Ψ_d^{⌊n/d⌋}`
with `Ψ_d = q^{-φ(d)}Φ_d(q²)`, and `[n]` adds `q^{1−n}`. One flat array, one allocation.
"""
function _psi_exponents(s::FactorialSum)
    Dmax=0; qpow=0
    for (n,c) in s.pre
        nn=Int(n); nn<=1 && continue
        qpow+=Int(c)*(nn*(1-nn)÷2)
        nn>Dmax && (Dmax=nn)
    end
    E=zeros(Int,Dmax)
    for (n,c) in s.pre
        nn=Int(n); nn<=1 && continue
        cc=Int(c)
        @inbounds for d in 2:nn
            E[d]+=cc*(nn÷d)
        end
    end
    return E,Dmax,qpow
end

"""
    _analytic_prefactor(s, tab) -> (value, phase_weight, nroots)

The prefactor, with what the running bound needs: `phase_weight = |pr + pd/2|`, the exponent of the `q`
power evaluated as `exp(phase·log q)` (error `|phase·log q|·u`), and `nroots`, the roots and balanced-factor
products the caller's operation count does not see.
"""
function _analytic_prefactor(s,tab)
    q=tab.q
    if s.sqrt_pre && q isa Complex && _negative_real_axis(q)
        # Ψ₂=q+1/q is the only negative balanced factor on this axis; all
        # Ψ_d, d>2, are positive. Thus the factor-by-factor root is
        # i^E₂ sqrt(|Π[n]!^c|), E₂=Σ c floor(n/2). No exponent vector,
        # transcendental phase or branch-cut precision escalation is needed.
        p = _ascaled(one(real(q))); e2 = 0
        for (n,c) in s.pre
            a = tab.facts[n+1]
            p = _amul_power(p,AnalyticScaled(abs(real(a.m)),a.e),Int(c))
            e2 += Int(c) * (Int(n) ÷ 2)
        end
        return _quarter_turn(_asqrt(p),e2),0.0,0,1.0
    end
    if !s.sqrt_pre || (q isa Real && q>0)
        p=_ascaled(one(q))
        for (n,c) in s.pre
            p=_amul_power(p,tab.facts[n+1],Int(c))
        end
        return (s.sqrt_pre ? _asqrt(p) : p),0.0,0,1.0
    end
    psi=_table_balanced(tab)
    E,Dmax,qpow=_psi_exponents(s)
    aq,bq=split_exp(qpow)                     # the residual q-power, split the same way
    pr=aq; pd=bq
    r=_ascaled(one(q))                        # the prefactor's sign is +1: q-factorials carry none
    nroots=0
    if q isa Real
        # (the real balanced branch; `cut` is irrelevant there and reported as 1)
        v=_ascaled(one(q))
        @inbounds for d in 2:Dmax
            e=E[d]; e==0 && continue
            a,bb=split_exp(e)
            if a!=0; r=_amul_power(r,psi[d],a); pr+=a*_totient(d); nroots+=1; end
            if bb!=0; v=_amul_power(v,psi[d],bb); pd+=bb*_totient(d); nroots+=1; end
        end
        r=_amul(r,_apow(_ascaled(q),pr))
        v=_amul(v,_apow(_ascaled(q),pd))
        # Integer powers, so the phase costs a count of multiplications, not a logarithm.
        return _amul(r,_asqrt(v)),0.0,nroots+abs(Int(pr))+abs(Int(pd)),1.0
    end
    # Complex q: root the radical factor by factor, never the assembled product. √ is not multiplicative
    # across its cut; one root per Ψ_d depends only on the set of factors, so radicals cancel exactly in
    # coherence identities (Biedenharn–Elliott 40/40, against 9–36/40 for a root of the product).
    # `rad` is square-free, so every exponent is ±1.
    sv=_ascaled(one(q))
    sq=_table_roots(tab)
    # Off the circle, a small imaginary part of Ψ_d decides its root, so the plain tier is refused when the
    # factor is too close to the cut (`CUT_MARGIN`). On the circle Ψ_d is real: drop the noise and fix the
    # branch at +i√|Ψ_d|, which is continuous in arg q, matches `Level(k)` as θ → π/h, and still cancels in
    # coherence identities.
    cut=1.0
    oncut=tab.circle
    @inbounds for d in 2:Dmax
        e=E[d]; e==0 && continue
        a,bb=split_exp(e)
        if a!=0; r=_amul_power(r,psi[d],a); pr+=a*_totient(d); nroots+=1; end
        bb==0 && continue
        pd+=bb*_totient(d); nroots+=1
        if !oncut
            z=psi[d].m
            if real(z)<0 && !iszero(z)
                cut=min(cut,Float64(abs(imag(z))/abs(z)))
            end
        end
        sv=_amul_power(sv,sq[d],bb)
    end
    ph=pr+pd//2
    return _amul(_amul(r,sv),_aexp_pow(q,ph)),abs(Float64(ph)),nroots,cut
end

"""
Weight of the phase term `exp(P·log q)` in the plain tier's running bound: three roundings reach the exponent
plus one for the scaling; validated numerically.
"""
const CPHASE = 4

"""
    _prefactor_cond(q, N) -> Float64

How accurate a single table entry is, in units of `u`: `max_n |q^{-2n}| / |1 − q^{-2n}|`, the digits lost in
building `[n]` from `1 − q^{-2n}`, which blows up near a root of unity. Real `q` is excluded (it would reject
near-classical `q` for no gain). The loop stops once `|q^{-2n}|` leaves a neighbourhood of 1.
"""
function _prefactor_cond(q,N::Int)
    (q isa Real || _negative_real_axis(q)) && return 1.0
    u=inv(q*q); t=one(u); c=1.0
    for _ in 1:N
        t*=u
        at=abs(t); dn=abs(one(t)-t)
        iszero(dn) && return Inf
        c=max(c,Float64(at/dn))
        (at>4 || at<0.25) && break
    end
    return c
end

"""
    _prefactor_cond(tab, N) -> Float64

The same number, from prefix maxima cached on the table (on the circle the loop otherwise runs its full length
on every call).
"""
function _prefactor_cond(tab::AnalyticRuleTable,N::Int)
    (tab.q isa Real || _negative_real_axis(tab.q)) && return 1.0
    N <= 0 && return 1.0
    cm = tab.condmax
    if isempty(cm)
        M = length(tab.ints)
        cm = Vector{Float64}(undef,M)
        q = tab.q; u = inv(q*q); t = one(u); c = 1.0; done = false
        for n in 1:M
            if !done
                t *= u
                at = abs(t); dn = abs(one(t)-t)
                c = iszero(dn) ? Inf : max(c,Float64(at/dn))
                (at > 4 || at < 0.25 || !isfinite(c)) && (done = true)
            end
            cm[n] = c
        end
        tab.condmax = cm
    end
    return @inbounds cm[min(N,length(cm))]
end

"""
    NEAR_ROOT_COND

Above this prefactor condition number, `q` is within rounding of a root of unity of order at most `2N`. The
value returned is still the one at the Float64 `q` passed, but a one-ulp change of `q` changes it by a factor
of order one (`q6j(1,1,1,1,1,1; q = cispi(1/3))` is `5.8e15 + 1.0e16im`; at `k = 1` it is 0). Deliberate
nearby probes such as `cispi(1/7)·(1 + 1e−10)` are far below the threshold and certified normally.
"""
const NEAR_ROOT_COND = ldexp(1.0, -8) / eps(Float64)

"""
Name the level, if `q` is a rounded `exp(iπ/h)`, and say what was returned. Reached only from the branch
that has already computed the condition number, so it costs nothing when it does not fire.
"""
function _warn_near_root(q, N::Int)
    m = 0
    u = q * q; t = one(u)                       # q^{2n}; a root of unity of order dividing 2n
    for n in 1:N
        t *= u
        if abs(one(t) - t) <= 1e-12
            m = 2n
            break
        end
    end
    m == 0 && return nothing
    a = round(Int, Float64(angle(complex(q))) * m / (2 * pi))
    target = (a == 1 && iseven(m) && m >= 4) ? "`Level($(m ÷ 2 - 2))` or `Exact($(m ÷ 2 - 2))`" :
                                               "a level target"
    @warn("q is within rounding of a root of unity of order $m; the value returned is the one at the " *
          "Float64 q that was passed, and a one-ulp change in q changes it by a factor of order one. " *
          "Use $target for an exact root of unity.", maxlog = 3)
    return nothing
end

"""
How far a rooted factor must sit from the branch cut, relative to its size, before a plain Float64 prefactor is
trusted. A guard for `q` near the unit circle (on it the branch is fixed by convention and `cut` is 1).
"""
const CUT_MARGIN = 1e-6

"Mantissas that carry guard digits of their own, so the sum needs no compensation."
_wide_arith(q) = real(q) isa BigFloat || real(q) isa DWNum

"""
`q^n` and `q^{n/2}` as scaled values, accurate to about one rounding: Float64 binary powering errs like
`n·u` and the CG weights reach `n ~ j²`, so a Float64 pass powers in double words and rounds once.
"""
_aqpow(q,n::Int) = _apow(_ascaled(q),n)
_aqpow(q::ComplexF64,n::Int) = _anarrow(_apow(_ascaled(_dwnum(q)),n))
_aqhalfpow(q::Real,n::Int) = iseven(n) ? _aqpow(q,n÷2) : _apow(_asqrt(_ascaled(q)),n)   # odd n: q > 0
_aqhalfpow(q::Complex,n::Int) = _negative_real_axis(q) ?
    _quarter_turn(_aqhalfpow(-real(q),n),n) : _aexp_pow(q,n//2)
_aqhalfpow(q::ComplexF64,n::Int) = _anarrow(_aqhalfpow(_dwnum(q),n))
# Real Float64: Base's `^` is compensated (measured ≤ 0.3 ulp for integer and half-integer exponents) and
# five times cheaper than a double-word powering, whenever the result cannot leave the exponent range.
@inline _pow_in_range(q::Float64,n::Int) = abs(n)*(abs(exponent(q))+1) < 1000
_aqpow(q::Float64,n::Int) =
    _pow_in_range(q,n) ? _ascaled(q^n) : _anarrow(_apow(_ascaled(_dwnum(q)),n))
_aqhalfpow(q::Float64,n::Int) = iseven(n) ? _aqpow(q,n÷2) :
    _pow_in_range(q,n) ? _ascaled(q^(0.5*n)) : _anarrow(_aqhalfpow(_dwnum(q),n))
_anarrow(a::AnalyticScaled) = _ascaled(_narrow(a.m),a.e)

"""
    _analytic_pass(s, tab, w = 0, e2 = 0)

One scaled pass over the rule `s` at `tab.q`. A nonzero `w` multiplies term `z` by `q^{w z}` (one constant
factor per ratio step) and `e2` multiplies the value by `q^{e2/2}`: the q-power weights of the quantum
Clebsch–Gordan coefficient, which a symmetric factorial rule cannot carry (see `qcg.jl`). Both default to
zero, and then the pass is the unweighted one.
"""
function _analytic_pass(s::FactorialSum,tab::AnalyticRuleTable,w::Int=0,e2::Int=0)
    q=tab.q
    ints=tab.ints; facts=tab.facts
    qw=iszero(w) ? _ascaled(one(q)) : _aqpow(q,w)
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
    # Σ_j j|t_j|, the index-weighted mass: term j has been through j ratio steps, so it carries j times
    # the per-step rounding. This accumulator supplies the term-generation error in the running bound.
    wmass=AnalyticScaled(zero(abs(one(q))),0); jstep=0
    for z in s.zlo:s.zhi-1
        # The weight q^w seeds the ratio, so it costs no product of its own.
        ratio=s.alternating ? _aneg(qw) : qw
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
        jstep+=1
        wmass=_aadd(wmass,_amul(AnalyticScaled(oftype(abs(one(q)),jstep),0),_aabs(t)))
    end
    pre,pw,nr,cut=_analytic_prefactor(s,tab)
    if !(iszero(w) && iszero(e2))
        # Inside the scaled representation: q^{e2/2} grows like q^{j²} and only the product is O(1). The
        # first term's weight q^{w·zlo} is common to every term, so it joins the same power: one rounding
        # for the power and one for the product.
        pz=e2+2*w*s.zlo
        iszero(pz) || (pre=_amul(pre,_aqhalfpow(q,pz)))
        nr+=2
    end
    v=_amul(pre,acc)
    s.sign0<0 && (v=_aneg(v))
    cond=iszero(acc.m) ? oftype(real(q),Inf) : _avalue(_adiv(mass,_aabs(acc)))
    wrel=iszero(acc.m) ? oftype(real(q),Inf) : _avalue(_adiv(wmass,_aabs(acc)))
    return v,cond,wrel,pw,nr,cut
end

# A split prime admits i ↦ ι in F_p, so real and complex inputs share word arithmetic.
# One nonzero image is a proof; a vanishing or undefined image still needs exact confirmation.
const ANALYTIC_ZERO_MODULUS = Montgomery(CLASSICAL_PRIMES[2]) # p ≡ 1 mod 4
const ANALYTIC_ZERO_IOTA = to_mont(ANALYTIC_ZERO_MODULUS,
                                  root_of_unity(ANALYTIC_ZERO_MODULUS.p, 4))

"Largest term factorial argument; a common prefactor does not enter the sum's zero decision."
function _sum_max_argument(s::FactorialSum)
    N = 0
    for f in s.fac
        N = max(N, _arg(f, s.zlo), _arg(f, s.zhi))
    end
    return N
end

function _analytic_rational_q(q::Number)
    a = Rational{BigInt}(real(q))
    return iszero(imag(q)) ? a : complex(a, Rational{BigInt}(imag(q)))
end

"Exact rational residue in Montgomery form, or nothing if its denominator vanishes."
function _analytic_residue(m::Montgomery, q::Rational{BigInt})
    a = to_mont(m, UInt64(mod(numerator(q), m.p)))
    b = to_mont(m, UInt64(mod(denominator(q), m.p)))
    iszero(b) && return nothing
    return b == m.one ? a : mont_mul(m, a, mont_inv(m, b))
end

"Prove a regular weighted sum nonzero at the supplied rational q; false means inconclusive."
function _analytic_sum_nonzero_mod(s::FactorialSum, q::Number, w::Int,
                                  m::Montgomery, iota::UInt64)
    is_empty_sum(s) && return false
    r = _analytic_residue(m, real(q))
    r === nothing && return false
    if !iszero(imag(q))
        b = _analytic_residue(m, imag(q))
        b === nothing && return false
        r = mont_add(m, r, mont_mul(m, b, iota))
    end
    iszero(r) && return false
    ri = mont_inv(m, r)
    x = mont_add(m, r, ri)
    qints = Vector{UInt64}(undef, _sum_max_argument(s))
    prev, cur = UInt64(0), m.one
    for n in eachindex(qints)
        # Never perturb q to escape a bad image: that would test a different value.
        iszero(cur) && return false
        qints[n] = cur
        prev, cur = cur, mont_sub(m, mont_mul(m, x, cur), prev)
    end
    # r is nonzero, so exponents may be reduced modulo p-1, including negative w.
    step = mont_pow(m, r, mod(w, Int(m.p - 1)))
    s.alternating && (step = mont_sub(m, UInt64(0), step))
    # With term ratio A/B, Horner gives H_z = 1 + (A/B)H_{z+1}.
    # Store H=P/Q: P ← BQ+AP, Q ← BQ. The checked q-integers keep Q and the
    # omitted first term nonzero, so P alone decides whether this image vanishes.
    P = Q = m.one
    for z in (s.zhi - 1):-1:s.zlo
        A, B = step, m.one
        for f in s.fac
            lo, hi, c = _factor_step(f, z)
            for n in lo:hi
                v = abs(c) == 1 ? qints[n] : mont_pow(m, qints[n], abs(c))
                c > 0 ? (A = mont_mul(m, A, v)) : (B = mont_mul(m, B, v))
            end
        end
        BQ = mont_mul(m, B, Q)
        P = mont_add(m, BQ, mont_mul(m, A, P))
        Q = BQ
    end
    return !iszero(P)
end

"Exact zero decision at the supplied q; modular nonzeros skip the rational fallback."
_analytic_sum_iszero(s::FactorialSum, q::Number, w::Int) =
    _analytic_sum_iszero(s, q, w, ANALYTIC_ZERO_MODULUS, ANALYTIC_ZERO_IOTA)

function _analytic_sum_iszero(s::FactorialSum, q::Number, w::Int,
                             m::Montgomery, iota::UInt64)
    is_empty_sum(s) && return true
    iszero(w) && pairwise_zero(s) && return true
    r = _analytic_rational_q(q)
    _analytic_sum_nonzero_mod(s, r, w, m, iota) && return false
    return _analytic_sum_iszero_exact(s, r, w)
end

"Rational fallback for a regular sum, including an integer q-power weight."
function _analytic_sum_iszero_exact(s::FactorialSum, q::Number, w::Int)
    is_empty_sum(s) && return true
    # Finite binary floats have exact rational coordinates. The nonzero prefactor and first term
    # do not affect the decision; Horner ratios avoid expanding a polynomial or choosing radical roots.
    r = _analytic_rational_q(q)
    x = r + inv(r)
    N = _sum_max_argument(s)
    qints = Vector{typeof(r)}(undef, N)
    prev, cur = zero(r), one(r)
    for n in 1:N
        qints[n] = cur
        prev, cur = cur, x * cur - prev
    end
    step = (s.alternating ? -one(r) : one(r)) * r^w
    total = one(r)
    for z in (s.zhi - 1):-1:s.zlo
        ratio = step
        for f in s.fac
            lo, hi, c = _factor_step(f, z)
            for n in lo:hi
                ratio *= qints[n]^c
            end
        end
        total = one(r) + ratio * total
    end
    return iszero(total)
end

function _analytic_q(q::Number)
    qq=float(q*1.0)
    isfinite(qq) && !iszero(qq) || throw(DomainError(q,"analytic q must be finite and nonzero"))
    # Do not perturb exactly representable singular roots through log/exp.
    # Their cancellations belong to the level/symbolic projection machinery.
    qq == -1 && throw(DomainError(q,"q = -1 is degenerate (q - 1/q = 0); use a nearby q"))
    qq in (im,-im) && throw(DomainError(q,"q = ±i is the root of unity of level 0; use `Level(0)` or `Exact(0)`"))
    # One convention for negative real q and the same point supplied as complex:
    # principal roots of balanced factors, and arg(q)=π for half powers.
    return _negative_real_axis(qq) ? complex(real(qq)) : qq
end

"""
    analytic_value(s, q, T; workspace) -> value

The same evaluation, with `T` as a **floor** on the working precision and output type: the result type is
`promote_type(T, typeof(value))`, so `BigFloat` widens `q`, and a `Complex{BigFloat}` or complex `q` is never
narrowed by the default `T = Float64`.
"""
function analytic_value(s::FactorialSum,q::Number,::Type{T};workspace=nothing,labels=nothing,
                        weight::Tuple{Int,Int}=(0,0),near_edge=nothing) where {T}
    R = real(T)
    qq = (R === BigFloat && !(real(typeof(float(q*1.0))) === BigFloat)) ?
         (q isa Real ? BigFloat(q) : Complex{BigFloat}(q)) : q
    v = analytic_value(s,qq;workspace=workspace,labels=labels,weight=weight,near_edge=near_edge)
    return convert(promote_type(T,typeof(v)),v)
end

"""
Direct, scaled analytic evaluation with compensated sums and precision escalation. For a bare 6j, `labels`
enable the near-edge tier (the column recurrence) after the double word fails; it declines when it cannot
certify, and the arbitrary-precision ladder answers.
"""
function analytic_value(s::FactorialSum,q::Number;workspace=nothing,labels=nothing,
                        weight::Tuple{Int,Int}=(0,0),near_edge=nothing)
    # Keep the positive-real specialization concrete. Promoting negative inputs
    # here avoids a Real/Complex union throughout every positive-real pass.
    q isa Real && q < 0 && return analytic_value(s,complex(q);workspace=workspace,
                                                labels=labels,weight=weight,near_edge=near_edge)
    qq = q isa Real ? real(_analytic_q(q)) : _analytic_q(q)
    T=typeof(qq); R=typeof(real(qq))
    input_q = _negative_real_axis(q) ? qq : q
    is_empty_sum(s) && return zero(T)
    w,e2=weight
    weighted=!(iszero(w) && iszero(e2))
    N=max_argument(s)
    # Conservative operation count: table and prefix-product roundings as well as
    # the accumulated ratio operations. Multiplied by the working precision's
    # unit roundoff and by the sum's own cancellation it bounds the result.
    ops=16*(N+1)*(1+sum(f->abs(Int(f.c)),s.fac;init=0)+
                   sum(p->abs(Int(p.second)),s.pre;init=0))*(s.zhi-s.zlo+1)
    # The weights: one product per step, and powers whose error in the pass's own arithmetic grows like
    # the exponent (`_aqpow`), charged on the same footing.
    weighted && (ops+=16*(s.zhi-s.zlo+1)*(abs(w)*(s.zhi+1)+abs(e2)+2))
    tol=R===BigFloat ? eps(one(real(qq)))*32 : R(64)*eps(R)
    if R===Float64
        # Plain tier, accepted on a running bound (an operation count was ~300× too loose and never passed):
        # term generation c·u·Σ j|t_j|, compensated summation 2u·Σ|t_j|, prefactor and first-term roundings,
        # and for complex q the phase term CPHASE·|P|·|log q| plus the roots and balanced factors (`nr`).
        # Accepting at RTOL_PLAIN gives the same contract as the level and classical kernels.
        mode=POLICY[]
        if !(mode === :strict || mode === :strict_lazy || mode === :compensated_only)
            tabp=_analytic_table(qq,N,workspace)
            vp,kp,wrel,pw,nr,cut=_analytic_pass(s,tabp,w,e2)
            # A weight adds one product per step, and q^w its own rounding.
            cstep=2*sum(f->abs(Int(f.c)),s.fac;init=0)+2+(iszero(w) ? 0 : 2)
            mpre=2*(sum(p->abs(Int(p.second)),s.pre;init=0)+sum(f->abs(Int(f.c)),s.fac;init=0))+2
            phase=iszero(pw) ? zero(R) : R(CPHASE)*R(pw)*abs(log(qq))
            # Table entries are only `pcond·u` accurate, so `pcond` multiplies the count of reads, separately
            # for the prefactor (up to its largest argument) and the sum (up to N). Tables rounded from double
            # words are charged nothing; the condition still drives the near-root warning.
            pc=_prefactor_cond(tabp,N)
            pc >= NEAR_ROOT_COND && _warn_near_root(qq,N)
            pcond=tabp.accurate ? one(R) : R(pc)
            pcpre=tabp.accurate ? one(R) : R(_prefactor_cond(tabp,_pre_max(s)))
            relb=(cstep*wrel*pcond+2*kp+(mpre+2*nr)*pcpre+phase)*eps(R)
            if isfinite(relb) && relb <= RTOL_PLAIN && cut >= CUT_MARGIN
                return T(_avalue(vp))
            end
        end
        # A zero never certifies itself, so a pairwise-cancelling sum would otherwise pay for the double
        # word and then for arbitrary precision. The test is a proof and a few integer comparisons.
        # The zero screens are statements about the unweighted rule; a weight moves terms apart.
        !weighted && pairwise_zero(s) && return zero(T)
        # Double-word tier: u² per operation, so this bound is met for every
        # well-conditioned q and the arbitrary-precision tiers below are reached
        # only near a singularity of the rule.
        vd,kd,_,_,_,_=_analytic_pass(s,_analytic_table(_dwnum(qq),N,workspace),w,e2)
        if isfinite(kd) && ops*eps(DWNum)*kd <= tol
            return T(_aldexp(_narrow(vd.m),vd.e))
        end
        # Near-edge tier before arbitrary precision: O(distance) along the column, asked for this path's own
        # promise RTOL_PLAIN/8 (a stricter target declined the whole unit circle for no promised gain).
        if labels !== nothing && !weighted
            vr = sixj_entry(labels,qq,_column_workspace(workspace);rtol=RTOL_PLAIN/8)
            vr === nothing || return T(vr)
        end
        # The same tier for a rule that brings its own: a coupling coefficient from its Casimir column
        # (`qcg_columns.jl`), which declines with `nothing` when it cannot certify the entry.
        if near_edge !== nothing
            vr = near_edge()
            vr === nothing || return T(vr)
        end
        v=AnalyticScaled(T(_narrow(vd.m)),vd.e)
    else
        v,kappa,_,_,_,_=_analytic_pass(s,_analytic_table(qq,N,workspace),w,e2)
        if isfinite(kappa) && ops*eps(one(real(qq)))*kappa <= tol
            return T(_avalue(v))
        end
    end
    # Every fixed-precision tier has failed. If the sum is the zero function of q there is no bound for
    # any precision to certify, and the ladder below would only double its bits until it gave up; the
    # generic-q filter and exact polynomial check decide that first.
    !weighted && is_generic_zero(s) && return zero(T)
    # Rebuild from the supplied q at each precision, never from rounded table
    # entries. A tier that certifies its own bound is accepted on its own: the
    # machine-precision value is not a reliable witness, so requiring the two to
    # agree only forced a second arbitrary-precision pass. A rounded zero, or a sum
    # still buried in rounding noise at two widths, triggers exact confirmation at the supplied q.
    zero_checked=false
    cancelled=false
    bits=max(128,precision(R)+32)
    for _ in 1:8
        next,condition,_,_,_,_=setprecision(BigFloat,bits) do
            qb=input_q isa Real ? BigFloat(input_q) : Complex{BigFloat}(input_q)
            _analytic_pass(s,_analytic_table(qb,N,workspace),w,e2)
        end
        if !iszero(next.m) && isfinite(condition) && ops*ldexp(one(condition),-bits)*condition < tol/4
            return T(_avalue(next))
        end
        small = !isfinite(condition) || exponent(condition) >= bits - 64
        if !zero_checked && (iszero(next.m) || (small && cancelled))
            _analytic_sum_iszero(s,input_q,w) && return zero(T)
            zero_checked=true
        end
        cancelled=small
        bits*=2
    end
    throw(ErrorException("analytic evaluation did not converge; increase input precision or use an exact level target"))
end
