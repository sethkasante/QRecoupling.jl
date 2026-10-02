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
    circle::Bool                      # |q| = 1 to within the input's own rounding; see `_on_unit_circle`
    ints::Vector{AnalyticScaled{T}}
    inverses::Vector{AnalyticScaled{T}}
    facts::Vector{AnalyticScaled{T}}  # index n+1 represents [n]!
    balanced::Vector{AnalyticScaled{T}}
    roots::Vector{AnalyticScaled{T}}  # √Ψ_d on the package's branch; complex q only
    condmax::Vector{Float64}          # prefix maxima of the q-integer condition; see `_prefactor_cond`
end

"""
    _on_unit_circle(q) -> Bool

Whether `q` is on the unit circle up to the rounding of its own Float64 representation.

This is a question about the *input*, not about the tier evaluating it: a `q` the caller wrote as
`cispi(t)` has `|q| = 1 ± 2⁻⁵³` because that is as close as a pair of Float64s gets, and widening it to a
double word or to BigFloat reproduces the same near-miss rather than curing it. So the same tolerance is
used at every precision, which is what makes the tiers agree with one another.

On the circle every balanced factor `Ψ_d` is *exactly real*, and a negative one sits on the branch cut of
the square root. Left to the arithmetic, which of `±i√|Ψ_d|` comes back is then decided by the sign of the
rounding error in `|q| − 1` — measured, **21 of 72** unit-circle 6j symbols flipped sign under a one-ulp
perturbation of `q` that does not change the point being asked about. `_analytic_prefactor` therefore
discards that imaginary part and roots the real value, which fixes the branch at `+i√|Ψ_d|`; see there.
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
`√Ψ_d` for every `d`, on the branch the package takes — one principal root per balanced factor, and the
imaginary part discarded first when `q` is on the unit circle (see `_on_unit_circle`).

Cached because it is the *only* transcendental left in a complex-q prefactor and it does not depend on the
rule: ten rooted factors at a complex `sqrt` apiece were about 300 ns of a 900 ns symbol, recomputed on
every call at the same `q`. Built on first use, like `balanced`, and never needed on the positive real
axis where the prefactor takes one root of the assembled product instead.
"""
function _table_roots(tab::AnalyticRuleTable)
    length(tab.roots) == length(tab.ints) && return tab.roots
    psi = _table_balanced(tab)
    r = tab.circle ? [_asqrt(_deimag(p)) for p in psi] : [_asqrt(p) for p in psi]
    tab.roots = r
    return r
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
    AnalyticRuleTable(q,bits,_on_unit_circle(q),ints,empty,facts,copy(empty),copy(empty),Float64[])
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
    AnalyticRuleTable(q,precision(DWNum),_on_unit_circle(q),ints,empty,facts,copy(empty),copy(empty),Float64[])
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
Process-wide analytic tables, for callers that pass no workspace.

Building the table is 20–35% of a generic-q symbol and it depends on `q` and the capacity alone, so a
caller that sweeps labels at one `q` — which is what `q6j(js...; q = z)` in a loop is — was paying for it
on every call. A workspace already avoided that; this gives the same saving to the scalar API, which is
where most calls come from. Measured 1.4–2.2× on a cold scalar call, the whole of the gap between the
`total` and `warm-ws` columns of `dev/results/analytic_path.md` §4.

Entries are shared between tasks. They are only ever *extended* — `inverses`, `balanced` and `condmax` are
built in a local and then assigned whole — so a task either sees the empty field and rebuilds it or sees a
finished one, and two tasks that race produce the same content. Nothing here needs a lock beyond the cache
itself, which `empty_caches!` clears.
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
`Π Ψ_d^{e_d}` for one balanced monomial, from the table's balanced factors.

A top-level function rather than the closure this used to be. The closure assigned to a variable that the
enclosing scope also assigned to, so Julia shared them and boxed it — **35 heap allocations per prefactor**
on the plain real path, which is a third of the cost of an entire generic-q symbol and does not exist on
the level or classical paths.
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

The Ψ-basis exponents of a rule's prefactor `Π [n]!^c`, and its residual q-power.

`[n]! = Π_{d≥2} Ψ_d^{⌊n/d⌋}` with `Ψ_d = q^{-φ(d)}Φ_d(q²)`, and `[n]` contributes `q^{1−n}` on top, so one
pass over the rule's factorials gives both — the same numbers a `CycloBuffer` accumulates, and the same
split into a square part and a square-free radical that `snapshot_square_root` makes, but as a flat
exponent array instead of a buffer plus two monomials carrying two sparse exponent vectors. Those were
**9 heap allocations in every complex-q prefactor**, on the path every symbol off the positive real axis
takes; this is one.
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

The prefactor, together with what a running error bound needs to know about *how* it was computed.
`phase_weight` is `|pr + pd/2|`, the exponent of the bare `q` power: on the complex path that power is
`exp(phase·log q)`, whose relative error is governed by `|phase·log q|·u` rather than by any count of
roundings, which is exactly the term the plain tier was missing. `nroots` counts the individual square
roots and balanced-factor multiplications, which the operation count of the caller does not see.
"""
function _analytic_prefactor(s,tab)
    q=tab.q
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
    # Complex q: the square root of the radical must be taken *factor by factor*, never of the
    # assembled product. √ is not multiplicative across its branch cut, so the principal root of
    # ∏ Ψ_d differs from ∏ √Ψ_d by a sign that depends on how the individual phases add up. With a
    # root of the product, the two sides of a coherence identity assemble different products and
    # their radicals no longer cancel; with one root per Ψ_d the choice depends only on the *set* of
    # factors, which both sides share, so the cancellation is exact. Measured over 40 label sets:
    # Biedenharn–Elliott 40/40 on and off the unit circle, against 9–36/40 for a root of the product
    # (`dev/results/user_facing_exact.md` §2, `dev/prototypes/branch_consistent_prefactor.jl`).
    # `rad` is square-free, so every exponent here is ±1.
    sv=_ascaled(one(q))
    sq=_table_roots(tab)
    # Distance of each rooted factor from the branch cut, relative to its own size — for `q` *off* the
    # unit circle, where a small imaginary part of Ψ_d is a fact about q and decides the root honestly.
    # The plain tier is refused when it cannot resolve that sign (`CUT_MARGIN`).
    #
    # **On** the circle there is no sign to resolve: every Ψ_d is exactly real, a negative one lies on
    # the cut, and what the arithmetic returns is decided by the rounding error in |q| − 1 rather than by
    # the point being asked about (21 of 72 symbols flipped under a one-ulp perturbation; see
    # `_on_unit_circle`). The imaginary part is therefore discarded and the real value rooted, which fixes
    # the branch at +i√|Ψ_d| for every tier at once. That choice is not arbitrary:
    #
    #  * it is continuous in θ = arg q — √Ψ_d runs down to 0 along the reals and back up along +i, so the
    #    symbol has no jumps between roots of unity, which the ulp-dependent branch did not manage;
    #  * it agrees with the level path, where the radicand of an admissible symbol is positive and no
    #    choice arises, so the limit θ → π/h reproduces `Level(k)`;
    #  * it depends only on the *set* of Ψ_d, which both sides of a coherence identity share, so the
    #    radicals still cancel exactly — the property the factor-by-factor rule exists for.
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
Weight of the phase term `exp(P·log q)` in the plain tier's running bound. Three roundings reach the
exponent — the logarithm, the multiplication by `P`, and the exponential's own argument reduction — and
one more is allowed for the scaled representation, so four is the constant that makes the bound hold with
margin; it is validated by measurement, not by the count alone (`dev/results/numeric_audit.md` §3).
"""
const CPHASE = 4

"""
    _prefactor_cond(q, N) -> Float64

How accurate a single table entry is, in units of `u`. Every q-integer is built from `1 − q^{-2n}`, and
subtracting two numbers of modulus one loses digits in proportion to `|q^{-2n}| / |1 − q^{-2n}|`. Near a
root of unity of order m that blows up at `n = m`, which is exactly where the plain tier's value went
wrong: the *sum's* condition number says nothing about a prefactor assembled from inaccurate q-integers.

Real `q` is excluded on purpose. Its tier was measured and accepted with this factor absent, and adding it
would reject near-classical real `q` (at `q = 0.99` the factor is ≈ 50) for no measured gain in accuracy —
a speed regression, not a correctness fix.

The loop is short: once `|q^{-2n}|` leaves a neighbourhood of 1 the ratio is bounded by 2 and stays there,
so only `|q| = 1` runs the full length, and there nothing can overflow.
"""
function _prefactor_cond(q,N::Int)
    q isa Real && return 1.0
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

The same number, read from the table instead of recomputed. Off the unit circle the loop above stops after
a few terms and either form is free; **on** it `|q^{-2n}|` never leaves the neighbourhood of 1, so the loop
runs its full length on every call — as expensive as building the table it is guarding. The prefix maxima
are a property of `q` alone, so they are computed once per table and indexed.
"""
function _prefactor_cond(tab::AnalyticRuleTable,N::Int)
    tab.q isa Real && return 1.0
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

Above this prefactor condition number, `q` is a root of unity that could not be written down exactly.

`_prefactor_cond` is `max_n |q^{-2n}| / |1 − q^{-2n}|`, so `c ≥ 2⁻⁸/u` means some `[n]` agrees with zero
to within a few of its own last bits: `q` is within rounding of a root of unity of order at most `2N`. The value
returned is still the value *at the Float64 `q` that was passed* — a 512-bit rerun of the same `q` agrees
with it to `2e−16` — but `q` itself is then only a `2⁻⁵³` approximation of the point the caller meant, and
perturbing it by one ulp changes the answer by a factor of order one. `q6j(1,1,1,1,1,1; q = cispi(1/3))`
returns `5.8e15 + 1.0e16im`; at the level it names, `q6j(1,1,1,1,1,1; k = 1)` is `0`.

A deliberate probe of the neighbourhood is not caught: at `cispi(1/7)·(1 + 1e−10)` the smallest `|1 − q^{-2n}|`
is `1.4e−9`, seven orders above the threshold, and those values are certified in the ordinary way.
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
How far a rooted factor must sit from the branch cut, relative to its own magnitude, before a plain
Float64 prefactor is trusted. Below this the principal root's *sign* is set by rounding rather than by the
value, and the tier falls through to the double word.

This is now a guard for `q` **near** the unit circle, not on it. On the circle the imaginary part of a
`Ψ_d` is rounding noise rather than information, so `_analytic_prefactor` discards it and the branch is
fixed by convention (see `_on_unit_circle`); there is no sign left for a wider mantissa to resolve, and
`cut` comes back as 1. Just off the circle a small imaginary part is a fact about `q` and does decide the
root, which is what this refuses to guess in Float64.
"""
const CUT_MARGIN = 1e-6

"Mantissas that carry guard digits of their own, so the sum needs no compensation."
_wide_arith(q) = real(q) isa BigFloat || real(q) isa DWNum

"""
`q^n` and `q^{n/2}` as scaled values, accurate to about one rounding of the pass's own precision.

Binary powering in Float64 carries a relative error that grows like `n·u` (each squaring doubles the error
already there), and the Clebsch–Gordan weights reach `n ~ j²`. A Float64 pass therefore powers in a double
word and rounds once; the double-word and arbitrary-precision passes power in their own arithmetic, where
`n·u` is far below what they certify.
"""
_aqpow(q,n::Int) = _apow(_ascaled(q),n)
_aqpow(q::ComplexF64,n::Int) = _anarrow(_apow(_ascaled(_dwnum(q)),n))
_aqhalfpow(q::Real,n::Int) = iseven(n) ? _aqpow(q,n÷2) : _apow(_asqrt(_ascaled(q)),n)   # odd n: q > 0
_aqhalfpow(q::Complex,n::Int) = _aexp_pow(q,n//2)
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
    # the per-step rounding. This is the accumulator Theorem 1 needs, and it costs one add per term.
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

"""
    analytic_value(s, q, T; workspace) -> value

The same evaluation, with `T` as a **floor** on the working precision and on the output type.

Without this method the `T` keyword was silently dropped at generic `q`: `q6j(js...; q = 0.7, T = BigFloat)`
returned a `Float64`, and a batch returned a `Vector{Float64}`, while the identical call at a level
honoured `T`. Asking for `BigFloat` now widens `q` and runs the whole ladder there.

`T` is a floor and not an exact output type, because at generic `q` the *parameter* also carries a
precision and a `q` the caller gave in `Complex{BigFloat}` must not be narrowed by the default
`T = Float64`. So the result type is `promote_type(T, typeof(value))`: asking for more gives more, asking
for less than `q` already carries gives what `q` carries. A real `T` with a complex `q` widens the same
way, which is why `q6j(js...; q = 0.8 + 0.3im)` keeps returning a `ComplexF64`.
"""
function analytic_value(s::FactorialSum,q::Number,::Type{T};workspace=nothing,labels=nothing,
                        weight::Tuple{Int,Int}=(0,0)) where {T}
    R = real(T)
    qq = (R === BigFloat && !(real(typeof(float(q*1.0))) === BigFloat)) ?
         (q isa Real ? BigFloat(q) : Complex{BigFloat}(q)) : q
    v = analytic_value(s,qq;workspace=workspace,labels=labels,weight=weight)
    return convert(promote_type(T,typeof(v)),v)
end

"""
Direct, scaled analytic evaluation with compensated sums and precision escalation.

`labels` are the doubled 6j labels when the rule *is* a bare 6j, and they buy the near-edge tier: the
symbol as one entry of its column, filled by the three-term recurrence in `O(distance)` with no Racah-sum
cancellation at all. It is offered the case only after the double word has failed, exactly as
`level_escalate` does at a level, and it declines whenever its own estimate cannot certify the entry — so
the arbitrary-precision ladder below is still what answers when the recurrence is unsafe.
"""
function analytic_value(s::FactorialSum,q::Number;workspace=nothing,labels=nothing,
                        weight::Tuple{Int,Int}=(0,0))
    qq=_analytic_q(q); T=typeof(qq); R=typeof(real(qq))
    is_empty_sum(s) && return zero(T)
    w,e2=weight
    weighted=!(iszero(w) && iszero(e2))
    # q^{e2/2} for odd e2 is not real at negative q; the weighted rules are evaluated there as complex.
    weighted && q isa Real && signbit(qq) && isodd(e2) &&
        return analytic_value(s,complex(q);workspace=workspace,weight=weight)
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
        # Plain tier, accepted on a *running* bound rather than an operation count.
        #
        # The operation-count bound below (`ops`) is ~300× looser than the work actually done, because it
        # charges every term the worst case: at j = 5 it is of order 10⁴ against an index-weighted
        # c·Σj|t_j|/|Σ| of order 50κ. Measured, it accepted a plain pass for *0 of 13* label sets in every
        # regime, so the plain tier was dead code and every generic-q call paid for double-word tables and
        # arithmetic — 2.7–4.4× on the pass and 4.9–7.7× on the table (`dev/results/numeric_audit.md` §3).
        #
        # The bound here has the shape of Theorem 1: the term-generation error c·u·Σ j|t_j|, the
        # summation error (the plain pass is compensated, so 2u·Σ|t_j| covers it), and the prefactor and
        # first-term roundings. Accepting at `RTOL_PLAIN` makes the contract the same one the level and
        # classical kernels already offer, rather than a second, stricter one that nothing could meet.
        # A value that fails it falls through to exactly the tiers that ran before.
        #
        # Complex q needs one more term. Off the positive real axis the prefactor ends in
        # `_aexp_pow(q, P) = exp(P·log q)` with a balanced phase `P = pr + pd/2` that grows with the
        # labels, and the error of that is not a count of roundings: perturbing the exponent by δ
        # multiplies the result by exp(δ), so the relative error is |P·log q|·u up to a small constant.
        # Charging `CPHASE·|P|·|log q|`, plus the individual square roots and balanced factors that the
        # operation count never saw (`nr`), makes the bound cover the complex path as well. On the
        # positive real axis `pw` and `nr` are zero and this is exactly the bound measured before, so
        # the real tier is unchanged. `dev/results/numeric_audit.md` §3.
        mode=POLICY[]
        if !(mode === :strict || mode === :strict_lazy || mode === :compensated_only)
            tabp=_analytic_table(qq,N,workspace)
            vp,kp,wrel,pw,nr,cut=_analytic_pass(s,tabp,w,e2)
            # A weight adds one product per step, and q^w its own rounding.
            cstep=2*sum(f->abs(Int(f.c)),s.fac;init=0)+2+(iszero(w) ? 0 : 2)
            mpre=2*(sum(p->abs(Int(p.second)),s.pre;init=0)+sum(f->abs(Int(f.c)),s.fac;init=0))+2
            phase=iszero(pw) ? zero(R) : R(CPHASE)*R(pw)*abs(log(qq))
            # Every table entry a rule reads is only `pcond·u` accurate, so `pcond` multiplies the
            # *count* of reads — the prefactor's factorials and balanced factors, and equally the ratio
            # steps of the sum. Charging it on the prefactor alone left the bound short: over 8,000
            # (labels, q) pairs it understated the measured error in 5 of them, worst by 3.1×, all on the
            # unit circle where `pcond` is the several-fold thing it is, and one accepted case came within
            # 3% of the promise. With it on both terms nothing understates. On the positive real axis
            # `pcond` is 1 and the bound is bit-for-bit the one measured before.
            #
            # The two terms get their *own* condition number, because they do not read the same entries.
            # The prefactor's factorials stop at the largest argument in `s.pre`; the sum's ratio steps go
            # to `N`. On the unit circle the worst `[n]` is often above the prefactor's reach, and
            # charging the whole bound at `condmax[N]` then rejects a pass whose prefactor never touched
            # the bad index — which is most of what kept the unit circle on the double word at large spin.
            pcond=R(_prefactor_cond(tabp,N))
            pcpre=R(_prefactor_cond(tabp,_pre_max(s)))
            pcond >= NEAR_ROOT_COND && _warn_near_root(qq,N)
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
        # The near-edge tier, before arbitrary precision: O(distance) along the column instead of a sum
        # whose digits have gone. It returns `nothing` rather than a value it cannot stand behind.
        #
        # It is asked for this path's own promise rather than its stricter default. `RTOL_PLAIN/8` is what
        # every other value returned here is certified to, with a factor of eight in hand for an estimate
        # that is a measured heuristic; asking for 1e−14 instead made the tier decline the whole unit
        # circle, where the table's conditioning puts a floor of about 5e−14 on it, and the arbitrary
        # precision tier then cost 46 µs at j = 30 to deliver digits nobody was promised.
        if labels !== nothing && !weighted
            vr = sixj_entry(labels,qq,_column_workspace(workspace);rtol=RTOL_PLAIN/8)
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
    # generic-q screen decides that first, as the modular screens do at a level and at q = 1.
    !weighted && is_generic_zero(s) && return zero(T)
    # Rebuild from the supplied q at each precision, never from rounded table
    # entries. A tier that certifies its own bound is accepted on its own: the
    # machine-precision value is not a reliable witness, so requiring the two to
    # agree only forced a second arbitrary-precision pass. An exactly cancelling
    # sum carries no bound, and there two tiers must agree instead.
    previous=v
    bits=max(128,precision(R)+32)
    for _ in 1:8
        next,condition,_,_,_,_=setprecision(BigFloat,bits) do
            qb=q isa Real ? BigFloat(q) : Complex{BigFloat}(q)
            _analytic_pass(s,_analytic_table(qb,N,workspace),w,e2)
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
