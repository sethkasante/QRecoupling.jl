# ---------------------------------------------------------------------------------
#  q-integers, q-factorials and q-binomials, on the same footing as the symbols
#
#  These are not sums: each is a bare product `Π [n]!^{c}`, which is a `FactorialSum` whose range is a
#  single point and whose factor list is empty. Saying so routes them through the same evaluator as
#  `q6j` — the level kernel with its running bound, the classical integer path, the analytic ladder at
#  generic q, and the real basis when an exact value is asked for — instead of through the monomial
#  projection they used to take. One consequence is worth stating: `qint(5)` is now a *number*, the
#  classical value, exactly as `q6j(…)` is. The cyclotomic monomial is still there, as
#  `qint(Symbolic(), 5).dcr`.
# ---------------------------------------------------------------------------------

"""
A product of q-factorials `Π [n]!^{c}` as a factorial rule: one term, a prefactor, and nothing else.

Every route the symbols use accepts this, so a q-binomial gets the level kernel's certified bound and the
generic-q ladder for free rather than having its own arithmetic.
"""
_qfact_rule(pairs) = FactorialSum(0:0; prefactor = pairs)

"The empty rule, whose value is the zero of whichever carrier is asked for."
const _ZERO_PRODUCT_RULE = FactorialSum((), false, (), false, Int8(0), 0, -1)

"""
    _product_value(pairs, k, q, exact, T; workspace) -> value

`Π [n]!^{c}` wherever the symbols evaluate. `pairs` may be empty, which is the value 1, and `nothing`
stands for the empty product, which is 0 in whatever carrier the target names.

Numeric products use direct monomial evaluation: there is no sum to certify, so table lookups avoid
the fixed overhead of the factorial-rule summation machinery. `Symbolic()` returns the rule, while an
exact level value is built from its q-factorial content in the real basis. The interface remains shared
with the symbols, with arithmetic suited to each target.
"""
function _product_value(pairs, k, q, exact::Bool, ::Type{T}; workspace = nothing) where {T}
    if exact && !isnothing(k)
        kk = Int(k); kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
        pairs === nothing && return zero(ExactX, kk)
        return _qfact_exactx(pairs, kk)
    end
    if _is_classical(q) && !exact
        # At q = 1 the rule's kernel works in integer ratios, so `qint(5)` is 5.0 and not
        # 4.999999999999999 — which is what the monomial route returns, since it multiplies floating
        # values of Φ_d(1). Eight times the cost of that route and still a microsecond; for the one
        # place a user checks a package's arithmetic by eye, the right trade.
        return classical_value(pairs === nothing ? _ZERO_PRODUCT_RULE : _qfact_rule(pairs), T;
                               workspace = workspace)
    end
    return qeval(_pairs_mono(pairs); k = k, q = q, exact = exact, T = T)
end

# ---------------------------------------------------------------------------------
#  Numbers straight from the tables, for Float64 targets
#
#  A q-integer is one table entry, a q-factorial one split entry, a q-binomial three. Reading them is what
#  the symbols do already, so the cheap route is also the accurate one: the level table holds [n] correctly
#  rounded, and the split factorials hold [n]! as a correctly rounded mantissa and an exact exponent. The
#  monomial route this replaces built a cyclotomic monomial on every call (7–10 allocations, and 21 with
#  1.2 µs for a classical `qint(5)`), multiplied floating values of Φ_d, and near q = 1 lost up to 4e-10.
#  Anything these paths do not cover — other element types, exact targets, a level factorial that reaches
#  [h] = 0 — returns `nothing` and takes the route it always took.
# ---------------------------------------------------------------------------------

"[n] at q = e^{iπ/(k+2)}, correctly rounded from the level table, for every integer n by periodicity."
@inline function _qint_level(n::Int, k::Int)
    h = k + 2
    r = mod(n, 2h)
    (r == 0 || r == h) && return 0.0
    tab = qint_tables(Float64, k)
    return r < h ? @inbounds(tab.q[r]) : -@inbounds(tab.q[2h - r])
end

"""
[n] at a real positive or complex q, as q^{1−n}(q^{2n} − 1)/(q² − 1).

When q^{2n} is near 1 that difference is the whole value and cancels in floating point, so it is taken
from `expm1` of 2n·log q instead; otherwise the integer powers are more accurate than exp(n·log q), whose
error grows like n·|log q|·u. `nothing` where neither form applies (a negative real q).
"""
function _qint_analytic(n::Int, q::Union{Float64,ComplexF64})
    q isa Float64 && !(q > 0) && return nothing
    L = log(q)
    if abs(2n * L) < 1
        return q^(1 - n) * (expm1(2n * L) / expm1(2L))
    end
    # Written so that the large power is the one outside: for |q| > 1 the mirrored form, since
    # q^{1−n}·q^{2n} would be 0·∞ long before the value itself overflows.
    abs(q) <= 1 && return q^(1 - n) * ((q^(2n) - 1) / ((q - 1) * (q + 1)))
    return q^(n - 1) * ((1 - inv(q^(2n))) / (((q - 1) * (q + 1)) / (q * q)))
end

# ---- q-factorials at generic q, without intermediate overflow ----
#
# [n]! = q^{−n(n−1)/2} ∏_{j=2}^{n} g_j,  g_j = (q^{2j} − 1)/(q² − 1),
#
# kept as a mantissa and a binary exponent throughout, so a value that fits in a Float64 is returned even
# when its factors or its prefactor do not: the monomial route overflowed `qfact(30; q = 2.7)` (true value
# 3.1·10¹⁸⁹) to Inf, and 39 representable cases of a test grid to Inf or 0. Each g_j is taken the way
# `_qint_analytic` takes [n] — from `expm1` when q^{2j} is near 1, from an exact-base power otherwise, and
# for |q| > 1 as q^{2j}(1 − q^{−2j}) with the large power split — so its error stays a few units whatever
# the size of q^{2j}.

"`x·2^n`, for complex `x` too."
@inline _xldexp(x::Float64, n::Int) = ldexp(x, n)
@inline _xldexp(x::ComplexF64, n::Int) = ComplexF64(ldexp(real(x), n), ldexp(imag(x), n))

"Scale `m` by a power of two into [1/2, 1) in its largest component; `e` absorbs the exponent."
@inline function _snorm(m::Float64, e::Int)
    iszero(m) && return m, e
    fr, ex = frexp(m)
    return fr, e + ex
end
@inline function _snorm(m::ComplexF64, e::Int)
    a = max(abs(real(m)), abs(imag(m)))
    (iszero(a) || !isfinite(a)) && return m, e
    ex = exponent(a) + 1
    return ComplexF64(ldexp(real(m), -ex), ldexp(imag(m), -ex)), e + ex
end

# ---- double-word values, real (hi, lo) or complex (re_hi, re_lo, im_hi, im_lo) ----
#
# Powers are where a float product loses most: binary powering doubles the accumulated relative error at
# every squaring, so z^t costs ~t·u — 1.6e-13 for cis(0.7)^4950, where Julia's compensated real `^` gives
# 2e-17 but its complex `^` does not. Carried in double words the same powering costs ~t·u², and so does
# the running q^{2j} of the factorial below; each is rounded to a double once, where it is used.

@inline _tprod(a::Float64, b::Float64) = (p = a * b; (p, fma(a, b, -p)))
@inline function _tsum(a::Float64, b::Float64)
    s = a + b; bb = s - a
    return s, (a - (s - bb)) + (b - bb)
end
@inline _qtsum(a::Float64, b::Float64) = (s = a + b; (s, b - (s - a)))
@inline function _dwmul(ah, al, bh, bl)
    p, e = _tprod(ah, bh)
    return _qtsum(p, e + (ah * bl + al * bh))
end
@inline function _dwadd(ah, al, bh, bl)
    s, e = _tsum(ah, bh)
    return _qtsum(s, e + (al + bl))
end

_dwof(q::Float64) = (q, 0.0)
_dwof(q::ComplexF64) = (real(q), 0.0, imag(q), 0.0)
_qdwone(::Float64) = (1.0, 0.0)
_qdwone(::ComplexF64) = (1.0, 0.0, 0.0, 0.0)
@inline _qdwmul(a::NTuple{2,Float64}, b::NTuple{2,Float64}) = _dwmul(a[1], a[2], b[1], b[2])
@inline function _qdwmul(a::NTuple{4,Float64}, b::NTuple{4,Float64})
    rr = _dwmul(a[1], a[2], b[1], b[2]); ii = _dwmul(a[3], a[4], b[3], b[4])
    ri = _dwmul(a[1], a[2], b[3], b[4]); ir = _dwmul(a[3], a[4], b[1], b[2])
    re = _dwadd(rr[1], rr[2], -ii[1], -ii[2]); im = _dwadd(ri[1], ri[2], ir[1], ir[2])
    return (re[1], re[2], im[1], im[2])
end
@inline _dwsub1(a::NTuple{2,Float64}) = _dwadd(a[1], a[2], -1.0, 0.0)
@inline _dwsub1(a::NTuple{4,Float64}) = ((r = _dwadd(a[1], a[2], -1.0, 0.0)); (r[1], r[2], a[3], a[4]))
@inline _dwscale(a::NTuple{N,Float64}, k::Int) where {N} = map(x -> ldexp(x, k), a)
@inline _dwround(a::NTuple{2,Float64}) = a[1] + a[2]
@inline _dwround(a::NTuple{4,Float64}) = ComplexF64(a[1] + a[2], a[3] + a[4])
@inline function _dwnorm(a::NTuple{N,Float64}, e::Int) where {N}
    mx = N == 2 ? abs(a[1]) : max(abs(a[1]), abs(a[3]))
    (iszero(mx) || !isfinite(mx)) && return a, e
    ex = exponent(mx) + 1
    return _dwscale(a, -ex), e + ex
end

# Renormalising at every step costs a frexp or two ldexps; the exponent range only needs protecting when a
# value has drifted far, so the loops below renormalise lazily, outside [2^−400, 2^400].
@inline function _lnorm(m::Union{Float64,ComplexF64}, e::Int)
    a = m isa Float64 ? abs(m) : max(abs(real(m)), abs(imag(m)))
    (0x1p-400 < a < 0x1p400) && return m, e
    return _snorm(m, e)
end
@inline function _ldwnorm(a::NTuple{N,Float64}, e::Int) where {N}
    mx = N == 2 ? abs(a[1]) : max(abs(a[1]), abs(a[3]))
    (0x1p-400 < mx < 0x1p400) && return a, e
    return _dwnorm(a, e)
end

"b^t (t ≥ 0) for a double-word base with exponent `be`, as a split double word: error ~t·u², not t·u."
function _dwpow(b::NTuple{N,Float64}, be::Int, t::Int) where {N}
    r = N == 2 ? (1.0, 0.0) : (1.0, 0.0, 0.0, 0.0); re = 0
    b, be = _dwnorm(b, be)
    while t > 0
        if isodd(t)
            r = _qdwmul(r, b); re += be; r, re = _ldwnorm(r, re)
        end
        t >>= 1
        if t > 0
            b = _qdwmul(b, b); be *= 2; b, be = _dwnorm(b, be)
        end
    end
    return r, re
end

"""
[n]! at a real positive or complex q as a split value (m, e), or `nothing` where the formula does not apply.

    [n]! = q^{−n(n−1)/2} ∏_{j=2}^{n} (q^{2j} − 1) / (q² − 1)^{n−1}

The running q^{2j}, q² − 1 and both powers are double words, so the only roundings of size u are one per
factor; a denominator rounded once and divided n times would make its error systematic, n·u.
"""
function _qfact_split(n::Int, q::T) where {T<:Union{Float64,ComplexF64}}
    q isa Float64 && !(q > 0) && return nothing
    (0x1p-300 < abs(q) < 0x1p300) || return nothing       # q² itself must be a double
    n <= 1 && return one(T), 0
    q2, q2e = _dwpow(_dwof(q), 0, 2)
    d, de = _dwnorm(_dwsub1(_dwscale(q2, q2e)), 0)          # q² − 1, still a double word
    iszero(d[1]) && return nothing
    p, pe = q2, q2e                                         # q^{2j}, starting at j = 1
    m, e = one(T), 0
    q2, q2e = _dwscale(q2, q2e), 0                          # in range: |q|² is a modest number
    for j in 2:n
        p = _qdwmul(p, q2); pe += q2e; p, pe = _ldwnorm(p, pe)
        g, ge = _pm1(p, pe)
        m, e = _lnorm(m * g, e + ge)
    end
    dm, dme = _dwpow(d, de, n - 1)
    m, e = _snorm(m / _dwround(dm), e - dme)
    qm, qme = _dwpow(_dwof(q), 0, n * (n - 1) ÷ 2)
    return _snorm(m / _dwround(qm), e - qme)
end

"""
The q-binomial [n choose m] at a real positive or complex q as a split value, as one product of m ratios

    q^{−m(n−m)} ∏_{j=1}^{m} (q^{2(n−m+j)} − 1) / (q^{2j} − 1),

with the running powers in double words. Dividing three factorials would cost three times as much and
round three times.
"""
function _qbinomial_split(n::Int, m::Int, q::T) where {T<:Union{Float64,ComplexF64}}
    q isa Float64 && !(q > 0) && return nothing
    (0x1p-300 < abs(q) < 0x1p300) || return nothing       # q² itself must be a double
    m = min(m, n - m)
    m <= 0 && return one(T), 0
    q2, q2e = _dwpow(_dwof(q), 0, 2)
    q2, q2e = _dwscale(q2, q2e), 0
    a, ae = _dwpow(q2, q2e, n - m)                       # q^{2(n−m)}, advanced to q^{2(n−m+j)}
    b, be = _qdwone(q), 0                                  # q^{2j}
    val, e = one(T), 0
    for j in 1:m
        a = _qdwmul(a, q2); ae += q2e; a, ae = _ldwnorm(a, ae)
        b = _qdwmul(b, q2); be += q2e; b, be = _ldwnorm(b, be)
        num, ne = _pm1(a, ae); den, dne = _pm1(b, be)
        iszero(den) && return nothing
        val, e = _lnorm(val * (num / den), e + ne - dne)
    end
    qm, qme = _dwpow(_dwof(q), 0, m * (n - m))
    return _snorm(val / _dwround(qm), e - qme)
end

"p − 1 for a split double word p = (p, pe), as a split double: exact enough in every regime of pe."
@inline function _pm1(p::NTuple{N,Float64}, pe::Int) where {N}
    pe == 0 && return _dwround(_dwsub1(p)), 0
    mx = N == 2 ? abs(p[1]) : max(abs(p[1]), abs(p[3]))
    lg = exponent(mx) + pe                                  # log2 |p·2^pe|, roughly
    lg > 60 && return _dwround(p), pe                       # 1 is below the last place
    lg < -60 && return (N == 2 ? -1.0 : ComplexF64(-1.0)), 0 # p is below the last place of 1
    return _dwround(_dwsub1(_dwscale(p, pe))), 0
end

"A split value raised to an integer power p, and returned as a number (Inf or 0 only if the value is)."
function _split_value(m::T, e::Int, p::Int) where {T}
    if p < 0
        m = inv(m); e = -e; p = -p
    end
    rm, re = one(T), 0
    while p > 0
        if isodd(p)
            rm *= m; re += e; rm, re = _snorm(rm, re)
        end
        p >>= 1
        if p > 0
            m *= m; e *= 2; m, e = _snorm(m, e)
        end
    end
    return _xldexp(rm, re)
end

"Π [n]!^c from the split tables at a level; `nothing` if a factorial reaches [h] = 0."
function _qfacts_level(k::Int, pairs::NTuple{N,Tuple{Int,Int}}) where {N}
    tab = qint_tables(Float64, k)
    m, e = 1.0, 0
    for (n, c) in pairs
        n <= k + 1 || return nothing
        n >= 1 && ((m, e) = _split_mul(m, e, tab, n, c))
    end
    return ldexp(m, e)
end

"""
The Float64 value of a q-integer power, q-factorial power or q-binomial, or `nothing` to take the general
route. `what` is `:int`, `:fact` or `:binomial`; `a`, `b` are (n, p) or (n, m).
"""
function _qnumber_float(what::Symbol, a::Int, b::Int, k, q)
    # the same conventions as the pair constructors: negative n is zero, [0] = [1] = 1, [0]! = [1]! = 1
    if what === :binomial
        (b < 0 || b > a) && return 0.0
        (b == 0 || b == a) && return 1.0
    else
        a < 0 && return 0.0
        (a <= 1 || b == 0) && return 1.0
    end
    if _is_classical(q)
        if what === :int
            return Float64(a)^b
        elseif what === :fact
            tab = classical_tables(Float64, a)
            m, e = _split_mul(1.0, 0, tab, a, b)
            return ldexp(m, e)
        else
            a <= 66 && return Float64(binomial(a, b))          # exact in Int64, rounded once
            tab = classical_tables(Float64, a)
            m, e = _split_mul(1.0, 0, tab, a, 1)
            m, e = _split_mul(m, e, tab, b, -1)
            m, e = _split_mul(m, e, tab, a - b, -1)
            return ldexp(m, e)
        end
    elseif k !== nothing
        kk = Int(k)
        what === :int && return _qint_level(a, kk)^b
        what === :fact && return _qfacts_level(kk, ((a, b),))
        return _qfacts_level(kk, ((a, 1), (b, -1), (a - b, -1)))
    elseif q isa Union{Float64,ComplexF64}
        if what === :int
            v = _qint_analytic(a, q)
            return v === nothing ? nothing : v^b
        elseif what === :fact
            f = _qfact_split(a, q)
            return f === nothing ? nothing : _split_value(f[1], f[2], b)
        else
            f = _qbinomial_split(a, b, q)
            return f === nothing ? nothing : _split_value(f[1], f[2], 1)
        end
    end
    return nothing
end

"The cyclotomic monomial of `Π [n]!^{c}` — what the numeric projections consume."
function _pairs_mono(pairs)
    pairs === nothing && return ZERO_MONOMIAL
    isempty(pairs) && return ONE_MONOMIAL
    buf = CycloBuffer(max(maximum(first, pairs), 1))
    for (m, c) in pairs
        add_qfact!(buf, Int(m), Int(c))
    end
    return snapshot(buf)
end

"`nothing` for a vanishing product, `Pair[]` for one that is identically 1, else the q-factorial content."
function _qint_pairs(n::Int, p::Int)
    n < 0 && return nothing
    (n <= 1 || p == 0) && return Pair{Int,Int}[]
    return [n => p, n - 1 => -p]
end
function _qfact_pairs(n::Int, p::Int)
    n < 0 && return nothing
    (n <= 1 || p == 0) && return Pair{Int,Int}[]
    return [n => p]
end
function _qbinomial_pairs(n::Int, m::Int)
    (m < 0 || m > n) && return nothing
    (m == 0 || m == n) && return Pair{Int,Int}[]
    return [n => 1, m => -1, n - m => -1]
end

"""
    qint(n, p = 1; k, q, exact, T)

The q-integer `[n]^p`, classical by default — `[n] = (qⁿ − q⁻ⁿ)/(q − q⁻¹)`, with the Turaev–Viro
convention `[0] = 1`. `n < 0` is zero.

The keywords are the ones every symbol takes: `k` for a level, `q` for a parameter, `exact = true` for
the exact classical rational or, with `k`, the real-basis level value. `qint(Symbolic(), n)` returns the
rule-backed symbolic x-form instead; its `.dcr` property constructs a one-term compatibility DCR.

```julia
qint(5)                 # 5.0
qint(5; k = 10)         # [5] at q = exp(iπ/12)
qint(5; q = 0.7)
qint(Exact(10), 5)      # exact, as a polynomial in x; `radical` for the surd
```
"""
function qint end

Base.@constprop :aggressive @inline function qint(n::Integer, p::Integer = 1; k = nothing, q = nothing, exact::Bool = false,
              T::Type{TT} = Float64, workspace = nothing) where {TT}
    q = _evaluation_q(k, q, exact)
    if !exact && T === Float64 && !(k isa AbstractVector)
        v = _qnumber_float(:int, Int(n), Int(p), k, q)
        v === nothing || return v
        # a Float64 target lands on a Float64 or, at complex q, a ComplexF64: saying so keeps the
        # return unboxed, since the general route below is not inferable
        return _product_value(_qint_pairs(Int(n), Int(p)), k, q, exact, T; workspace = workspace)::Union{Float64,ComplexF64}
    end
    return _product_value(_qint_pairs(Int(n), Int(p)), k, q, exact, T; workspace = workspace)
end

"""
    qfact(n, p = 1; k, q, exact, T)

The q-factorial `([n]!)^p = ([1][2]⋯[n])^p`, classical by default, with the same keywords as [`qint`](@ref).
`n < 0` is zero and `[0]! = [1]! = 1`.
"""
function qfact end

Base.@constprop :aggressive @inline function qfact(n::Integer, p::Integer = 1; k = nothing, q = nothing, exact::Bool = false,
               T::Type{TT} = Float64, workspace = nothing) where {TT}
    q = _evaluation_q(k, q, exact)
    if !exact && T === Float64 && !(k isa AbstractVector)
        v = _qnumber_float(:fact, Int(n), Int(p), k, q)
        v === nothing || return v
        # a Float64 target lands on a Float64 or, at complex q, a ComplexF64: saying so keeps the
        # return unboxed, since the general route below is not inferable
        return _product_value(_qfact_pairs(Int(n), Int(p)), k, q, exact, T; workspace = workspace)::Union{Float64,ComplexF64}
    end
    return _product_value(_qfact_pairs(Int(n), Int(p)), k, q, exact, T; workspace = workspace)
end

"""
    qbinomial(n, m; k, q, exact, T)

The Gaussian binomial `[n]! / ([m]! [n−m]!)`, classical by default, with the same keywords as
[`qint`](@ref). Outside `0 ≤ m ≤ n` it is zero.
"""
function qbinomial end

Base.@constprop :aggressive @inline function qbinomial(n::Integer, m::Integer; k = nothing, q = nothing, exact::Bool = false,
                   T::Type{TT} = Float64, workspace = nothing) where {TT}
    q = _evaluation_q(k, q, exact)
    if !exact && T === Float64 && !(k isa AbstractVector)
        v = _qnumber_float(:binomial, Int(n), Int(m), k, q)
        v === nothing || return v
        # a Float64 target lands on a Float64 or, at complex q, a ComplexF64: saying so keeps the
        # return unboxed, since the general route below is not inferable
        return _product_value(_qbinomial_pairs(Int(n), Int(m)), k, q, exact, T; workspace = workspace)::Union{Float64,ComplexF64}
    end
    return _product_value(_qbinomial_pairs(Int(n), Int(m)), k, q, exact, T; workspace = workspace)
end

"The rule itself, which is what `Symbolic()` asks for; the `Symbolic` methods live in `targets.jl`."
_product_symbolic(pairs) = SymbolicValue(pairs === nothing ? _ZERO_PRODUCT_RULE : _qfact_rule(pairs))
