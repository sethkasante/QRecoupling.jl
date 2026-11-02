# ---------------------------------------------------------------------------------
#  Error-free transformations and double-word arithmetic
#
#  TwoProd/TwoSum and the double-word products, quotients and square roots used by compensated Horner
#  (`factorial_kernels.jl`), the K-word tiers (`multiword.jl`) and the family recurrences (`families.jl`).
# ---------------------------------------------------------------------------------

# ---- error-free transformations and double-word arithmetic ----

@inline _two_prod(a, b) = (p = a * b; (p, fma(a, b, -p)))
@inline function _two_sum(a, b)
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)
end
@inline _fast_two_sum(a, b) = (s = a + b; (s, b - (s - a)))

"(ah + al)(bh + bl) as a double word."
@inline function _dw_mul(ah, al, bh, bl)
    p, e = _two_prod(ah, bh)
    e = fma(ah, bl, fma(al, bh, e))
    return _fast_two_sum(p, e)
end

"√(h + l) as a double word, h > 0."
@inline function _dw_sqrt(h, l)
    r = sqrt(h)
    p, e = _two_prod(r, r)
    corr = ((h - p) - e + l) / (2r)
    return _fast_two_sum(r, corr)
end

"Multiply (mh + ml)·2^e by ([n]!)^c with double-word split entries."
@inline function _split_mul_dw(mh::T, ml::T, e::Int, tab::QIntTables{T}, n::Int, c::Integer) where {T}
    @inbounds if c == 1
        mh, ml = _dw_mul(mh, ml, tab.fm[n+1], tab.fml[n+1])
        return mh, ml, e + tab.fe[n+1]
    elseif c == -1
        mh, ml = _dw_mul(mh, ml, tab.gm[n+1], tab.gml[n+1])
        return mh, ml, e + tab.ge[n+1]
    elseif c > 0
        for _ in 1:c
            mh, ml = _dw_mul(mh, ml, tab.fm[n+1], tab.fml[n+1]); e += tab.fe[n+1]
        end
    else
        for _ in 1:-c
            mh, ml = _dw_mul(mh, ml, tab.gm[n+1], tab.gml[n+1]); e += tab.ge[n+1]
        end
    end
    return mh, ml, e
end

@inline function _renorm_dw(mh, ml, e::Int)
    fr, ex = frexp(mh)
    return fr, ldexp(ml, -ex), e + ex
end


# ---- tables of q-integers in double words ----
#
# Building [n] in `BigFloat` costs about a microsecond per entry and dominates any sweep over levels or q, where each
# point needs a table of its own. Here a wide number is used only for about 2√N seeds. Every [n] is then a few double-word
# operations away from them, so its error does not grow with n: no recurrence runs through the table.

"Below this many entries a wide build is cheap, and its low parts are a fraction of u² rather than a few u²."
const WIDE_TABLE_BELOW = 64

_word(x) = (h = Float64(x); (h, Float64(x - h)))

"s(i·t) and c(i·t) for i = 0, 1, …, as double words, stepping by the addition formula from (s(t), c(t))."
function _fill_angles!(sv, cv, st::BigFloat, ct::BigFloat, σ::Int)
    s, c = zero(BigFloat), one(BigFloat)
    @inbounds for i in eachindex(sv)
        sv[i] = _word(s); cv[i] = _word(c)
        s, c = s * ct + c * st, c * ct + σ * s * st
    end
    return nothing
end

"""
    _qints_by_addition(N, σ, sc) -> (q, ql, qi, qil)

[n] = s(n)/s(1) and its inverse for n in 1:N, as double words, with `sc(m)` the pair (s(m), c(m)) in `BigFloat`:
sine and cosine with σ = -1, or the hyperbolic pair with σ = +1, up to an argument where both stay nonnegative. Writing
n = aM + b, s(n) = s(aM)c(b) + c(aM)s(b) adds two nonnegative products, so nothing cancels, and every entry is a few u²
from the exact value whatever its n (`_read_errors`).
"""
function _qints_by_addition(N::Int, σ::Int, sc::F) where {F}
    M = isqrt(N) + 1
    A = N ÷ M + 1
    sb = Vector{NTuple{2,Float64}}(undef, M); cb = similar(sb)
    sa = Vector{NTuple{2,Float64}}(undef, A); ca = similar(sa)
    inv1 = setprecision(BigFloat, 256) do
        s1, c1 = sc(1)
        sM, cM = sc(M)
        _fill_angles!(sb, cb, s1, c1, σ)
        _fill_angles!(sa, ca, sM, cM, σ)
        _word(inv(s1))
    end
    q = Vector{Float64}(undef, N); ql = similar(q); qi = similar(q); qil = similar(q)
    q[1] = qi[1] = 1.0; ql[1] = qil[1] = 0.0
    @inbounds for n in 2:N
        a, b = divrem(n, M)
        s = iszero(a) ? sb[b+1] : _dwa(_dwm(sa[a+1], cb[b+1]), _dwm(ca[a+1], sb[b+1]))
        v = _dwm(s, inv1)
        w = _dwdiv((1.0, 0.0), v)
        q[n], ql[n] = v
        qi[n], qil[n] = w
    end
    return q, ql, qi, qil
end

"Running products of (q, ql) as split double words, from the empty product. The error grows with n (`_read_errors`)."
function _split_factorials_dw(q::Vector{Float64}, ql::Vector{Float64})
    N = length(q)
    fm = Vector{Float64}(undef, N + 1); fml = similar(fm); fe = Vector{Int}(undef, N + 1)
    h, l, e = 1.0, 0.0, 0
    @inbounds for n in 0:N
        n > 0 && ((h, l) = _dw_mul(h, l, q[n], ql[n]))
        h, d = frexp(h)
        l = ldexp(l, -d); e += d
        fm[n+1], fml[n+1], fe[n+1] = h, l, e
    end
    return fm, fml, fe
end

"The level tables without a wide number per entry: [h - n] = [n], so only the first quadrant is computed."
function QIntTables(::Type{Float64}, k::Int)
    h = k + 2
    N = k + 1
    N < WIDE_TABLE_BELOW && return _qint_tables_from(Float64, () -> _qints_wide(h, N))
    half = h ÷ 2
    q1, ql1, qi1, qil1 = _qints_by_addition(half, -1, m -> sincos(m * big(pi) / h))
    mirror(v) = (w = Vector{Float64}(undef, N); copyto!(w, v); for n in half+1:N; w[n] = v[h-n]; end; w)
    q, ql, qi, qil = mirror(q1), mirror(ql1), mirror(qi1), mirror(qil1)
    fm, fml, fe = _split_factorials_dw(q, ql)
    gm, gml, ge = _split_factorials_dw(qi, qil)
    return QIntTables{Float64}(q, qi, ql, qil, fm, fe, fml, gm, ge, gml, false)
end

"""
Errors of a table read in double words, in units of u²: a q-integer or its inverse, and a factorial or its inverse.
A wide build is covered by the constants of the bounds. In a table by addition the largest error measured on an entry
is 5.1 u², over every n up to k = 10⁵, and on a factorial about n/3 u². The analysis allows 26 u² and 31n u², all
aligned, which they do not do: these charges already add every read's error linearly along a sum.
"""
@inline _read_errors(tab::QIntTables) =
    tab.classical || length(tab.q) < WIDE_TABLE_BELOW ? (0, 0) : (6, 2 * length(tab.q))
