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

