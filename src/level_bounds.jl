# Small outward enclosures for exact-level certification. MPFR operations are correctly rounded;
# expanding each result by one adjacent representable number encloses the exact operation.
# No global rounding-mode changes are needed. This is not a public interval-arithmetic API.
_enclose(x::Integer) = iszero(x) ? (BigFloat(0), BigFloat(0)) :
    (prevfloat(BigFloat(x)), nextfloat(BigFloat(x)))
_iadd(a, b) = (prevfloat(a[1] + b[1]), nextfloat(a[2] + b[2]))
function _imul(a, b)
    p = (a[1] * b[1], a[1] * b[2], a[2] * b[1], a[2] * b[2])
    return prevfloat(minimum(p)), nextfloat(maximum(p))
end
function _idiv(a, b)
    b[1] > 0 || throw(DomainError(b, "the enclosing denominator must be positive"))
    return _imul(a, (prevfloat(inv(b[2])), nextfloat(inv(b[1]))))
end
_isqrt(a) = (prevfloat(sqrt(max(a[1], BigFloat(0)))), nextfloat(sqrt(a[2])))
_ilog2(a) = (prevfloat(log2(a[1])), nextfloat(log2(a[2])))

"Enclose sin(pi*n/h) for 0 < n <= h/2, using its monotonicity on this interval."
function _sin_interval(n::Int, h::Int)
    2n == h && return (BigFloat(1), BigFloat(1))
    a = _idiv(_enclose(n), _enclose(h))
    return prevfloat(sinpi(max(BigFloat(0), a[1]))),
           nextfloat(sinpi(min(BigFloat(0.5), a[2])))
end

"Enclose the absolute value of [n] at q^a, with a coprime to 2h and h not dividing n."
function _qint_interval(n::Int, a::Int, h::Int)
    n == 1 && return (BigFloat(1), BigFloat(1))
    r = mod(widemul(a, n), h)
    r = Int(min(r, h - r))
    r > 0 || throw(DomainError(n, "a vanishing q-integer has no nonzero enclosure"))
    b = mod(a, h); b = min(b, h - b)
    return _idiv(_sin_interval(r, h), _sin_interval(b, h))
end

"Fractional bits of the fixed-point logarithms in the prime-count bound."
const LOG_FRACTION_BITS = 20

"""
Integers `(up, lo)` with lo[r] ≤ 2^20 log2 sin(πr/h) ≤ up[r] for r ≤ h/2. A conjugate of a q-integer is a ratio of two
of these sines, so once they are enclosed the bound on a whole sum needs integer additions only.
"""
function _log_sine_table(h::Int)
    up = Vector{Int}(undef, h ÷ 2); lo = similar(up)
    setprecision(BigFloat, 128) do
        for r in 1:h÷2
            a = _ilog2(_sin_interval(r, h))
            up[r] = ceil(Int, ldexp(a[2], LOG_FRACTION_BITS))
            lo[r] = floor(Int, ldexp(a[1], LOG_FRACTION_BITS))
        end
    end
    return up, lo
end

function _interval_ball(a)
    mid = (a[1] + a[2]) / 2
    err = nextfloat(max(mid - a[1], a[2] - mid))
    return mid, err
end
