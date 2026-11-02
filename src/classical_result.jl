# The exact value at q = 1: a sign and a rational under a root.

"""
    ClassicalResult
Represents σ * √(sq_val) where sq_val is a Rational{BigInt}.
Preserves precision for the classical Ponzano-Regge limit.
"""
struct ClassicalResult <: Real
    sign::Int
    sq_val::Rational{BigInt} 
end

# --- ClassicalResult: Base & Properties ---

Base.:(==)(a::ClassicalResult, b::ClassicalResult) = (a.sign == b.sign) && (a.sq_val == b.sq_val)
Base.isapprox(a::ClassicalResult, b::ClassicalResult; atol=1e-15) = abs(Float64(a) - Float64(b)) < atol

Base.zero(::Type{ClassicalResult}) = ClassicalResult(0, 0//1)
Base.zero(res::ClassicalResult)    = ClassicalResult(0, 0//1)
Base.iszero(res::ClassicalResult)  = res.sign == 0 || res.sq_val == 0//1

Base.one(::Type{ClassicalResult})  = ClassicalResult(1, 1//1)
Base.one(res::ClassicalResult)     = ClassicalResult(1, 1//1)
Base.isone(res::ClassicalResult)   = res.sign == 1 && res.sq_val == 1//1

Base.sign(res::ClassicalResult)    = res.sign
Base.abs(res::ClassicalResult)     = ClassicalResult(abs(res.sign), res.sq_val)


# --- ClassicalResult: math operations ---

Base.:-(a::ClassicalResult) = ClassicalResult(-a.sign, a.sq_val)
Base.inv(a::ClassicalResult) = iszero(a) ? throw(DivideError()) : ClassicalResult(a.sign, 1 // a.sq_val)

Base.:*(a::ClassicalResult, b::ClassicalResult) = ClassicalResult(a.sign * b.sign, a.sq_val * b.sq_val)
Base.:/(a::ClassicalResult, b::ClassicalResult) = a * inv(b)

function Base.:^(res::ClassicalResult, p::Int)
    p == 0 && return one(ClassicalResult)
    iszero(res) && return zero(ClassicalResult)
    # also handles negative powers (inverts the fraction)
    new_sign = iseven(p) ? 1 : res.sign
    return ClassicalResult(new_sign, res.sq_val^p)
end

function Base.:*(c::Union{Int, Rational}, res::ClassicalResult)
    iszero(c) && return zero(ClassicalResult)
    iszero(res) && return res
    return ClassicalResult(res.sign * sign(c), res.sq_val * c^2)
end
Base.:*(res::ClassicalResult, c::Union{Int, Rational}) = c * res

function Base.:/(res::ClassicalResult, c::Union{Int, Rational})
    iszero(c) && throw(DivideError())
    iszero(res) && return res
    return ClassicalResult(res.sign * sign(c), res.sq_val / c^2)
end

function _is_perfect_square(r::Rational)
    num, den = numerator(r), denominator(r)
    return isqrt(num)^2 == num && isqrt(den)^2 == den
end

function Base.:+(a::ClassicalResult, b::ClassicalResult)
    iszero(a) && return b
    iszero(b) && return a
    
    # merge exact perfect squares exactly
    if _is_perfect_square(a.sq_val) && _is_perfect_square(b.sq_val)
        val = Rational{BigInt}(a) + Rational{BigInt}(b)
        return ClassicalResult(sign(val), val^2)
    end
    
    # merge identical radicals exactly
    if a.sq_val == b.sq_val
        if a.sign == b.sign
            return ClassicalResult(a.sign, a.sq_val * 4) 
        else
            return zero(ClassicalResult)
        end
    end
    
    # fallback to BigFloat 
    return BigFloat(a) + BigFloat(b)
end

Base.:-(a::ClassicalResult, b::ClassicalResult) = a + (-b)

# --- ClassicalResult: mixed type fallbacks ---

function Base.:+(a::ClassicalResult, c::Union{Integer, Rational})
    iszero(c) && return a
    return a + ClassicalResult(sign(c), Rational{BigInt}(c^2))
end
Base.:+(c::Union{Integer, Rational}, a::ClassicalResult) = a + c
Base.:-(a::ClassicalResult, c::Union{Integer, Rational}) = a + (-c)
Base.:-(c::Union{Integer, Rational}, a::ClassicalResult) = c + (-a)

Base.:+(a::ClassicalResult, b::AbstractFloat) = BigFloat(a) + BigFloat(b)
Base.:+(a::AbstractFloat, b::ClassicalResult) = BigFloat(a) + BigFloat(b)
Base.:-(a::ClassicalResult, b::AbstractFloat) = BigFloat(a) - BigFloat(b)
Base.:-(a::AbstractFloat, b::ClassicalResult) = BigFloat(a) - BigFloat(b)

Base.:*(res::ClassicalResult, x::AbstractFloat) = BigFloat(res) * BigFloat(x)
Base.:*(x::AbstractFloat, res::ClassicalResult) = BigFloat(x) * BigFloat(res)
Base.:/(res::ClassicalResult, x::AbstractFloat) = BigFloat(res) / BigFloat(x)


# --- Conversions ---

Base.Float64(res::ClassicalResult) = res.sign * Float64(sqrt(BigFloat(res.sq_val)))
Base.BigFloat(res::ClassicalResult) = res.sign * sqrt(BigFloat(res.sq_val))
Base.AbstractFloat(res::ClassicalResult) = Float64(res)
Base.float(::Type{ClassicalResult}) = Float64
Base.big(res::ClassicalResult) = BigFloat(res)

function Base.Rational{T}(res::ClassicalResult) where {T <: Integer}
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    if s_num^2 == num && s_den^2 == den
        return res.sign * Rational{T}(T(s_num), T(s_den))
    else
        throw(InexactError(:Rational, Rational{T}, res))
    end
end

"""
A classical exact value compares with the number it is.

`q6j(1,1,1,1,1,1; exact = true) == 1//6` used to be a `promote` error: the type is `√(sq_val)` with a
sign, so equality against a rational is an exact question (is the rational's square the stored one?) and
`Rational(res)` is the value when the root is rational. Without these a documented result could only be
checked by converting to `Float64` first, which is the opposite of the point.
"""
Base.Rational(res::ClassicalResult) = Rational{BigInt}(res)

function Base.:(==)(a::ClassicalResult, b::Union{Integer,Rational})
    iszero(a) && return iszero(b)
    Int(sign(b)) == a.sign || return false
    return Rational{BigInt}(b)^2 == a.sq_val
end
Base.:(==)(b::Union{Integer,Rational}, a::ClassicalResult) = a == b

"A value that is rational hashes as that rational, since the two compare equal."
function Base.hash(a::ClassicalResult, h::UInt)
    iszero(a.sign) && return hash(0, h)
    n = numerator(a.sq_val); d = denominator(a.sq_val)
    rn = isqrt(n); rd = isqrt(d)
    rn * rn == n && rd * rd == d && return hash(a.sign * rn // rd, h)
    return hash((:ClassicalResult, a.sign, n, d), h)
end
Base.:(==)(a::ClassicalResult, b::AbstractFloat) = Float64(a) == b
Base.:(==)(b::AbstractFloat, a::ClassicalResult) = Float64(a) == b
Base.isapprox(a::ClassicalResult, b::Real; kwargs...) = isapprox(Float64(a), Float64(b); kwargs...)
Base.isapprox(a::Real, b::ClassicalResult; kwargs...) = isapprox(Float64(a), Float64(b); kwargs...)

# ---- REPL print ----

function Base.show(io::IO, res::ClassicalResult)
    iszero(res) && return print(io, "0//1")
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    if s_num^2 == num && s_den^2 == den
        print(io, res.sign < 0 ? "-" : "", Rational(s_num, s_den))
    else
        print(io, res.sign < 0 ? "-" : "", "√(", res.sq_val, ")")
    end
end
