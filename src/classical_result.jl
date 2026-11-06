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

"x^n for n ≥ 1 by repeated squaring, with no product by one."
function _square_power(x, n::Integer)
    while iseven(n)
        x = x * x; n >>= 1
    end
    out = x
    while (n >>= 1) > 0
        x = x * x
        isodd(n) && (out = out * x)
    end
    return out
end

"""
x^p for a half-integer p, through `sqrt(x)`. Exact values are a field element times one square root, so a root of a root
or a cube root is not one of them, and says so.
"""
function _rational_power(x, p::Rational)
    d = denominator(p); n = Int(numerator(p))
    d == 1 && return x^n
    d == 2 && return sqrt(x)^n
    throw(ArgumentError("the power $p of an exact value is not exact in this form, since only square roots are kept: " *
                        "`float(x)^p` gives its number"))
end
function _rational_power(x, p::AbstractFloat)
    (isfinite(p) && isinteger(2p)) || throw(ArgumentError("the power $p of an exact value is not exact: only integer " *
                                                          "and half-integer powers are, and `float(x)^p` gives its number"))
    return _rational_power(x, Int(2p) // 2)
end

"""
    sqrt(x::ClassicalResult)

Exact when the root is again ±√(rational), so when `x` is rational, such as the square of a symbol. A root of an
irrational √r is a fourth root, which is not kept, and says so.
"""
function Base.sqrt(a::ClassicalResult)
    iszero(a) && return zero(ClassicalResult)
    a.sign < 0 && throw(DomainError(a, "sqrt of a negative value"))
    _is_perfect_square(a.sq_val) || throw(ArgumentError("the square root of $a is a fourth root, which an exact value " *
                                                        "does not keep: `sqrt(float(x))` gives its number"))
    return ClassicalResult(1, isqrt(numerator(a.sq_val)) // isqrt(denominator(a.sq_val)))
end
Base.:^(a::ClassicalResult, p::Rational) = _rational_power(a, p)          # two methods, as Base has ^(::Number, ::Rational)
Base.:^(a::ClassicalResult, p::AbstractFloat) = _rational_power(a, p)

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
    
    # unlike roots stay exact as a sum of them
    return _classical_value(_add_terms(_root_terms(a), _root_terms(b)))
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

# ---- sums of unlike roots ----

"""
    ClassicalSum

An exact sum of unlike square roots at q = 1, Σ c_r √r with each r a square-free integer and each c_r rational, such as
`1//2 + √15/6`. It is what adding two `ClassicalResult`s gives when their roots differ, and it holds at least two roots: a
sum that collapses to one root, or to zero, is a `ClassicalResult` again.

Square roots of distinct square-free integers are linearly independent over ℚ, so this form is unique. `iszero` and `==`
compare coefficients and are proofs. `+`, `-`, `*`, `/`, `inv` and integer powers stay exact, and `float` gives a number.
"""
struct ClassicalSum <: Real
    terms::Dict{BigInt,Rational{BigInt}}
    ClassicalSum(terms::Dict{BigInt,Rational{BigInt}}) = new(terms)
end

const _SMALL_PRIMES = let sieve = trues(1 << 16)
    sieve[1] = false
    for i in 2:isqrt(length(sieve))
        sieve[i] && (sieve[i*i:i:end] .= false)
    end
    findall(sieve)
end

"""
`n = o²·r` with `r` square-free, by trial division. A value at q = 1 is built from factorials, so its primes are at most its
largest factorial argument, and what is left after these primes is a single prime unless an argument passes 2¹⁶.
"""
function _square_free(n::BigInt)
    o = big(1); r = big(1); m = n
    for p in _SMALL_PRIMES
        p * p > m && break
        e = 0
        while iszero(m % p)
            m ÷= p; e += 1
        end
        iszero(e) && continue
        o *= big(p)^(e ÷ 2)
        isodd(e) && (r *= p)
    end
    if m > 1
        s = isqrt(m)
        s * s == m ? (o *= s) : (r *= m)
    end
    return o, r
end

"σ√(n/d) as `r => c`: with n and d coprime their square-free parts are too, so √(n/d) = (o_n/(o_d r_d)) √(r_n r_d)."
function _root_terms(a::ClassicalResult)
    iszero(a) && return Dict{BigInt,Rational{BigInt}}()
    on, rn = _square_free(numerator(a.sq_val))
    od, rd = _square_free(denominator(a.sq_val))
    return Dict(rn * rd => a.sign * (on // (od * rd)))
end
_root_terms(a::ClassicalSum) = a.terms
_root_terms(c::Union{Integer,Rational}) =
    iszero(c) ? Dict{BigInt,Rational{BigInt}}() : Dict(big(1) => Rational{BigInt}(c))

"The simplest type for these terms: a `ClassicalResult` for at most one root."
function _classical_value(t::Dict{BigInt,Rational{BigInt}})
    length(t) >= 2 && return ClassicalSum(t)
    isempty(t) && return zero(ClassicalResult)
    r, c = only(t)
    return ClassicalResult(Int(sign(c)), c^2 * r)
end

function _add_terms(a::Dict{BigInt,Rational{BigInt}}, b::Dict{BigInt,Rational{BigInt}})
    t = copy(a)
    for (r, c) in b
        v = get(t, r, 0 // big(1)) + c
        iszero(v) ? delete!(t, r) : (t[r] = v)
    end
    return t
end

function _multiply_terms(a::Dict{BigInt,Rational{BigInt}}, b::Dict{BigInt,Rational{BigInt}})
    t = Dict{BigInt,Rational{BigInt}}()
    for (r, c) in a, (s, e) in b
        g = gcd(r, s)                                        # √r √s = g √(rs/g²), and rs/g² is square-free
        k = (r ÷ g) * (s ÷ g)
        v = get(t, k, 0 // big(1)) + c * e * g
        iszero(v) ? delete!(t, k) : (t[k] = v)
    end
    return t
end

const _ClassicalExact = Union{ClassicalResult,ClassicalSum}

Base.:+(a::ClassicalSum, b::ClassicalSum) = _classical_value(_add_terms(a.terms, b.terms))
Base.:+(a::ClassicalSum, b::ClassicalResult) = _classical_value(_add_terms(a.terms, _root_terms(b)))
Base.:+(a::ClassicalResult, b::ClassicalSum) = b + a
Base.:+(a::ClassicalSum, c::Union{Integer,Rational}) = _classical_value(_add_terms(a.terms, _root_terms(c)))
Base.:+(c::Union{Integer,Rational}, a::ClassicalSum) = a + c
Base.:-(a::ClassicalSum) = ClassicalSum(Dict(r => -c for (r, c) in a.terms))
Base.:-(a::ClassicalSum, b::ClassicalSum) = a + (-b)
Base.:-(a::ClassicalSum, b::ClassicalResult) = a + (-b)
Base.:-(a::ClassicalResult, b::ClassicalSum) = a + (-b)
Base.:-(a::ClassicalSum, c::Union{Integer,Rational}) = a + (-c)
Base.:-(c::Union{Integer,Rational}, a::ClassicalSum) = c + (-a)

Base.:*(a::ClassicalSum, b::ClassicalSum) = _classical_value(_multiply_terms(a.terms, b.terms))
Base.:*(a::ClassicalSum, b::ClassicalResult) = _classical_value(_multiply_terms(a.terms, _root_terms(b)))
Base.:*(a::ClassicalResult, b::ClassicalSum) = b * a
Base.:*(c::Union{Int,BigInt,Rational}, a::ClassicalSum) = _classical_value(_multiply_terms(a.terms, _root_terms(c)))
Base.:*(a::ClassicalSum, c::Union{Int,BigInt,Rational}) = c * a

"""
1/a by conjugates. Flipping the sign of √p for one prime p is a field automorphism, so a·σ(a) is nonzero, and it has no p in
its roots: each step removes a prime, until one root is left to invert.
"""
function Base.inv(a::ClassicalSum)
    r = maximum(keys(a.terms))
    i = findfirst(p -> iszero(r % p), _SMALL_PRIMES)
    p = i === nothing ? r : big(_SMALL_PRIMES[i])
    conjugate = ClassicalSum(Dict(s => iszero(s % p) ? -c : c for (s, c) in a.terms))
    return conjugate * inv(a * conjugate)
end
Base.:/(a::_ClassicalExact, b::ClassicalSum) = a * inv(b)
Base.:/(c::Union{Integer,Rational}, b::ClassicalSum) = c * inv(b)
Base.:/(a::ClassicalSum, b::ClassicalResult) = a * inv(b)
function Base.:/(a::ClassicalSum, c::Union{Int,BigInt,Rational})
    iszero(c) && throw(DivideError())
    return a * (1 // big(c))
end
function Base.:^(a::ClassicalSum, p::Integer)
    p < 0 && return inv(a)^(-p)
    return iszero(p) ? one(ClassicalResult) : _square_power(a, p)
end
Base.sqrt(a::ClassicalSum) = throw(ArgumentError("the square root of a sum of unlike roots is a nested root, which an " *
                                                 "exact value does not keep: `sqrt(float(x))` gives its number"))
Base.:^(a::ClassicalSum, p::Rational) = _rational_power(a, p)          # two methods, as Base has ^(::Number, ::Rational)
Base.:^(a::ClassicalSum, p::AbstractFloat) = _rational_power(a, p)

Base.zero(::Type{ClassicalSum}) = zero(ClassicalResult)
Base.zero(::ClassicalSum) = zero(ClassicalResult)
Base.one(::Type{ClassicalSum}) = one(ClassicalResult)
Base.one(::ClassicalSum) = one(ClassicalResult)
Base.iszero(::ClassicalSum) = false                          # it holds at least two independent roots
Base.isone(::ClassicalSum) = false

Base.:(==)(a::ClassicalSum, b::ClassicalSum) = a.terms == b.terms
Base.:(==)(::ClassicalSum, ::ClassicalResult) = false
Base.:(==)(::ClassicalResult, ::ClassicalSum) = false
Base.:(==)(::ClassicalSum, ::Union{Integer,Rational}) = false
Base.:(==)(::Union{Integer,Rational}, ::ClassicalSum) = false
Base.:(==)(a::ClassicalSum, b::AbstractFloat) = Float64(a) == b
Base.:(==)(b::AbstractFloat, a::ClassicalSum) = Float64(a) == b
Base.hash(a::ClassicalSum, h::UInt) = hash((:ClassicalSum, sort!(collect(a.terms); by = first)), h)

"The sum at this many bits or more: enough that its terms' cancellation leaves the digits asked for. It is never zero."
function _sum_value(a::ClassicalSum, bits::Int)
    work = bits + 32
    while true
        v, mass = setprecision(BigFloat, work) do
            total = big(0.0); size = big(0.0)
            for (r, c) in a.terms
                t = BigFloat(numerator(c)) / BigFloat(denominator(c)) * sqrt(BigFloat(r))
                total += t; size += abs(t)
            end
            total, size
        end
        if !iszero(v) && exponent(mass) - exponent(v) + bits + 8 <= work
            return setprecision(() -> BigFloat(v), BigFloat, bits)
        end
        work *= 2
    end
end

Base.Float64(a::ClassicalSum) = Float64(_sum_value(a, 64))
Base.BigFloat(a::ClassicalSum) = _sum_value(a, precision(BigFloat))
Base.AbstractFloat(a::ClassicalSum) = Float64(a)
Base.float(::Type{ClassicalSum}) = Float64
Base.big(a::ClassicalSum) = BigFloat(a)
Base.sign(a::ClassicalSum) = Int(sign(_sum_value(a, 64)))
Base.abs(a::ClassicalSum) = sign(a) < 0 ? -a : a

for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(a::ClassicalSum, b::AbstractFloat) = $op(BigFloat(a), BigFloat(b))
        Base.$op(a::AbstractFloat, b::ClassicalSum) = $op(BigFloat(a), BigFloat(b))
    end
end
Base.isapprox(a::ClassicalSum, b::ClassicalSum; kwargs...) = isapprox(Float64(a), Float64(b); kwargs...)
Base.isapprox(a::ClassicalSum, b::ClassicalResult; kwargs...) = isapprox(Float64(a), Float64(b); kwargs...)
Base.isapprox(a::ClassicalResult, b::ClassicalSum; kwargs...) = isapprox(Float64(a), Float64(b); kwargs...)
Base.isapprox(a::ClassicalSum, b::Real; kwargs...) = isapprox(Float64(a), Float64(b); kwargs...)
Base.isapprox(a::Real, b::ClassicalSum; kwargs...) = isapprox(Float64(a), Float64(b); kwargs...)

# The REPL shows the first terms, as for a form. `print` and `string` give them all.
function Base.show(io::IO, a::ClassicalSum)
    terms = sort!(collect(a.terms); by = first)
    most = get(io, :limit, false) ? QFORM_SHOWN_TERMS : length(terms)
    for (i, (r, c)) in enumerate(first(terms, most))
        i > 1 && print(io, c < 0 ? " - " : " + ")
        i == 1 && c < 0 && print(io, "-")
        n, d = abs(numerator(c)), denominator(c)
        if isone(r)
            print(io, n // d)
        else
            print(io, isone(n) ? "" : string(n), "√", r, isone(d) ? "" : "/" * string(d))
        end
    end
    most < length(terms) && print(io, " + … (", length(terms), " terms)")
end

# ---- order: exact values at q = 1 are real numbers ----
#
# Exact among exact values and rationals. Against a float the exact value is rounded to `Float64` first, as `==` does,
# so that `==`, `<` and `isless` agree. With a float or an irrational it promotes to `BigFloat`, as its float arithmetic
# gives; with a rational a `ClassicalResult` stays exact. Base's `Real` fallbacks then supply the rest (`isless`,
# `isinf`, `signbit`, `max`, sorting).

const _ExactRational = Union{Integer,Rational}
_exact(c::_ExactRational) = iszero(c) ? zero(ClassicalResult) : ClassicalResult(Int(sign(c)), Rational{BigInt}(c)^2)
_cmp_exact(a::ClassicalResult, b::ClassicalResult) =
    a.sign != b.sign ? cmp(a.sign, b.sign) : a.sign * cmp(a.sq_val, b.sq_val)
_cmp_exact(a, b) = Int(sign(a - b))                         # the difference is exact

for op in (:<, :<=)
    @eval Base.$op(a::_ClassicalExact, b::Union{_ClassicalExact,_ExactRational}) = $op(_cmp_exact(a, b), 0)
    @eval Base.$op(a::_ExactRational, b::_ClassicalExact) = $op(_cmp_exact(a, b), 0)
    @eval Base.$op(a::_ClassicalExact, b::AbstractFloat) = $op(Float64(a), b)
    @eval Base.$op(a::AbstractFloat, b::_ClassicalExact) = $op(a, Float64(b))
end
Base.isfinite(::_ClassicalExact) = true
Base.promote_rule(::Type{ClassicalResult}, ::Type{<:_ExactRational}) = ClassicalResult
Base.convert(::Type{ClassicalResult}, c::_ExactRational) = _exact(c)
Base.promote_rule(::Type{<:_ClassicalExact}, ::Type{<:Union{AbstractFloat,AbstractIrrational}}) = BigFloat
# a sum holds two roots at least, so it has no common type with a rational for Base's `max` to promote to
Base.max(a::ClassicalSum, b::Union{ClassicalResult,_ExactRational}) = a < b ? b : a
Base.max(a::Union{ClassicalResult,_ExactRational}, b::ClassicalSum) = a < b ? b : a
Base.min(a::ClassicalSum, b::Union{ClassicalResult,_ExactRational}) = b < a ? b : a
Base.min(a::Union{ClassicalResult,_ExactRational}, b::ClassicalSum) = b < a ? b : a

# Rounding is exact: one root by `isqrt` (an irrational root is never an integer), a sum by certified digits.
Base.isinteger(a::ClassicalResult) = isone(denominator(a.sq_val)) && _is_perfect_square(a.sq_val)
Base.isinteger(::ClassicalSum) = false
function _floor_exact(a::ClassicalResult)
    _is_perfect_square(a.sq_val) && return floor(BigInt, Rational{BigInt}(a))
    f = isqrt(numerator(a.sq_val) ÷ denominator(a.sq_val))   # ⌊√x⌋ = isqrt(⌊x⌋)
    return a.sign > 0 ? f : -f - 1
end
function _floor_exact(a::ClassicalSum)
    for bits in (64 .<< (0:30))
        f = setprecision(BigFloat, bits + 16) do
            v = BigFloat(_sum_value(a, bits)); e = abs(v) * BigFloat(2)^(2 - bits)
            floor(BigInt, v - e) == floor(BigInt, v + e) ? floor(BigInt, v) : nothing
        end
        f === nothing || return f
    end
end
function _round_exact(a::_ClassicalExact, r::RoundingMode)
    a isa ClassicalResult && _is_perfect_square(a.sq_val) && return round(BigInt, Rational{BigInt}(a), r)
    f = _floor_exact(r in (RoundNearest, RoundNearestTiesAway, RoundNearestTiesUp) ? a + 1//2 : a)   # no ties
    return r === RoundUp || (r === RoundToZero && a < 0) ? f + 1 : f
end
Base.round(a::_ClassicalExact, r::RoundingMode = RoundNearest) = _exact(_round_exact(a, r))
Base.round(::Type{T}, a::_ClassicalExact, r::RoundingMode = RoundNearest) where {T<:Integer} = T(_round_exact(a, r))
for (f, r) in ((:floor, RoundDown), (:ceil, RoundUp), (:trunc, RoundToZero))    # Julia 1.10 does not route these
    @eval Base.$f(::Type{T}, a::_ClassicalExact) where {T<:Integer} = round(T, a, $r)   # to `round(T, x, r)`
end
