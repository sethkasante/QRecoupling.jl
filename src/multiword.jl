# ---------------------------------------------------------------------------------
#  K-word floating point, for the precision tier beyond compensated Horner
#
#  A value is an unevaluated sum of K doubles (≈ 53K bits), built only from error-free transformations
#  (TwoSum, FMA-based TwoProd) and a renormalisation by VecSum passes, so it allocates nothing. The
#  summation kernel is generic in its number type, so the ratio loop run in `MW{K}` has relative error
#  ≈ u_K·κ, u_K ≈ 2^(−53K) — the accuracy class of K-fold compensation. K = 3 gives full double output up
#  to κ ≈ 1e30, K = 4 up to κ ≈ 1e46, 3.5–7× faster than the BigFloat fallback it precedes
#  (`dev/results/kfold_lazy_families.md`).
# ---------------------------------------------------------------------------------

struct MW{K} <: AbstractFloat
    w::NTuple{K,Float64}
end

@inline Base.getindex(a::MW, i::Int) = a.w[i]
MW{K}(x::Float64) where {K} = MW{K}(ntuple(i -> i == 1 ? x : 0.0, K))
MW{K}(x::Integer) where {K} = MW{K}(Float64(x))            # exact below 2^53
MW{K}(x::MW{K}) where {K} = x
function MW{K}(x::BigFloat) where {K}
    w = ntuple(_ -> 0.0, K)
    r = x
    for i in 1:K
        h = Float64(r)
        w = Base.setindex(w, h, i)
        r = r - h
    end
    return MW{K}(w)
end
Base.BigFloat(a::MW) = sum(BigFloat(x) for x in a.w)
Base.Float64(a::MW) = foldr(+, a.w)
Base.convert(::Type{MW{K}}, x::Real) where {K} = MW{K}(x)
Base.promote_rule(::Type{MW{K}}, ::Type{<:Union{Integer,Float64}}) where {K} = MW{K}

Base.zero(::Type{MW{K}}) where {K} = MW{K}(0.0)
Base.one(::Type{MW{K}}) where {K} = MW{K}(1.0)
Base.zero(::MW{K}) where {K} = zero(MW{K})
Base.one(::MW{K}) where {K} = one(MW{K})
Base.eps(::Type{MW{K}}) where {K} = 2.0^(-50K)          # conservative: renormalisation is not exactly rounded
Base.precision(::Type{MW{K}}) where {K} = 53K
Base.iszero(a::MW) = iszero(a[1])
Base.isfinite(a::MW) = isfinite(a[1])
Base.signbit(a::MW) = signbit(a[1])
Base.:-(a::MW{K}) where {K} = MW{K}(map(-, a.w))
Base.abs(a::MW) = signbit(a[1]) ? -a : a
Base.ldexp(a::MW{K}, e::Integer) where {K} = MW{K}(map(x -> ldexp(x, e), a.w))
function Base.frexp(a::MW{K}) where {K}
    _, e = frexp(a[1])
    return ldexp(a, -e), e
end
Base.:<(a::MW, b::MW) = (a - b)[1] < 0
Base.:<=(a::MW, b::MW) = (a - b)[1] <= 0
Base.:(==)(a::MW, b::MW) = iszero((a - b)[1])
Base.:<(a::MW, b::Real) = a < MW{length(a.w)}(Float64(b))
Base.:<(a::Real, b::MW) = MW{length(b.w)}(Float64(a)) < b
Base.:<=(a::MW, b::Real) = a <= MW{length(a.w)}(Float64(b))


"One VecSum pass: Σt = s + Σ errs exactly. Recursion on the tuple tail, so it unrolls at compile time."
@inline _vecsum(t::Tuple{Float64}) = (t[1], ())
@inline function _vecsum(t::NTuple{M,Float64}) where {M}
    s, errs = _vecsum(Base.tail(t))
    s2, e = _two_sum(t[1], s)
    return s2, (e, errs...)
end

@inline _sumall(t::Tuple{}) = 0.0
@inline _sumall(t::Tuple) = t[1] + _sumall(Base.tail(t))

"Extract K leading words from a list of terms by K VecSum passes (type-stable: each pass is one call)."
@inline _extract(::Val{0}, t::Tuple) = ()
@inline _extract(::Val{K}, t::Tuple{}) where {K} = ntuple(_ -> 0.0, Val(K))
@inline function _extract(::Val{1}, t::NTuple{M,Float64}) where {M}
    s, rest = _vecsum(t)
    return (s + _sumall(rest),)                     # leftover ≈ u^K relative folds into the last word
end
@inline function _extract(::Val{K}, t::NTuple{M,Float64}) where {K,M}
    s, rest = _vecsum(t)
    return (s, _extract(Val(K - 1), rest)...)
end
"Renormalise to K words with a well-formed leading word: a second, cheap pass over the K words fixes the
case where cancellation leaves the first word zero and the value in a later one."
@inline _renorm(::Val{K}, t::Tuple) where {K} = MW{K}(_extract(Val(K), _extract(Val(K), t)))

Base.:+(a::MW{K}, b::MW{K}) where {K} = _renorm(Val(K), (a.w..., b.w...))
Base.:-(a::MW{K}, b::MW{K}) where {K} = a + (-b)
Base.:+(a::MW{K}, b::Real) where {K} = a + MW{K}(Float64(b))
Base.:+(a::Real, b::MW{K}) where {K} = MW{K}(Float64(a)) + b
Base.:-(a::MW{K}, b::Real) where {K} = a - MW{K}(Float64(b))
Base.:-(a::Real, b::MW{K}) where {K} = MW{K}(Float64(a)) - b

function Base.:*(a::MW{3}, b::MW{3})
    p11, e11 = _two_prod(a[1], b[1])
    p12, e12 = _two_prod(a[1], b[2])
    p21, e21 = _two_prod(a[2], b[1])
    return _renorm(Val(3), (p11, p12, p21, e11, e12, e21, a[1] * b[3], a[2] * b[2], a[3] * b[1]))
end
function Base.:*(a::MW{4}, b::MW{4})
    p11, e11 = _two_prod(a[1], b[1])
    p12, e12 = _two_prod(a[1], b[2]); p21, e21 = _two_prod(a[2], b[1])
    p13, e13 = _two_prod(a[1], b[3]); p22, e22 = _two_prod(a[2], b[2]); p31, e31 = _two_prod(a[3], b[1])
    return _renorm(Val(4), (p11, p12, p21, e11, p13, p22, p31, e12, e21,
                            e13, e22, e31, a[1] * b[4], a[2] * b[3], a[3] * b[2], a[4] * b[1]))
end
Base.:*(a::MW{K}, b::Real) where {K} = a * MW{K}(Float64(b))
Base.:*(a::Real, b::MW{K}) where {K} = MW{K}(Float64(a)) * b

"Division by Newton correction: K rounds of q ← q + (a − q·b)/b₁."
function Base.:/(a::MW{K}, b::MW{K}) where {K}
    q = MW{K}(a[1] / b[1])
    for _ in 2:K
        r = a - q * b
        q = q + MW{K}(r[1] / b[1])
    end
    return q
end
Base.:/(a::MW{K}, b::Real) where {K} = a / MW{K}(Float64(b))
Base.:/(a::Real, b::MW{K}) where {K} = MW{K}(Float64(a)) / b
Base.inv(a::MW{K}) where {K} = one(MW{K}) / a

function Base.sqrt(a::MW{K}) where {K}
    iszero(a) && return a
    s = MW{K}(sqrt(a[1]))
    for _ in 2:K
        r = a - s * s
        s = s + MW{K}(r[1] / (2 * s[1]))
    end
    return s
end
Base.fma(a::MW{K}, b::MW{K}, c::MW{K}) where {K} = a * b + c
Base.:^(a::MW{K}, n::Integer) where {K} = n >= 0 ? prod(a for _ in 1:n; init = one(MW{K})) : inv(a^(-n))
Base.log10(a::MW) = log10(Float64(a))
Base.show(io::IO, a::MW{K}) where {K} = print(io, "MW{$K}(", BigFloat(a), ")")

_table_precision(::Type{MW{K}}) where {K} = 53K + 64
Base.:-(x::BigFloat, y::MW) = x - BigFloat(y)

const MW3_TABLES = LevelCache{QIntTables{MW{3}}}()
const MW4_TABLES = LevelCache{QIntTables{MW{4}}}()

"Level tables in K-word precision (K = 3 or 4), cached per level like the Float64 tables."
mw_tables(::Val{3}, k::Int) = get_level!(() -> _qint_tables_from(MW{3}, () -> _qints_wide(k + 2, k + 1)), MW3_TABLES, k)
mw_tables(::Val{4}, k::Int) = get_level!(() -> _qint_tables_from(MW{4}, () -> _qints_wide(k + 2, k + 1)), MW4_TABLES, k)
