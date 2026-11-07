# A form as a value: what `==`, `iszero` and arithmetic mean for what `Exact()` and `Exact(k)` return. Arithmetic goes
# through the plain-Julia kernel at a level, except for forms that are one product, which stay forms. A generic form has
# only that arithmetic. `==` expands nothing.


"Whether a form has arithmetic: one at a level, without powers of q."
function _check_arithmetic(f::QForm)
    f.k === nothing && throw(ArgumentError("a generic form has arithmetic only as one product of q-integers: use `Exact(k)` for a level, or `f(1)` for the exact value at q = 1"))
    f.weight == (0, 0) || f.rule === nothing ||
        throw(ArgumentError("exact arithmetic at a level is not available for forms with powers of q"))
    return nothing
end

"What a form means in arithmetic: the exact value at its level."
function _form_value(f::QForm)
    _check_arithmetic(f)
    iszero(f) && return zero(LevelValue, f.k)
    return level_value(f.rule, f.k)
end

"The zero that holds at every q, as `zero(QForm)` gives it (matrix products and `sum` start from it): it adds at any level."
_any_zero(f::QForm) = f.k === nothing && f.rule === nothing
_form_value(f::QForm, k) = _any_zero(f) && k !== nothing ? zero(LevelValue, k) : _form_value(f)

_same_form_level(f::QForm, g::QForm) = f.k == g.k || _any_zero(f) || _any_zero(g) ||
    throw(ArgumentError("these forms are at different levels: $(f.k) and $(g.k)"))

"A rule up to the order of its factors. Two values with one key are one value."
_rule_key(f::QForm) = (f.weight, Int(f.rule.sign0), f.rule.sqrt_pre, f.rule.alternating, f.rule.zlo, f.rule.zhi,
                       sort!(collect(f.rule.pre)), sort!(collect(f.rule.fac); by = x -> (x.a, x.b, x.c)))

"""
Whether two forms are equal: at every q for generic forms, and at the level otherwise. Both read the squares of the values
from their rules, one residue first. Equal squares are proved at a level by places and, at every q, on the written forms.
"""
function Base.:(==)(f::QForm, g::QForm)
    _same_form_level(f, g)
    (f.rule === nothing || g.rule === nothing) && return iszero(f) && iszero(g)
    if f.k === nothing
        _generic_squares_differ(f, g) && return false
        (iszero(f) || iszero(g)) && return iszero(f) && iszero(g)
        _rule_key(f) == _rule_key(g) && return true
        return _generic_equal(_written(f), _written(g))
    end
    V = LevelSquare(f.rule, f.k, f.weight); W = LevelSquare(g.rule, f.k, g.weight)
    _first_place_differs(V, W, f.k) && return false
    (iszero(f) || iszero(g)) && return iszero(f) && iszero(g)
    _rule_key(f) == _rule_key(g) && return true
    _level_squares_equal(V, W, f.k) || return false
    return _same_root(f, numeric_value(g; bits = 64))
end
function Base.:(==)(f::QForm, c::Union{Integer,Rational})
    (iszero(c) || iszero(f)) && return iszero(f) && iszero(c)
    if f.k === nothing                                   # a number is a constant function of q
        _generic_constant_differs(f, c) && return false
        s = sign(c)
        return _generic_equal(_written(f), QWritten{Rational{BigInt}}(nothing, Int(s), abs(Rational{BigInt}(c)), Pair{Int,Int}[], big(1)//big(1),
                                                                      Pair{Int,Int}[], [(big(1)//big(1), Pair{Int,Int}[])]))
    end
    _level_squares_equal(LevelSquare(f.rule, f.k, f.weight), LevelSquare(FactorialSum(0:0), f.k, (0, 0), c), f.k) || return false
    return _same_root(f, BigFloat(c))
end
Base.:(==)(c::Union{Integer,Rational}, f::QForm) = f == c
Base.:(==)(f::QForm, x::AbstractFloat) = isfinite(x) && f == Rational{BigInt}(x)
Base.:(==)(x::AbstractFloat, f::QForm) = f == x

# A form and a number have no hash in common, so for `isequal` a form equals only a form at the same level with the same value.
Base.isequal(f::QForm, g::QForm) = f.k == g.k && f == g
Base.isequal(::QForm, ::Union{Number,LevelValue,LevelValueSum}) = false
Base.isequal(::Union{Number,LevelValue,LevelValueSum}, ::QForm) = false

"Equal values have equal squares, so one residue of the square hashes them without writing anything."
function Base.hash(f::QForm, h::UInt)
    f.rule === nothing && return hash((:QForm, f.k, zero(UInt64)), h)
    if f.k === nothing
        V = LevelSquare(f.rule, collect(f.rule.zlo:f.rule.zhi), f.weight..., big(1)//big(1))
        return hash((:QForm, nothing, _generic_square(V)), h)
    end
    V = LevelSquare(f.rule, f.k, f.weight); n = f.k + 2
    t = _places(_level_field(n), 1, _squares_top(V, V))
    return hash((:QForm, f.k, _place_square(t, V, 1, t.as[1], n)), h)
end

"Equal squares leave v = ±w, both nonzero, so the certified numerical value tells which."
function _same_root(f::QForm, w)
    v = numeric_value(f; bits = 64)
    return setprecision(BigFloat, 64) do
        abs(v - w) < abs(v + w)
    end
end

Base.:(==)(f::QForm, v::Union{LevelValue,LevelValueSum}) = _form_value(f) == v
Base.:(==)(v::Union{LevelValue,LevelValueSum}, f::QForm) = f == v

Base.:-(f::QForm) = f.rule === nothing ? f : QForm(f.k, negated(f.rule), f.weight, f.reduce, nothing)

# ---- forms that are one product ----
#
# A form that is one product of q-factorials, sign · √(∏[n]!^c), stays a form under products, quotients, powers and
# roots, since these only add or scale its exponents. Nothing is expanded, so √[550] at k = 1098 costs what [550] does,
# and a generic form has this arithmetic too. The prefactor is a vector, so products of any length share one compiled
# rule shape. A sum goes through the field at a level, as before.

"Exponents n => e merged and sorted, without zeros or the labels 0 and 1, whose factorials are one."
function _merge_exponents(c::Vector{Pair{Int,Int}})
    sort!(c; by = first)
    out = Pair{Int,Int}[]
    for (n, e) in c
        n <= 1 && continue
        if !isempty(out) && out[end].first == n
            out[end] = n => out[end].second + e
        else
            push!(out, n => e)
        end
    end
    return filter!(p -> !iszero(p.second), out)
end

"Label exponents ∏[n]^e as factorial exponents ∏[n]!^c, since [n] = [n]!/[n − 1]!."
_label_factorials(L) = _merge_exponents(Pair{Int,Int}[p for (n, e) in L for p in (n => e, n - 1 => -e)])

"Whether the rule of a form is one term with no powers of q, which is what `_product_exponents` reads."
_is_product(f::QForm) = f.rule !== nothing && f.rule.sign0 != 0 && f.rule.zlo == f.rule.zhi && f.weight == (0, 0)

"""
A form that is one product of q-factorials, as `(sign, c)` with value sign · √(∏[n]!^c), or `nothing`. At a level every
argument is below h, so each factorial is positive, and a form with a larger one goes through the field.
"""
_product_exponents(f::QForm) = _is_product(f) ? _product_exponents(f.rule, f.k) : nothing
# The field of a form has an abstract type, so the loops run behind this barrier, compiled for the shape of the rule.
function _product_exponents(s::FactorialSum, k)
    z = s.zlo; m = s.sqrt_pre ? 1 : 2
    c = sizehint!(Pair{Int,Int}[], length(s.pre) + length(s.fac))
    for (n, e) in s.pre
        push!(c, Int(n) => m * Int(e))
    end
    for t in s.fac
        push!(c, _arg(t, z) => 2 * Int(t.c))
    end
    c = _merge_exponents(c)
    k === nothing || all(p -> p.first < k + 2, c) || return nothing
    return Int(s.sign0) * (s.alternating && isodd(z) ? -1 : 1), c
end

"The form sign · √(∏[n]!^c), with no root when every exponent is even."
function _product_form(k, reduce::Bool, sign::Int, c::Vector{Pair{Int,Int}})
    root = any(p -> isodd(p.second), c)
    pre = Pair{Int32,Int32}[Int32(n) => Int32(root ? e : e ÷ 2) for (n, e) in c]
    return QForm(k, FactorialSum(pre, root, (), false, Int8(sign), 0, 0), (0, 0), reduce, nothing)
end

"f · g^p for p = ±1 when both are one product, or `nothing`."
function _product_of(f::QForm, g::QForm, p::Int)
    (_is_product(f) && _is_product(g)) || return nothing
    a = _product_exponents(f); a === nothing && return nothing
    b = _product_exponents(g); b === nothing && return nothing
    c = append!(copy(a[2]), (n => p * e for (n, e) in b[2]))
    return _product_form(f.k, f.reduce & g.reduce, a[1] * b[1], _merge_exponents(c))
end

function Base.inv(f::QForm)
    a = _product_exponents(f)
    a === nothing && return inv(_form_value(f))
    return _product_form(f.k, f.reduce, a[1], [n => -e for (n, e) in a[2]])
end
function Base.:^(f::QForm, n::Integer)
    a = _product_exponents(f)
    a === nothing && return _form_value(f)^n
    return _product_form(f.k, f.reduce, isodd(n) ? a[1] : 1, [m => Int(n) * e for (m, e) in a[2]])
end
Base.:^(f::QForm, p::Rational) = _rational_power(f, p)          # two methods, as Base has ^(::Number, ::Rational)
Base.:^(f::QForm, p::AbstractFloat) = _rational_power(f, p)

"""
    sqrt(f::QForm)

The root of a form that is one product, ±c ∏[n]^e with no sum and no root of its own, such as a quantum dimension or a
6j symbol written as a single product. It is a form, at a level and at generic q, and nothing is expanded:
`sqrt(qint(Exact(1098), 550))` is `√([550])`. Products, quotients and powers of such roots stay forms. A sum of
q-integers, or a root already present, has no exact root here, and says so. So does the result of arithmetic with a sum,
since it is held as a sum: take the root of each symbol, as in `sqrt(qdim(Exact(k), a)) * sqrt(qdim(Exact(k), b))`.
"""
function Base.sqrt(f::QForm)
    f.weight == (0, 0) || f.rule === nothing ||
        throw(ArgumentError("exact roots are not available for forms with powers of q: " * _LEVEL_ROOT_HINT))
    a = _product_exponents(f)
    if a !== nothing
        a[1] < 0 && throw(DomainError(f, "sqrt of a negative value"))
        all(p -> iseven(p.second), a[2]) && return _product_form(f.k, f.reduce, 1, [n => e ÷ 2 for (n, e) in a[2]])
        # A root of its own. At generic q the labels are independent, so this is not a square. At a level the
        # relations among the labels can make it one, and the written form shows that.
        f.k === nothing && throw(ArgumentError("only a form that is one product of q-integers has an exact root: " * _LEVEL_ROOT_HINT))
    elseif iszero(f)
        return f
    end
    return _written_root(f)
end

"The root read off the written form, for a form whose rule is a sum or whose radicand is a square only at its level."
function _written_root(f::QForm)
    w = _written(f)
    product = length(w.terms) == 1 && iszero(w.qpow) && isempty(w.powers) && isempty(w.rad) &&
              _is_perfect_square(w.radcoef)
    product || throw(ArgumentError("only a form that is one product of q-integers has an exact root: " * _LEVEL_ROOT_HINT))
    c, labels = only(w.terms)
    c *= w.coef * isqrt(numerator(w.radcoef)) // isqrt(denominator(w.radcoef))
    w.sign * sign(c) < 0 && throw(DomainError(f, "sqrt of a negative value"))
    L = _merge_exponents(collect(Pair{Int,Int}, Iterators.flatten((w.out, labels))))
    f.k === nothing || all(p -> p.first < f.k + 2, L) ||
        throw(ArgumentError("only a form that is one product of q-integers has an exact root: " * _LEVEL_ROOT_HINT))
    isone(abs(c)) && return _product_form(f.k, f.reduce, 1, _label_factorials(L))
    # A rational coefficient has no place in a rule, so the root is a value of the field.
    f.k === nothing && throw(ArgumentError("the root of a generic form is kept only for a product of q-integers alone: " * _LEVEL_ROOT_HINT))
    E = zeros(Int, (f.k + 2) ÷ 2)
    for (n, e) in L
        E[min(n, f.k + 2 - n)] += e
    end
    return _level_root(f.k, abs(Rational{BigInt}(c)), E)
end
for (op, p) in ((:*, 1), (:/, -1))
    @eval function Base.$op(f::QForm, g::QForm)
        _same_form_level(f, g)
        _any_zero(f) | _any_zero(g) && return _any_zero(f) & _any_zero(g) ? f : $op(_form_value(f, g.k), _form_value(g, f.k))
        h = _product_of(f, g, $p)
        return h === nothing ? $op(_form_value(f), _form_value(g)) : h
    end
end

const _NO_MIXED_ARITHMETIC = "a form and a value at q = 1 have no arithmetic together: a level form is a number at another q, and a formula has no number"
for op in (:+, :-)
    @eval function Base.$op(f::QForm, g::QForm)
        _same_form_level(f, g)
        _any_zero(f) & _any_zero(g) && return f
        return $op(_form_value(f, g.k), _form_value(g, f.k))
    end
end
for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(f::QForm, c::Union{Integer,Rational}) = $op(_form_value(f), c)
        Base.$op(c::Union{Integer,Rational}, f::QForm) = $op(c, _form_value(f))
        Base.$op(f::QForm, v::Union{LevelValue,LevelValueSum}) = $op(_form_value(f, v.k), v)
        Base.$op(v::Union{LevelValue,LevelValueSum}, f::QForm) = $op(v, _form_value(f, v.k))
        Base.$op(::QForm, ::_ClassicalExact) = throw(ArgumentError(_NO_MIXED_ARITHMETIC))
        Base.$op(::_ClassicalExact, ::QForm) = throw(ArgumentError(_NO_MIXED_ARITHMETIC))
    end
end

# ---- order: a value at a level is a real number ----
#
# Real unless it carries powers of q (then it may be complex) or is generic (then it is not one number). Order is
# exact: certified Float64 values separate almost every pair; otherwise exact `==`; otherwise, the values differ,
# so certified values to more bits must separate them. Against a float the float's exact value is used, as `==` does.

const _LevelOrdered = Union{QForm,LevelValue,LevelValueSum}
const _PlainReal = Union{Integer,Rational,AbstractFloat}

function _ordered(f::QForm)
    f.k !== nothing && f.weight == (0, 0) && return f
    throw(ArgumentError("this form " * (f.k === nothing ? "holds at every q" : "carries powers of q and may be complex") *
                        ": order, `real`, `imag` and `conj` need a real value at a level"))
end
_ordered(x) = x
_certified(f::QForm) = _form_value(f)                 # LevelValue: certified to any number of bits
_certified(x) = x
_bits_value(x::Union{LevelValue,LevelValueSum}, bits) = numeric_value(x; bits)
_bits_value(x::_PlainReal, bits) = setprecision(() -> BigFloat(x), BigFloat, bits)

"sign(a − b), exactly."
function _level_cmp(a, b)
    x, y = Float64(_ordered(a)), Float64(_ordered(b))  # each within a few units of its last place
    m = max(abs(x), abs(y))
    2.0^-960 < m < Inf && abs(x - y) > 8eps() * m && return cmp(x, y)
    a == b && return 0
    a, b = _certified(a), _certified(b)
    bits = 128
    while true
        x, y = _bits_value(a, bits), _bits_value(b, bits)
        d = setprecision(() -> abs(x - y) - (abs(x) + abs(y)) * BigFloat(2)^(4 - bits), BigFloat, bits + 16)
        d > 0 && return cmp(x, y)
        bits *= 2
    end
end

Base.:<(a::_LevelOrdered, b::_LevelOrdered) = _level_cmp(a, b) < 0
Base.:<(a::_LevelOrdered, b::_PlainReal) = !isnan(b) && _level_cmp(a, b) < 0
Base.:<(a::_PlainReal, b::_LevelOrdered) = !isnan(a) && _level_cmp(a, b) < 0
Base.isless(a::_LevelOrdered, b::Union{_LevelOrdered,Integer,Rational}) = a < b
Base.isless(a::Union{Integer,Rational}, b::_LevelOrdered) = a < b
Base.isless(a::QForm, b::AbstractFloat) = isnan(b) || a < b            # NaN last, as Base does for `Real`
Base.isless(a::AbstractFloat, b::QForm) = !isnan(a) && a < b
Base.:<=(a::_LevelOrdered, b::_LevelOrdered) = _level_cmp(a, b) <= 0
Base.sign(x::_LevelOrdered) = _level_cmp(x, 0)
Base.signbit(f::QForm) = f < 0
Base.abs(f::QForm) = signbit(f) ? -f : f
Base.isfinite(::_LevelOrdered) = true
Base.isinf(::_LevelOrdered) = false
Base.real(f::QForm) = _ordered(f)
Base.conj(f::QForm) = _ordered(f)
Base.imag(f::QForm) = zero(_ordered(f))
# The zero that holds at every q, for code that asks a type for its zero (Plots does); it mixes with no level.
Base.zero(::Type{<:_LevelOrdered}) = _qform_zero(nothing)
for (f, pick) in ((:max, :(isless(b, a) ? a : b)), (:min, :(isless(b, a) ? b : a)))   # a rational has no level to
    @eval Base.$f(a::Union{LevelValue,LevelValueSum}, b::_ExactRational) = $pick       # promote to, so Base's `Real`
    @eval Base.$f(a::_ExactRational, b::Union{LevelValue,LevelValueSum}) = $pick       # fallback cannot
end

Base.promote_rule(::Type{<:_LevelOrdered}, ::Type{T}) where {T<:AbstractFloat} = T
Base.promote_rule(::Type{<:_LevelOrdered}, ::Type{<:AbstractIrrational}) = Float64
Base.promote_rule(::Type{<:AbstractIrrational}, ::Type{<:_LevelOrdered}) = Float64   # Base's own order asks for zero(T)
Base.convert(::Type{T}, x::Union{LevelValue,LevelValueSum}) where {T<:Union{Base.IEEEFloat,BigFloat}} =
    T(numeric_value(x; bits = max(128, precision(T) + 64)))
Base.BigFloat(x::Union{LevelValue,LevelValueSum}) = numeric_value(x)
Base.big(x::Union{LevelValue,LevelValueSum}) = BigFloat(x)
# a matrix product of forms holds forms and the values their arithmetic gives
Base.promote_rule(::Type{QForm}, ::Type{T}) where {T<:Union{LevelValue,LevelValueSum}} = T
Base.promote_rule(::Type{LevelValue}, ::Type{LevelValueSum}) = LevelValueSum
Base.convert(::Type{LevelValue}, f::QForm) = _form_value(f)
Base.convert(::Type{LevelValueSum}, f::QForm) = LevelValueSum(_form_value(f))
Base.convert(::Type{LevelValueSum}, v::LevelValue) = LevelValueSum(v)
Base.isequal(::QForm, ::Union{LevelValue,LevelValueSum}) = false
Base.isequal(::Union{LevelValue,LevelValueSum}, ::QForm) = false
