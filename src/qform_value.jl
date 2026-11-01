# A form as a value: what `==`, `iszero` and arithmetic mean for what `Exact()` and `Exact(k)` return. Arithmetic goes
# through the plain-Julia kernel at a level, and a generic form has none. `==` expands nothing.


"What a form means in arithmetic: the exact value at its level."
function _form_value(f::QForm)
    f.k === nothing && throw(ArgumentError("a generic form has no arithmetic: use `Exact(k)` for a level, or `f(1)` for the exact value at q = 1"))
    f.weight == (0, 0) || f.rule === nothing ||
        throw(ArgumentError("exact arithmetic at a level is not available for forms with powers of q"))
    iszero(f) && return zero(LevelValue, f.k)
    return level_value(f.rule, f.k)
end

_same_form_level(f::QForm, g::QForm) = f.k == g.k ||
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
Base.inv(f::QForm) = inv(_form_value(f))
Base.:^(f::QForm, n::Integer) = _form_value(f)^n
for op in (:+, :-, :*, :/)
    @eval begin
        function Base.$op(f::QForm, g::QForm)
            _same_form_level(f, g)
            return $op(_form_value(f), _form_value(g))
        end
        Base.$op(f::QForm, c::Union{Integer,Rational}) = $op(_form_value(f), c)
        Base.$op(c::Union{Integer,Rational}, f::QForm) = $op(c, _form_value(f))
        Base.$op(f::QForm, v::Union{LevelValue,LevelValueSum,ClassicalResult}) = $op(_form_value(f), v)
        Base.$op(v::Union{LevelValue,LevelValueSum,ClassicalResult}, f::QForm) = $op(v, _form_value(f))
    end
end
