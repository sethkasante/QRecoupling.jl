# A form as a value: what `==`, `iszero` and arithmetic mean for what `Exact()` and `Exact(k)` return. Arithmetic goes
# through the plain-Julia kernel at a level and the classical evaluator at q = 1. `==` at a level expands nothing.


"What a form means in arithmetic: its value at q = 1, or the exact value at its level."
function _form_value(f::QForm)
    f.k === nothing && return classical(f)
    (isempty(f.powers) && iszero(f.qpow)) ||
        throw(ArgumentError("exact arithmetic at a level is not available for forms with powers of q"))
    iszero(f) && return zero(LevelValue, f.k)
    f.rule === nothing && throw(ArgumentError("this form holds no rule to compute with"))
    return level_value(f.rule, f.k)
end

_same_form_level(f::QForm, g::QForm) = f.k == g.k ||
    throw(ArgumentError("these forms are at different levels: $(f.k) and $(g.k)"))

function Base.:(==)(f::QForm, g::QForm)
    _same_form_level(f, g)
    (f.sign, f.coef, f.out, f.radcoef, f.rad, f.terms, f.qpow, f.powers) ==
        (g.sign, g.coef, g.out, g.radcoef, g.rad, g.terms, g.qpow, g.powers) && return true
    f.k === nothing && return _generic_equal(f, g)
    (iszero(f) || iszero(g)) && return iszero(f) && iszero(g)
    _level_squares_equal(_form_square(f), _form_square(g), f.k) || return false
    return _same_root(f, numeric_value(g; bits = 64))
end
function Base.:(==)(f::QForm, c::Union{Integer,Rational})
    f.k === nothing && return _form_value(f) == c
    (iszero(f) || iszero(c)) && return iszero(f) && iszero(c)
    _level_squares_equal(_form_square(f), LevelSquare(FactorialSum(0:0), f.k, (0, 0), c), f.k) || return false
    return _same_root(f, BigFloat(c))
end
Base.:(==)(c::Union{Integer,Rational}, f::QForm) = f == c

"""
A form at a level as a rule to read at places. `==` compares squares there, which lie in the field of q and need no
radicals, so a form decides through its rule and never through what it prints.
"""
function _form_square(f::QForm)
    f.rule === nothing && throw(ArgumentError("this form holds no rule to compute with"))
    return LevelSquare(f.rule, f.k, f.weight)
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

Base.:-(f::QForm) = typeof(f)(f.k, -f.sign, f.coef, f.out, f.radcoef, f.rad, f.terms,
                              f.rule === nothing ? nothing : negated(f.rule), f.weight, f.qpow, f.powers)
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
