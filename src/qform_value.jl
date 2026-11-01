# A form as a value: what `==`, `iszero` and arithmetic mean for what `Exact()` and `Exact(k)` return. Read through
# the plain-Julia kernel at a level and through the classical evaluator at q = 1.


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
    if !(isempty(f.powers) && isempty(g.powers) && iszero(f.qpow) && iszero(g.qpow))
        # complex values: different to the 256 bits they are certified to, or undecided
        setprecision(BigFloat, 320) do
            abs(numeric_value(f; bits = 320) - numeric_value(g; bits = 320)) > big(2.0)^-256
        end && return false
        throw(ArgumentError("deciding equality of different forms with powers of q at a level is not available"))
    end
    return _form_value(f) == _form_value(g)
end
Base.:(==)(f::QForm, c::Union{Integer,Rational}) = _form_value(f) == c
Base.:(==)(c::Union{Integer,Rational}, f::QForm) = f == c
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
