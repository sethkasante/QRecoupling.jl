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
Base.:^(f::QForm, p::Union{Rational,AbstractFloat}) = _rational_power(f, p)

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
        h = _product_of(f, g, $p)
        return h === nothing ? $op(_form_value(f), _form_value(g)) : h
    end
end

const _NO_MIXED_ARITHMETIC = "a form and a value at q = 1 have no arithmetic together: a level form is a number at another q, and a formula has no number"
for op in (:+, :-)
    @eval function Base.$op(f::QForm, g::QForm)
        _same_form_level(f, g)
        return $op(_form_value(f), _form_value(g))
    end
end
for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(f::QForm, c::Union{Integer,Rational}) = $op(_form_value(f), c)
        Base.$op(c::Union{Integer,Rational}, f::QForm) = $op(c, _form_value(f))
        Base.$op(f::QForm, v::Union{LevelValue,LevelValueSum}) = $op(_form_value(f), v)
        Base.$op(v::Union{LevelValue,LevelValueSum}, f::QForm) = $op(v, _form_value(f))
        Base.$op(::QForm, ::_ClassicalExact) = throw(ArgumentError(_NO_MIXED_ARITHMETIC))
        Base.$op(::_ClassicalExact, ::QForm) = throw(ArgumentError(_NO_MIXED_ARITHMETIC))
    end
end
