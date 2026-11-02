module QRecouplingNemoExt

using Nemo, LRUCache
import QRecoupling
import QRecoupling: _qfact_exactx, psi_exponents, _halve_exponents, At,
    Exact,
    ExactX,
    ExactXSum,
    FactorialSum,
    GENERIC_MAX_CHARS,
    GENERIC_MAX_DEGREE,
    NoRadical,
    QFORM_SHOWN_TERMS,
    QForm,
    RadExpr,
    Symbolic,
    LevelXValue,
    XValue,
    level_x_value,
    _arg,
    _collapse,
    _factor_step,
    _fold_label,
    _q_sum_str,
    _qform_body,
    _same_level,
    _unwrap,
    doubled,
    exact_x,
    fsymbol,
    generic_sixj,
    generic_value,
    gsymbol,
    has_radical_form,
    is_empty_sum,
    level_admissible,
    numeric_value,
    prove_identity,
    q3j_factorial,
    q6j,
    qform,
    qint,
    radical,
    radical_form,
    radical_levels,
    radicand,
    sixj_sum,
    symbol_of,
    symbol_rule,
    to_subscript,
    to_superscript,
    xform,
    xpolynomial,
    xvalue

include("algebra/generic_x.jl")
include("algebra/exact_x.jl")
include("algebra/exact_radicals.jl")
include("algebra/exact_display.jl")
include("algebra/exact_arithmetic.jl")
include("algebra/symbolic_display.jl")

# Radicals and polynomials are requests on a form. They go through the expanded value, which only this backend holds.
function _expanded(f::QForm)
    f.k === nothing && throw(ArgumentError("a radical needs a level: use `Exact(k)`"))
    iszero(f) && return zero(ExactX, f.k)
    f.weight == (0, 0) || throw(ArgumentError("radicals are not available for forms with powers of q"))
    return exact_x(f.rule, f.k)
end
radical(f::QForm; kwargs...) = radical(_expanded(f); kwargs...)
radical_form(f::QForm; kwargs...) = radical_form(_expanded(f); kwargs...)
has_radical_form(f::QForm) = has_radical_form(_expanded(f))
xpolynomial(f::QForm) = xpolynomial(_expanded(f))
radicand(f::QForm) = radicand(_expanded(f))

level_x_value(f::QForm) = LevelXValue(_expanded(f))
xpolynomial(v::LevelXValue) = xpolynomial(v.value)
radicand(v::LevelXValue) = radicand(v.value)
numeric_value(v::LevelXValue; kwargs...) = numeric_value(v.value; kwargs...)
Base.float(v::LevelXValue) = float(v.value)
Base.iszero(v::LevelXValue) = iszero(v.value.p)
Base.:(==)(a::LevelXValue, b::LevelXValue) = a.value == b.value
Base.hash(v::LevelXValue, h::UInt) = hash(v.value, hash(:LevelXValue, h))

end
