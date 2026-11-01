module QRecouplingNemoExt

using Nemo, LRUCache
import QRecoupling
import QRecoupling: _qfact_exactx, psi_exponents, _halve_exponents, At,
    CycloBuffer,
    CyclotomicMonomial,
    DCR,
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
    SymbolicValue,
    XValue,
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
    mul!,
    numeric_value,
    phi_form,
    prove_identity,
    q3j_factorial,
    q6j,
    qform,
    qint,
    radical,
    radical_form,
    radical_levels,
    radicand,
    reset!,
    sixj_sum,
    snapshot,
    splits_completely,
    symbol_of,
    symbol_rule,
    to_subscript,
    to_superscript,
    x_form,
    xpolynomial,
    xvalue

include("algebra/generic_x.jl")
include("algebra/phi_form.jl")
include("algebra/exact_x.jl")
include("algebra/exact_radicals.jl")
include("algebra/exact_display.jl")
include("algebra/exact_arithmetic.jl")
include("algebra/symbolic_display.jl")

# Radicals and polynomials are requests on a form. They go through the expanded value, which only this backend holds.
function _expanded(f::QForm)
    f.k === nothing && throw(ArgumentError("a radical needs a level: use `Exact(k)`"))
    (isempty(f.powers) && iszero(f.qpow)) || throw(ArgumentError("radicals are not available for forms with powers of q"))
    iszero(f) && return zero(ExactX, f.k)
    f.rule === nothing && throw(ArgumentError("this form holds no rule to expand"))
    return exact_x(f.rule, f.k)
end
radical(f::QForm; kwargs...) = radical(_expanded(f); kwargs...)
radical_form(f::QForm; kwargs...) = radical_form(_expanded(f); kwargs...)
has_radical_form(f::QForm) = has_radical_form(_expanded(f))
xpolynomial(f::QForm) = xpolynomial(_expanded(f))
radicand(f::QForm) = radicand(_expanded(f))

end
