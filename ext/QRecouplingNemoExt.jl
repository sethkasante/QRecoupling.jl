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

end
