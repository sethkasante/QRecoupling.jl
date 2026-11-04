"""
Numerator, denominator and radical of the x-form, rendered — or described when too large to read.

Lives here because `show(::XValue)` is the only thing that needs it; losing it is what made `xform(v)`
throw `UndefVarError: _xvalue_str` rather than print.
"""
function _xvalue_str(xv::XValue)
    iszero(xv.num) && return "0"
    n = _xpoly_show(_toqq(xv.num); maxdeg = GENERIC_MAX_DEGREE, maxchars = GENERIC_MAX_CHARS)
    d = _xpoly_show(_toqq(xv.den); maxdeg = GENERIC_MAX_DEGREE, maxchars = GENERIC_MAX_CHARS)
    wrap(t) = occursin(" ", t) && !startswith(t, "(") ? "(" * t * ")" : t
    over = d == "1" ? "" : "/" * wrap(d)
    isempty(xv.rad) && return (over == "" ? n : wrap(n) * over)
    # the root last, as a form prints it
    root = "√(" * join("ψ" * to_subscript(e) for e in xv.rad) * ")"
    n == "1" && return root * over
    n == "−1" && return "−" * root * over
    return (over == "" ? n : wrap(n) * over) * " " * root
end

# printed as a form is, with what x and the ψ mean after it
Base.show(io::IO, v::XValue) = print(io, _xvalue_str(v), "   (x = q + q⁻¹",
                                     isempty(v.rad) ? "" : ", ψₑ the minimal polynomial of 2cos(2π/e)", ")")

