"""
Numerator, denominator and radical of the x-form, rendered — or described when too large to read.

Lives here because `show(::XValue)` is the only thing that needs it; losing it is what made `xform(v)`
throw `UndefVarError: _xvalue_str` rather than print.
"""
function _xvalue_str(xv::XValue)
    iszero(xv.num) && return "0"
    n = _xpoly_show(_toqq(xv.num); maxdeg = GENERIC_MAX_DEGREE, maxchars = GENERIC_MAX_CHARS)
    d = _xpoly_show(_toqq(xv.den); maxdeg = GENERIC_MAX_DEGREE, maxchars = GENERIC_MAX_CHARS)
    body = d == "1" ? n : (occursin(" ", n) && !startswith(n, "(") ? "(" * n * ")" : n) * " / " *
                          (occursin(" ", d) && !startswith(d, "(") ? "(" * d * ")" : d)
    isempty(xv.rad) && return body
    rad = join(["ψ" * to_subscript(e) for e in xv.rad], " · ")
    return "√(" * rad * ") · " * body
end

Base.show(io::IO,v::XValue) = print(io,"XValue(",_xvalue_str(v),")")
function Base.show(io::IO,::MIME"text/plain",v::XValue)
    println(io,"Exact expression in x = q + q⁻¹")
    print(io,"  = ",_xvalue_str(v))
    isempty(v.rad) || print(io,"\n  ψ_e(x) is the minimal polynomial of 2cos(2π/e)")
end

