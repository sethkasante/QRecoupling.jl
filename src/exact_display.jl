# Exact level values, continued from exact_radicals.jl: rendering of `RadExpr` and `ExactX`.

# ---------------------------------------------------------------------------------
#  Rendering
# ---------------------------------------------------------------------------------

_num_str(r::Rational{BigInt}) = denominator(r) == 1 ? string(abs(numerator(r))) :
                                string(abs(numerator(r))) * "/" * string(denominator(r))

"`c√m`, or `c√(…)` for a nested radicand, with a unit coefficient left implicit."
function _term_str(c::Rational{BigInt}, inner::RadExpr)
    body = is_rational(inner) ? "√" * _num_str(inner.rat) : "√(" * _rad_str(inner) * ")"
    a = abs(c)
    isone(a) && return body
    n = numerator(a); d = denominator(a)
    lead = isone(n) ? "" : string(n)
    return isone(d) ? lead * body : lead * body * "/" * string(d)
end

"""
A radical expression over a common denominator, so `3/2 − (1/2)√5` prints as `(3 − √5)/2` rather than as a
sum of fractions, and with the positive parts first, so it reads `(√5 − 3)/2` and not `(−3 + √5)/2`.
"""
function _rad_str(e::RadExpr)
    isempty(e.terms) && return (e.rat < 0 ? "−" : "") * _num_str(e.rat)
    den = denominator(e.rat)
    for (c, _) in e.terms
        den = lcm(den, denominator(c))
    end
    items = Tuple{Bool,String}[]                     # (negative?, body)
    r = e.rat * den
    iszero(r) || push!(items, (r < 0, _num_str(r)))
    for (c, inner) in e.terms
        push!(items, (c < 0, _term_str(c * den, inner)))
    end
    sort!(items, by = t -> t[1])                     # positives first; stable, so order is otherwise kept
    body = (items[1][1] ? "−" : "") * items[1][2]
    for (neg, str) in items[2:end]
        body *= (neg ? " − " : " + ") * str
    end
    den == 1 && return body
    return length(items) == 1 ? body * "/" * string(den) : "(" * body * ")/" * string(den)
end

Base.show(io::IO, e::RadExpr) = print(io, _rad_str(e))

"A polynomial in x over a common denominator: `(2x² − 7)/3`."
function _xpoly_str(f)
    degree(f) < 0 && return "0"
    den = BigInt(1)
    for i in 0:degree(f)
        den = lcm(den, denominator(Rational{BigInt}(coeff(f, i))))
    end
    out = IOBuffer()
    first = true
    for i in degree(f):-1:0
        c = Rational{BigInt}(coeff(f, i)) * den
        iszero(c) && continue
        n = numerator(c)
        if first
            n < 0 && print(out, "−")
            first = false
        else
            print(out, n < 0 ? " − " : " + ")
        end
        a = abs(n)
        (isone(a) && i > 0) || print(out, a)
        i == 0 && continue
        print(out, "x")
        i == 1 || print(out, to_superscript(i))
    end
    body = String(take!(out))
    isempty(body) && (body = "0")
    den == 1 && return body
    return occursin(" ", body) ? "(" * body * ")/" * string(den) : body * "/" * string(den)
end

"""
Longest polynomial the display writes out in full. Past it the line is a size, not a formula: at k = 420
the value is a degree-208 polynomial with 240-bit coefficients, and printing it is three screens of
digits that tell a reader nothing `xpolynomial(v)` would not tell them better.
"""
const XPOLY_MAX_DEGREE = 24
const XPOLY_MAX_CHARS = 240

"A polynomial too large to read, described instead: the same courtesy `phi_form` extends."
function _xpoly_summary(f)
    d = degree(f)
    cs = [Rational{BigInt}(coeff(f, i)) for i in 0:d]
    den = BigInt(1)
    for c in cs
        den = lcm(den, denominator(c))
    end
    bits = maximum(c -> ndigits(numerator(c * den); base = 2), cs; init = 0)
    nz = count(!iszero, cs)
    body = "an integer polynomial in x of degree $d, $nz nonzero coefficients, ≤$bits bits"
    return den == 1 ? "(" * body * ")" : "(" * body * ")/" * string(den)
end

"`_xpoly_str` when it fits the budget, `nothing` when it does not."
function _xpoly_fits(f; maxdeg::Int = XPOLY_MAX_DEGREE, maxchars::Int = XPOLY_MAX_CHARS)
    degree(f) > maxdeg && return nothing
    str = _xpoly_str(f)
    return length(str) > maxchars ? nothing : str
end

"`_xpoly_str`, or a description of it when writing it out would be worse than useless."
function _xpoly_show(f; maxdeg::Int = XPOLY_MAX_DEGREE, maxchars::Int = XPOLY_MAX_CHARS)
    str = _xpoly_fits(f; maxdeg = maxdeg, maxchars = maxchars)
    return str === nothing ? _xpoly_summary(f) : str
end

"`degree 208, 105 terms, ≤144 bits`, the shape of a polynomial without the polynomial."
function _xpoly_dims(f)
    d = degree(f)
    d < 0 && return "0"
    cs = [Rational{BigInt}(coeff(f, i)) for i in 0:d]
    den = BigInt(1)
    for c in cs
        den = lcm(den, denominator(c))
    end
    bits = maximum(c -> ndigits(numerator(c * den); base = 2), cs; init = 0)
    out = "degree $d, $(count(!iszero, cs)) terms, ≤$bits bits"
    return den == 1 ? out : "(" * out * ")/" * string(den)
end

"Longest square class the sizes line writes out; past it the count says as much and fits."
const SQCLASS_MAX_SHOWN = 6

"The class under the root, as ψ indices or — past a handful of them — as how many there are."
function _sqclass_str(cls::Vector{Int})
    isempty(cls) && return "1"
    length(cls) > SQCLASS_MAX_SHOWN && return string(length(cls)) * " ψ factors"
    return join(["ψ" * to_subscript(e) for e in cls], "·")
end

"Is either half of `P(x)·√(R(x))` past what the display will write out?"
function _xform_long(v::ExactX)
    iszero(v.p) && return false
    _xpoly_fits(v.p) === nothing && return true
    return !isempty(v.sqclass) && _xpoly_fits(v.r) === nothing
end

"The sizes of `P` and `R`, one line each, for the lines under a schematic form."
function _xform_dims(v::ExactX)
    isempty(v.sqclass) && return ["P: " * _xpoly_dims(v.p) * ",  R = 1"]
    return ["P: " * _xpoly_dims(v.p),
            "R: " * _xpoly_dims(v.r) * ",  " * _sqclass_str(v.sqclass)]
end

"""
The `P(x)·√(R(x))` line: what printing a value shows, always available and never computing anything —
no minimal polynomial, no Lagrange descent, nothing that can fail or take unbounded time. A radicand that
reduced to a bare constant keeps its class (it is not a rational square, or the constructor would have
folded it) and is written `√2` rather than `√(2)`.
"""
function _xform_str(v::ExactX)
    iszero(v.p) && return "0"
    # Past the budget the honest line is the *shape* of the value, `P(x)·√(R(x))`, with the sizes on the
    # line below: a reader takes nothing from three screens of digits, and `v.x_value` gives both halves
    # in full to anyone who wants them.
    _xform_long(v) && return isempty(v.sqclass) ? "P(x)" : "P(x) · √(R(x))"
    ps = _xpoly_show(v.p)
    isempty(v.sqclass) && return ps
    rs = _xpoly_show(v.r)
    root = occursin(" ", rs) || occursin("x", rs) ? "√(" * rs * ")" : "√" * rs
    ps == "1" && return root
    ps == "-1" && return "−" * root
    need = occursin(" ", ps) && !startswith(ps, "(")
    return (need ? "(" * ps * ")" : ps) * " · " * root
end

function Base.show(io::IO, v::ExactX)
    print(io, "ExactX(k = ", v.k, ", ", _xform_str(v), ")")
end

"The value's degree when that is cheap to know, `0` when it is only worth a sentence."
_cheap_degree(u, h::Int, d::Int) = d <= RADICAL_MINPOLY_MAX_DEGREE ? _value_degree(u, h) : 0

"""
    _radical_view(v; maxlen, degree_limit, allowed) -> (kind, expr, degree)

What [`radical`](@ref) can say about this value, decided in one place so that the answer and the sentence
explaining it cannot disagree. `:zero` and `:ok` carry an expression; `:none`, `:untried`, `:long` and
`:failed` become a [`NoRadical`](@ref) and are documented there.
"""
function _radical_view(v::ExactX; maxlen::Int = 80,
                      degree_limit::Int = RADICAL_DESCENT_MAX_DEGREE, allowed::Bool = true)
    iszero(v.p) && return (:zero, nothing, 0)
    h = v.k + 2
    d = euler_phi(2h) ÷ 2
    # A rational value, or a rational multiple of one surd, is its own radical: no field work, no
    # descent, no cap. This must come first, or a degenerate value at a big level is refused for the
    # size of a field it does not live in.
    if degree(v.p) <= 0 && (isempty(v.sqclass) || degree(v.r) <= 0)
        e = radical_form(v; maxlen = maxlen, degree_limit = 0)
        e === nothing && return (maxlen > 0 ? :long : :failed, nothing, 1)
        return (:ok, e, 1)
    end
    allowed || return (:untried, nothing, 0)
    vd = d
    # Three questions, in increasing cost, and none of them asked unless the answer is needed:
    #   1. does a radical exist?     `G_odd ⊆ Stab(v²)` — a couple of conjugations
    #   2. is the degree affordable? the orbit, stopped as soon as it passes the limit
    #   3. what is the degree?       the whole orbit, and only when it is small enough to be cheap
    # The old route asked (3) always, through a minimal polynomial of the field's degree: 164 ms at
    # k = 420 where (1) settles it in 12 ms.
    if !has_radical_form(v.k) || d > degree_limit
        u = _square(v)
        # Existence first, and cheaply: it is the question, and the degree is only the sentence.
        _degree_is_2power(u, h) || return (:none, nothing, _cheap_degree(u, h, d))
        d <= RADICAL_MINPOLY_MAX_DEGREE || return (:untried, nothing, 0)
        vd = _value_degree(u, h)
        vd <= degree_limit || return (:untried, nothing, vd)
    end
    vd <= degree_limit || return (:untried, nothing, vd)
    e = radical_form(v; maxlen = maxlen, degree_limit = 0)   # the degree is already decided
    # with no length budget the only way back is a descent that could not certify a sign, which is not
    # the same statement as "longer than you asked for"
    e === nothing && return (maxlen > 0 ? :long : :failed, nothing, vd)
    return (:ok, e, vd)
end

"""
Printing an `ExactX` does no arithmetic of any kind.

The stored `P(x)·√(R(x))` is exact at every level and free to produce; the radical form is a question
the reader asks (`radical(v)`, which prints the number with it), and the certified `≈` is a Horner pass
this display used to pay for at every value — `IOContext(io, :approximate => true)` asks for it back,
and `Float64(v)` was always the direct way. Past the length budget the line becomes the *shape* of the
value with its sizes underneath, which is what a degree-208 polynomial has to say for itself. Nothing
else is printed: the display is the value, and the rest is what the reader calls for.
"""
function Base.show(io::IO, ::MIME"text/plain", v::ExactX)
    h = v.k + 2
    d = euler_phi(2h) ÷ 2
    lines = ["Exact value at level k = " * string(v.k) * "   (x = 2cos(π/" * string(h) *
             "), degree " * string(d) * " over ℚ)"]
    if iszero(v.p)
        print(io, lines[1], "\n  = 0")
        return
    end
    push!(lines, "  = " * _xform_str(v))
    if _xform_long(v)
        for l in _xform_dims(v)
            push!(lines, "    " * l)
        end
        push!(lines, "    `v.x_value` gives (P, R) in full")
    end
    if get(io, :approximate, false)
        nv = numeric_value(v)
        push!(lines, nv === nothing ?
              "  (the value cancels beyond " * string(EVAL_MAX_BITS) * " bits at x = 2cos(π/" *
              string(h) * "); for a number use q6j(…; k = " * string(v.k) * "))" :
              "  ≈ " * string(Float64(nv)))
    end
    print(io, join(lines, "\n"))
end
