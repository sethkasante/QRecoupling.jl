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
Longest polynomial the display writes out in full; past it the line shows the shape and sizes instead.
"""
const XPOLY_MAX_DEGREE = 24
const XPOLY_MAX_CHARS = 240

"A polynomial too large to read, described instead."
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

"Largest field degree at which printing looks for a radical, and the longest radical it will print."
const QFORM_RADICAL_DEGREE = 16
const QFORM_RADICAL_CHARS = 80

"Roughly how long a form reads before its sum is cut short."
const QFORM_DISPLAY_MAX = 80

"""
A polynomial in x = [2] as `(coefficient, labels)` terms: the leading x^m goes with [m+1], which is monic of
degree m, and [n+1] = x[n] − [n−1] gives the rest. Labels fold to n ≤ h/2 and [h] = 0 drops out.
"""
function _xpoly_terms(f, h::Int)
    d = degree(f)
    cs = [Rational{BigInt}(coeff(f, i)) for i in 0:d]
    us = Vector{Vector{BigInt}}(undef, d + 1)                # us[m+1] holds the coefficients of [m+1]
    for m in 0:d
        us[m+1] = m == 0 ? BigInt[1] : m == 1 ? BigInt[0, 1] :
                  [(i >= 2 ? us[m][i-1] : big(0)) - (i <= m - 1 ? us[m-1][i] : big(0)) for i in 1:m+1]
    end
    labels = Dict{Int,Rational{BigInt}}()
    for m in d:-1:0
        a = cs[m+1]
        iszero(a) && continue
        cs[1:m+1] .-= a .* us[m+1]
        (m + 1) % h == 0 && continue
        n, t = _fold_label(m + 1, h)
        labels[n] = get(labels, n, 0 // big(1)) + t * a
    end
    return [(labels[n], n == 1 ? Pair{Int,Int}[] : [n => 1]) for n in sort!(collect(keys(labels)); rev = true)
            if !iszero(labels[n])]
end



"`P(x)·√(R(x))` in q-integers, for a value with no rule behind it, such as the result of arithmetic."
function _xvalue_str(v::ExactX; most::Int)
    h = v.k + 2
    text = _q_sum_str(_xpoly_terms(v.p, h); most = most)
    isempty(v.sqclass) && return _unwrap(text)
    return text * " √(" * _unwrap(_q_sum_str(_xpoly_terms(v.r, h); most = most)) * ")"
end

"A polynomial over ℚ as a string with integer coefficients, and the common denominator taken out."
function _integer_poly_str(f)
    d = denominator(f)
    return _xpoly_show(f * d; maxdeg = GENERIC_MAX_DEGREE, maxchars = GENERIC_MAX_CHARS), d
end

"`P(x)√(R(x))`, written with integer coefficients over one denominator."
function _level_x_str(v::ExactX)
    iszero(v.p) && return "0"
    n, d = _integer_poly_str(v.p)
    body = isone(d) ? n : (occursin(" ", n) ? "(" * n * ")" : n) * "/" * string(d)
    isempty(v.sqclass) && return body
    r, rd = _integer_poly_str(v.r)
    root = "√(" * (isone(rd) ? r : (occursin(" ", r) ? "(" * r * ")" : r) * "/" * string(rd)) * ")"
    body == "1" && return root
    body == "−1" && return "−" * root
    return body * " " * root
end

# printed as a form is, with what x means after the level
Base.show(io::IO, v::LevelXValue) =
    print(io, _level_x_str(v.value), "   (k = ", v.value.k, ", x = 2cos(π/", v.value.k + 2, "))")

"""
The value in q-integers: from its rule when it has one, which merges the labels at the level, and from
`P(x)·√(R(x))` otherwise. A long sum is cut after as many terms as fit.
"""
function _qform_str(v::ExactX)
    iszero(v.p) && return "0"
    f = nothing
    if v.rule !== nothing
        try
            f = qform(v.rule; k = v.k)
        catch e
            e isa ArgumentError || rethrow()                 # terms of different valuations at this level
        end
    end
    render = f === nothing || iszero(f) ? (m -> _xvalue_str(v; most = m)) : (m -> _qform_body(f; most = m))
    most = QFORM_SHOWN_TERMS
    text = render(most)
    while most > 1 && length(text) > QFORM_DISPLAY_MAX
        most -= 1
        text = render(most)
    end
    return text
end

"The radical to print beside a value, when one exists at a field degree worth trying and it is short."
function _radical_str(v::ExactX)
    euler_phi(2 * (v.k + 2)) ÷ 2 <= QFORM_RADICAL_DEGREE || degree(v.p) <= 0 || return nothing
    kind, e, _ = _radical_view(v; maxlen = QFORM_RADICAL_CHARS, degree_limit = QFORM_RADICAL_DEGREE)
    kind === :ok || return nothing
    text = string(e)
    return length(text) <= QFORM_RADICAL_CHARS ? text : nothing
end

function Base.show(io::IO, v::ExactX)
    print(io, "ExactX(k = ", v.k, ", ", _qform_str(v), ")")
end

"The value's degree when that is cheap to know, `0` when it is only worth a sentence."
_cheap_degree(u, h::Int, d::Int) = d <= RADICAL_MINPOLY_MAX_DEGREE ? _value_degree(u, h) : 0

"""
    _radical_view(v; maxlen, degree_limit, allowed) -> (kind, expr, degree)

What [`radical`](@ref) can say about this value, decided in one place so answer and explanation agree.
`:zero` and `:ok` carry an expression; the other kinds become a [`NoRadical`](@ref).
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
    # Cheapest first: does a radical exist (G_odd ⊆ Stab(v²)), is the degree affordable (orbit, stopped at
    # the limit), and only then what the degree is.
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
Printing an `ExactX` shows its value in q-integers, cut short when it is long, and the radical when one exists
and is short. `radical(v)`, `xform(v)` and `v.x_value` give the other forms in full.
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
    nv = numeric_value(v)
    push!(lines, nv === nothing ?
          "  (the value cancels beyond " * string(EVAL_MAX_BITS) * " bits at x = 2cos(π/" *
          string(h) * "); for a number use q6j(…; k = " * string(v.k) * "))" :
          "  ≈ " * string(Float64(nv)))
    push!(lines, "  = " * _qform_str(v))
    r = _radical_str(v)
    r === nothing || push!(lines, "  = " * r)
    print(io, join(lines, "\n"))
end
