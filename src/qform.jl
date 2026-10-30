# ---------------------------------------------------------------------------------
#  A value written in quantum integers
#
#  A rational is a ratio of integers; a recoupling coefficient is built the same way from q-integers:
#
#      sign · c · ∏[n]^e · √(r ∏[n]^f) · Σ_terms c_t ∏[n]^(e_t)
#
#  Generic q. The integers [n], n ≥ 2, are multiplicatively independent, since [n] = ∏_{d|n, d>1} Ψ_d and the
#  Ψ_d are coprime. A product of them therefore has one set of exponents, read off the rule by counting.
#
#  A level, h = k + 2. Three things change:
#    the labels fold,        [n + h] = −[n], [h − n] = [n], so 1 ≤ n ≤ h/2 and every [n] is positive;
#    the vanishing factors   leave their limit [jh]/[h] = (−1)^(j−1) j, and balance across a contributing term;
#    the rest obey           [pb] ∏_{j=1}^{p−1} [jh/p] = p ∏_{j=0}^{p−1} [b + jh/p]    for each prime p dividing h.
#  The last is the distribution relation of sin(pθ) = 2^(p−1) ∏ sin(θ + jπ/p), with the power of two removed by
#  its limit b → 0. Reducing by the lattice of these relations writes large labels through small ones.
#  The relations are valid identities, so reducing never changes a value. Two equal values can still print
#  differently, so equality belongs to `Exact(k)`, not to this form.
# ---------------------------------------------------------------------------------

"""
    QForm

A value written in q-integers, as returned by [`qform`](@ref). `f(q)` evaluates a generic form at a real
positive `q`, `q = 1` included; `numeric_value(f)` and `float(f)` evaluate a form at its level.
"""
struct QForm{R<:Rational}
    k::Union{Nothing,Int}
    sign::Int
    coef::R
    out::Vector{Pair{Int,Int}}                                # ∏[n]^e
    radcoef::R                                                # r under the square root
    rad::Vector{Pair{Int,Int}}                                # ∏[n]^f under the square root
    terms::Vector{Tuple{R,Vector{Pair{Int,Int}}}}
end

_qform_zero(k) = QForm{Rational{Int}}(k, 0, 0 // 1, Pair{Int,Int}[], 1 // 1, Pair{Int,Int}[],
                                      Tuple{Rational{Int},Vector{Pair{Int,Int}}}[])

Base.iszero(f::QForm) = iszero(f.sign)

# ---- the relations at a level ----

"Is `p` prime, for the small numbers that divide a level's h."
_small_prime(p::Int) = p >= 2 && all(d -> p % d != 0, 2:isqrt(p))

"""
Echelon form of the distribution relations at h, from the largest label down: rows `(pivot label, exponents,
constant)`, meaning ∏[n]^(exponents[n]) = constant, with the exponent of the pivot label positive and every
larger label absent. Labels are 2 ≤ n ≤ h/2.
"""
function _level_relations(h::Int)
    top = h ÷ 2
    rows = Tuple{Vector{Int},Rational{BigInt}}[]
    for p in 2:h
        (h % p == 0 && _small_prime(p)) || continue
        step = h ÷ p
        for b in 1:top
            (p * b) % h == 0 && continue
            any(j -> (b + j * step) % h == 0, 0:p-1) && continue
            e = zeros(Int, top)
            for (i, c) in Iterators.flatten(((p * b => 1,), ((j * step => 1) for j in 1:p-1),
                                             ((b + j * step => -1) for j in 0:p-1)))
                m, _ = _fold_label(i, h)
                m >= 2 && (e[m] += c)
            end
            any(!iszero, e) && push!(rows, (e, p // big(1)))
        end
    end
    out = Tuple{Int,Vector{Int},Rational{BigInt}}[]
    for col in top:-1:2
        while true
            nz = [r for r in rows if !iszero(r[1][col])]
            length(nz) <= 1 && break
            piv = argmin(r -> abs(r[1][col]), nz)
            for i in eachindex(rows)
                r = rows[i]
                (r === piv || iszero(r[1][col])) && continue
                q = div(r[1][col], piv[1][col])
                rows[i] = (r[1] .- q .* piv[1], r[2] / piv[2]^q)
            end
        end
        k = findfirst(r -> !iszero(r[1][col]), rows)
        k === nothing && continue
        r = popat!(rows, k)
        r[1][col] < 0 && (r = (-r[1], inv(r[2])))
        push!(out, (col, r[1], r[2]))
    end
    return out
end

const LEVEL_RELATIONS = LRU{Int,Vector{Tuple{Int,Vector{Int},Rational{BigInt}}}}(maxsize = 64)
_relations_at(h::Int) = get!(() -> _level_relations(h), LEVEL_RELATIONS, h)

"Reduce the exponents `e` (labels ≤ h/2) by the relations, returning how many times each row was taken out."
function _reduce_exponents!(e::Vector{Int}, rows)
    qs = zeros(Int, length(rows))
    for (r, (piv, row, _)) in enumerate(rows)
        q = fld(e[piv], row[piv])
        iszero(q) && continue
        @inbounds for i in eachindex(e)
            e[i] -= q * row[i]
        end
        qs[r] = q
    end
    return qs
end

"The constant a product picked up from `_reduce_exponents!`."
function _relation_constant(rows, qs::Vector{Int}, ::Type{R}) where {R}
    c = one(R)
    for (r, q) in enumerate(qs)
        iszero(q) || (c *= R(rows[r][3])^q)
    end
    return c
end

"Terms and total exponent of a form, read from the exponents alone: enough to tell which of two reads shorter."
function _exponent_cost(pre::Vector{Int}, es::Vector{Vector{Int}})
    ts = unique(es)
    g = reduce((a, b) -> min.(a, b), ts)
    total = sum(sum(abs(a - b) for (a, b) in zip(t, g); init = 0) for t in ts; init = 0)
    return (length(ts), total + sum(abs(a + 2b) for (a, b) in zip(pre, g); init = 0) ÷ 2)
end

# ---- exponents from a rule ----

"Label and sign of [i] at a level, for i not a multiple of h: 1 ≤ label ≤ h/2 and [i] = sign·[label]."
function _fold_label(i::Int, h::Int)
    i = mod(i, 2h); s = 1
    if i > h
        s = -1; i -= h
    end
    i > h ÷ 2 && (i = h - i)
    return i, s
end

"The exponents of [2], …, [N] in ∏_f [arg_f(z)]!^(c_f), by a difference array."
function _term_exponents!(e::Vector{Int}, s::FactorialSum, z::Int)
    fill!(e, 0)
    for f in s.fac
        n = _arg(f, z)
        n >= 2 || continue
        e[2] += f.c
        n < length(e) && (e[n+1] -= f.c)
    end
    for i in 3:length(e)
        e[i] += e[i-1]
    end
    return e
end

"""
Folded labels of 1..N at a level: `lab[i]` is the label, 0 for a multiple of h, and `neg[i]` whether the sign is
negative. Filled once per call, because the exponents of every term run over them.
"""
function _fold_table(N::Int, h::Int)
    lab = Vector{Int}(undef, N); neg = Vector{Bool}(undef, N)
    @inbounds for i in 1:N
        if i % h == 0
            lab[i] = 0; neg[i] = iseven(i ÷ h)
        else
            lab[i], t = _fold_label(i, h)
            neg[i] = t < 0
        end
    end
    return lab, neg
end

"Fold raw exponents to a level: `(exponents on labels ≤ h/2, rational factor, sign, power of [h])`."
function _fold_exponents(raw::Vector{Int}, h::Int, lab, neg, ::Type{R}) where {R}
    e = zeros(Int, h ÷ 2); r = one(R); s = 1; nh = 0
    @inbounds for i in 2:length(raw)
        c = raw[i]
        iszero(c) && continue
        m = lab[i]
        if m == 0
            r *= R(i ÷ h)^c
            nh += c
        elseif m >= 2
            e[m] += c
        end
        neg[i] && isodd(c) && (s = -s)
    end
    return e, r, s, nh
end

_nonzero_pairs(e::AbstractVector{Int}) = [i => e[i] for i in eachindex(e) if !iszero(e[i])]

"Move square factors out of a positive rational: returns `(outside, what stays under the root)`."
function _split_square(r::R) where {R<:Rational}
    isone(r) && return (r, r)
    function part(n)
        o = one(n); rest = n; p = 2
        while p * p <= rest && p < 100_000
            while rest % (p * p) == 0
                rest ÷= p * p; o *= p
            end
            p += 1
        end
        s = isqrt(rest)
        s * s == rest && (o *= s; rest = one(n))
        return o, rest
    end
    a, b = part(numerator(r)); c, d = part(denominator(r))
    return R(a // c), R(b // d)
end

"""
    qform(rule; k = nothing, reduce = true) -> QForm
    qform(symbol, labels...; k = nothing, reduce = true) -> QForm
    qform(q6j(Symbolic(), labels...); k = nothing, reduce = true) -> QForm

A value written in quantum integers `[n]`, the form recoupling theory is written in. Without `k` it holds at
every generic q, and `q = 1` turns it into integers. With `k` the labels fold to `n ≤ (k+2)/2` and vanishing
factors leave their limit. With `reduce = true` the distribution relations of the level are applied too: terms
can then merge or cancel, and a zero they reveal is returned as zero. The shorter of the two forms is returned.
Rules with `q`-power weights (`qcg`, `q3j`) are not covered.

Two equal values can print differently; `Exact(k)` decides equality. A sum that does not visibly cancel here can
still vanish: `iszero_at` decides that.

```julia
qform(q6j, 1, 1, 1, 1, 1, 1)          # ([5] − 1)/[4]!, which is 1/6 at q = 1
qform(q6j, 1, 1, 1, 1, 1, 1; k = 3)   # −1/[2]²
qform(q6j, 1, 1, 1, 2, 1, 1; k = 5)   # [2]/[3]²
```
"""
function qform(s::FactorialSum; k = nothing, reduce::Bool = true)
    try
        return _qform(s, k, reduce, Rational{Int})
    catch e
        (e isa OverflowError || e isa InexactError) || rethrow()
        return _qform(s, k, reduce, Rational{BigInt})
    end
end

"How long a form reads: its terms first, then its exponents and the digits of its numbers."
function _qform_cost(f::QForm)
    digits(c) = ndigits(numerator(c)) + ndigits(denominator(c)) - 2
    exps(p) = sum(abs(e) for (_, e) in p; init = 0)
    return (length(f.terms), exps(f.out) + exps(f.rad) + sum(exps(p) + digits(c) for (c, p) in f.terms; init = 0) +
                             digits(f.coef) + digits(f.radcoef))
end

function _qform(s::FactorialSum, k, reduce::Bool, ::Type{R}) where {R}
    level = k === nothing ? nothing : Int(k)
    is_empty_sum(s) && return _qform_zero(level)
    if level !== nothing
        st, segs = classify_at_level(s, level)
        st === :pole && throw(ArgumentError("the rule is singular at level $level"))
        st === :finite || return _qform_zero(level)
        zs = Int[z for seg in segs for z in seg]
    else
        zs = collect(s.zlo:s.zhi)
    end
    N = 2
    for (n, _) in s.pre
        N = max(N, Int(n))
    end
    for f in s.fac, z in (first(zs), last(zs))
        N = max(N, _arg(f, z))
    end
    raw = zeros(Int, N)
    # the prefactor squared, as doubled exponents
    pre2 = zeros(Int, N)
    for (n, c) in s.pre
        n >= 2 || continue
        pre2[2] += (s.sqrt_pre ? 1 : 2) * c
        n < N && (pre2[Int(n)+1] -= (s.sqrt_pre ? 1 : 2) * c)
    end
    for i in 3:N
        pre2[i] += pre2[i-1]
    end
    es = Vector{Vector{Int}}(undef, length(zs)); cs = Vector{R}(undef, length(zs))
    if level === nothing
        pre = pre2; prerat = one(R)
        for (t, z) in enumerate(zs)
            es[t] = copy(_term_exponents!(raw, s, z))
            cs[t] = s.alternating && isodd(z) ? -one(R) : one(R)
        end
        return _qform_assemble(level, Int(s.sign0), pre, prerat, es, cs)
    end
    h = level + 2
    lab, neg = _fold_table(N, h)
    pre, prerat, ps, nhp = _fold_exponents(pre2, h, lab, neg, R)
    ps < 0 && throw(ArgumentError("the radicand is negative at level $level"))
    for (t, z) in enumerate(zs)
        e, r, ts, nh = _fold_exponents(_term_exponents!(raw, s, z), h, lab, neg, R)
        nhp + 2nh == 0 || throw(ArgumentError("the terms of this rule have different valuations at level $level"))
        c = ts * r
        s.alternating && isodd(z) && (c = -c)
        es[t] = e; cs[t] = c
    end
    plain = _qform_assemble(level, Int(s.sign0), pre, prerat, es, cs)
    reduce || return plain
    rows = _relations_at(h)
    (isempty(rows) || iszero(plain)) && return plain         # h prime: the folded labels are independent
    # Exponents first. Constants are only worth computing for a form that merges terms or reads shorter.
    es2 = copy.(es); pre2r = copy(pre)
    qs = [_reduce_exponents!(e, rows) for e in es2]
    qp = _reduce_exponents!(pre2r, rows)
    _exponent_cost(pre2r, es2) < _exponent_cost(pre, es) || return plain
    cs2 = [cs[t] * _relation_constant(rows, qs[t], R) for t in eachindex(cs)]
    reduced = _qform_assemble(level, Int(s.sign0), pre2r, prerat * _relation_constant(rows, qp, R), es2, cs2)
    iszero(reduced) && return reduced
    return _qform_cost(reduced) < _qform_cost(plain) ? reduced : plain
end

"Merge equal terms, move what every term shares in front, and take even exponents out of the root."
function _qform_assemble(level, sign::Int, pre::Vector{Int}, prerat::R, es::Vector{Vector{Int}}, cs::Vector{R}) where {R}
    seen = Dict{Vector{Int},Int}(); keep = Int[]; coef = R[]
    for t in eachindex(es)
        i = get!(seen, es[t], length(keep) + 1)
        if i > length(keep)
            push!(keep, t); push!(coef, cs[t])
        else
            coef[i] += cs[t]
        end
    end
    live = [i for i in eachindex(keep) if !iszero(coef[i])]
    isempty(live) && return _qform_zero(level)
    score = [sum(abs, es[keep[i]]) for i in live]
    order = sortperm(score; rev = true)                      # stable: terms of equal length keep their order
    ts = live[order]
    n = length(pre)
    g = copy(es[keep[first(ts)]])
    for i in ts
        e = es[keep[i]]
        @inbounds for j in 1:n
            e[j] < g[j] && (g[j] = e[j])
        end
    end
    top = zero(numerator(prerat)); bot = one(top)
    for i in ts
        top = gcd(top, numerator(coef[i])); bot = lcm(bot, denominator(coef[i]))
    end
    scale = R(top // bot)
    coef[first(ts)] < 0 && (scale = -scale)
    terms = Vector{Tuple{R,Vector{Pair{Int,Int}}}}(undef, length(ts))
    for (a, i) in enumerate(ts)
        e = es[keep[i]]
        p = Pair{Int,Int}[]
        @inbounds for j in 1:n
            d = e[j] - g[j]
            iszero(d) || push!(p, j => d)
        end
        terms[a] = (coef[i] / scale, p)
    end
    out = Pair{Int,Int}[]; rad = Pair{Int,Int}[]
    @inbounds for j in 1:n
        a = pre[j] + 2 * g[j]
        iszero(a) && continue
        o = div(a, 2); f = a - 2o
        iszero(o) || push!(out, j => o)
        iszero(f) || push!(rad, j => f)
    end
    outside, under = _split_square(prerat)
    c = abs(scale) * outside
    return QForm{R}(level, scale < 0 ? -sign : sign, c, out, under, rad, terms)
end

qform(v::SymbolicValue; kw...) = qform(v.rule; kw...)

function qform(f::Function, args...; k = nothing, reduce::Bool = true)
    k !== nothing && !_admissible_at(f, Int(k), args...) && return _qform_zero(Int(k))
    return qform(_rule_for(f, args...); k = k, reduce = reduce)
end

# ---- numbers ----

function _qform_value(f::QForm, qn)
    iszero(f) && return zero(qn(2))
    mono(p) = (v = one(qn(2)); for (n, e) in p; v *= qn(n)^e; end; v)
    total = zero(qn(2))
    for (c, p) in f.terms
        total += oftype(total, numerator(c)) / oftype(total, denominator(c)) * mono(p)
    end
    w = oftype(total, numerator(f.coef)) / oftype(total, denominator(f.coef))
    u = oftype(total, numerator(f.radcoef)) / oftype(total, denominator(f.radcoef))
    return f.sign * w * mono(f.out) * sqrt(u * mono(f.rad)) * total
end

function (f::QForm)(q::Real)
    f.k === nothing || throw(ArgumentError("this form is at level $(f.k); use `numeric_value` or `float`"))
    q > 0 || throw(DomainError(q, "a QForm is evaluated at positive real q"))
    Q = float(q)
    return _qform_value(f, isone(Q) ? (n -> oftype(Q, n)) : (n -> (Q^n - Q^(-n)) / (Q - inv(Q))))
end

function _qform_level(f::QForm)
    f.k === nothing && throw(ArgumentError("a generic form is evaluated with `f(q)`"))
    h = f.k + 2
    return _qform_value(f, n -> sinpi(BigFloat(n) / h) / sinpi(inv(BigFloat(h))))
end

function numeric_value(f::QForm; bits::Integer = precision(BigFloat))
    v = setprecision(() -> _qform_level(f), BigFloat, Int(bits) + 32)
    return setprecision(() -> BigFloat(v), BigFloat, Int(bits))
end
Base.float(f::QForm) = Float64(numeric_value(f; bits = 128))
Base.Float64(f::QForm) = float(f)

# ---- display ----

_q_power(n::Int, e::Int) = "[" * string(n) * "]" * (abs(e) == 1 ? "" : to_superscript(abs(e)))
_q_factorial(n::Int, e::Int) = "[" * string(n) * "]!" * (abs(e) == 1 ? "" : to_superscript(abs(e)))

"""
A product of q-integers as `(above, below)` text. A run of consecutive labels with one exponent reads as a
factorial when that is shorter, [2]…[b] = [b]! and [a]…[b] = [b]!/[a−1]!, and factorials from different runs
combine, so the two ends of neighbouring runs telescope. Single labels come first.
"""
function _q_product(p::Vector{Pair{Int,Int}})
    single = Dict{Int,Int}(); fact = Dict{Int,Int}()
    add!(d, n, e) = (d[n] = get(d, n, 0) + e)
    i = 1
    while i <= length(p)
        a, e = p[i]
        j = i
        while j < length(p) && p[j+1].first == p[j].first + 1 && p[j+1].second == e
            j += 1
        end
        b = p[j].first; len = j - i + 1
        if a == 2 && len >= 2
            add!(fact, b, e)
        elseif a > 2 && len >= 4
            add!(fact, b, e)
            a == 3 ? add!(single, 2, -e) : add!(fact, a - 1, -e)
        else
            for t in i:j
                add!(single, p[t].first, e)
            end
        end
        i = j + 1
    end
    up = String[]; dn = String[]
    for (d, text) in ((single, _q_power), (fact, _q_factorial)), n in sort!(collect(keys(d)))
        e = d[n]
        iszero(e) || push!(e > 0 ? up : dn, text(n, e))
    end
    return join(up), join(dn)
end

function _q_fraction(c::Rational, p::Vector{Pair{Int,Int}}, inner::String = "")
    up, down = _q_product(p)
    a = numerator(c); b = denominator(c)
    top = (isone(a) && !(isempty(up) && isempty(inner)) ? "" : string(a)) * up * inner
    bot = (isone(b) ? "" : string(b)) * down
    isempty(bot) && return top
    factors = count(==('['), bot) + (isone(b) ? 0 : 1)
    return top * "/" * (factors > 1 ? "(" * bot * ")" : bot)
end

const QFORM_SHOWN_TERMS = 6

"Drop a pair of parentheses that encloses the whole text."
function _unwrap(text::String)
    startswith(text, "(") && endswith(text, ")") || return text
    depth = 0
    for (i, c) in pairs(text)
        depth += c == '(' ? 1 : c == ')' ? -1 : 0
        depth == 0 && i < lastindex(text) && return text
    end
    return text[nextind(text, 1):prevind(text, lastindex(text))]
end

function _q_sum(terms; most::Int = QFORM_SHOWN_TERMS)
    length(terms) == 1 && isone(first(terms)[1]) && isempty(first(terms)[2]) && return ""
    shown = first(terms, most)
    parts = String[]
    for (i, (c, p)) in enumerate(shown)
        t = _q_fraction(abs(c), p)
        push!(parts, i == 1 ? (c < 0 ? "−" : "") * t : (c < 0 ? " − " : " + ") * t)
    end
    extra = length(terms) - length(shown)
    extra > 0 && push!(parts, " + … ($(length(terms)) terms)")
    return length(terms) == 1 ? join(parts) : "(" * join(parts) * ")"
end

"The form as text, without its level, showing at most `most` terms of its sum."
function _qform_body(f::QForm; most::Int = QFORM_SHOWN_TERMS)
    iszero(f) && return "0"
    text = (f.sign < 0 ? "−" : "") * _q_fraction(f.coef, f.out, _q_sum(f.terms; most = most))
    isempty(f.rad) && isone(f.radcoef) && return _unwrap(text)
    return text * " √(" * _q_fraction(f.radcoef, f.rad) * ")"
end

function Base.show(io::IO, f::QForm)
    print(io, _qform_body(f))
    f.k === nothing || print(io, "   (k = ", f.k, ")")
end
