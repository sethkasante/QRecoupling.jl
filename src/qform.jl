# ---------------------------------------------------------------------------------
#  A value written in quantum integers
#
#  A rational is a ratio of integers; a recoupling coefficient is built the same way from q-integers:
#
#      sign · c · q^(a/2) · ∏[n]^e · √(r ∏[n]^f) · Σ_terms c_t q^(b_t) ∏[n]^(e_t)
#
#  The powers of q come from the weights of `qcg` and `q3j`, and are absent from every other symbol.
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
#  Powers of q reduce by q^h = −1. The relations hold among the [n], so they act on each term whatever its power.
#  The relations are valid identities, so reducing never changes a value. Two equal values can still print
#  differently, so equality is decided by `==`, not by reading the form.
# ---------------------------------------------------------------------------------

"""
    QForm

A value written in q-integers, as returned by `Exact()`, `Exact(k)` and [`qform`](@ref).

A generic form holds at every q. `f(q)` evaluates it at a real positive `q`, and used as a number it means its value at
q = 1: [`classical`](@ref) gives that exactly, and `==` with a number, `float`, `numeric_value` and arithmetic read it.
Between two forms `==` means equal at every q, decided exactly.

A form at a level is one number, complex when it carries powers of q. `float` and `numeric_value` give it, and
`==`, `iszero` and arithmetic are exact.

A form holds its rule, which is short and is all that deciding and evaluating read. What it prints is written from the
rule the first time it is shown, so a form costs no more than its rule until then.
"""
mutable struct QForm
    const k::Union{Nothing,Int}
    const rule::Union{Nothing,FactorialSum}             # `nothing` for a form known to be zero
    const weight::Tuple{Int,Int}                        # the weights (w, e2) of the rule
    const reduce::Bool                                  # whether writing it reduces by the relations of the level
    written::Any                                        # its `QWritten`, once it is asked for
end

"What a form prints, in the shape described at the head of this file."
struct QWritten{R<:Rational}
    k::Union{Nothing,Int}
    sign::Int
    coef::R
    out::Vector{Pair{Int,Int}}                                # ∏[n]^e
    radcoef::R                                                # r under the square root
    rad::Vector{Pair{Int,Int}}                                # ∏[n]^f under the square root
    terms::Vector{Tuple{R,Vector{Pair{Int,Int}}}}
    qpow::Int                                                 # a, doubled: q^(a/2) in front
    powers::Vector{Int}                                       # q^(b_t) of each term; empty when there are none
end
QWritten{R}(k, sign, coef, out, radcoef, rad, terms) where {R<:Rational} =
    QWritten{R}(k, sign, coef, out, radcoef, rad, terms, 0, Int[])

_written_zero(k) = QWritten{Rational{Int}}(k, 0, 0 // 1, Pair{Int,Int}[], 1 // 1, Pair{Int,Int}[],
                                           Tuple{Rational{Int},Vector{Pair{Int,Int}}}[])
_qform_zero(k) = QForm(k, nothing, (0, 0), true, _written_zero(k))

"What this form prints, written now if it has not been shown before. Two threads at worst write it twice."
function _written(f::QForm)
    w = f.written
    w === nothing || return w
    f.written = _write(f.rule, f.k, f.reduce, f.weight)
    return f.written
end

Base.iszero(f::QWritten) = iszero(f.sign)

"""
Zero or not, from the rule. One nonzero residue settles nearly every value in microseconds, at the level or, for a
generic form, at a level above its factorial arguments, where a sum vanishing at every q would vanish too. Writing the
form then proves what is left.
"""
function Base.iszero(f::QForm)
    f.written === nothing || return iszero(f.written)
    s = f.rule; w = f.weight[1]
    is_empty_sum(s) && return true
    K = f.k === nothing ? nextpow(2, max(_rule_top(s, [s.zlo, s.zhi]), 16)) : f.k
    if classify_at_level(s, K)[1] === :finite
        (iszero(w) ? _level_zero_screen(s, K, level_zero_table(K)) === false : !_weighted_residue_vanishes(s, w, K)) &&
            return false
    end
    return iszero(_written(f))
end

# ---- deciding a generic form ----

"""
Whether Σ c_t q^(b_t) ∏[n]^(e_t) vanishes at every q, for terms `(c_t, b_t, labels => e_t)`. Nothing is expanded.

Clearing denominators gives an integer polynomial in q, of degree at most D and with coefficients summing to at most H,
since q^(n−1)[n] has n coefficients equal to 1. A nonzero one cannot vanish at an integer Q > H. So this asks whether that
polynomial vanishes at Q = 2^B, a number of at most log₂H + D·B bits, by its residues at word primes. A nonzero residue
settles the question at once, and otherwise there are enough primes to make the value zero.
"""
function _generic_sum_vanishes(input::Vector{Tuple{Rational{BigInt},Int,Vector{Pair{Int,Int}}}})
    # terms that two forms share cancel here, and the bound below never sees them
    merged = Dict{Tuple{Int,Vector{Pair{Int,Int}}},Rational{BigInt}}()
    for (c, b, p) in input
        key = (b, sort(p; by = first))
        merged[key] = get(merged, key, 0 // big(1)) + c
    end
    terms = [(c, key[1], key[2]) for (key, c) in merged if !iszero(c)]
    isempty(terms) && return true
    least = Dict{Int,Int}()
    for (_, _, p) in terms, (n, e) in p
        least[n] = min(get(least, n, 0), e)
    end
    for (_, _, p) in terms, n in keys(least)
        any(x -> x.first == n, p) || (least[n] = min(least[n], 0))
    end
    den = foldl(lcm, (denominator(t[1]) for t in terms); init = big(1))
    H = big(0); lo = typemax(Int); hi = typemin(Int)
    cleared = Vector{Tuple{BigInt,Int,Dict{Int,Int}}}(undef, length(terms))
    for (i, (c, b, p)) in enumerate(terms)
        e = Dict(n => -a for (n, a) in least if !iszero(a))
        for (n, a) in p
            e[n] = get(e, n, 0) + a
        end
        filter!(x -> !iszero(x.second), e)
        C = numerator(c * den)
        H += abs(C) * prod((big(n)^a for (n, a) in e); init = big(1))
        span = sum(a * (n - 1) for (n, a) in e; init = 0)
        lo = min(lo, b - span); hi = max(hi, b + span)
        cleared[i] = (C, b, e)
    end
    B = ndigits(H; base = 2) + 1
    need = ndigits(H; base = 2) + (hi - lo) * B + 2
    N = maximum((maximum(keys(e); init = 1) for (_, _, e) in cleared); init = 1)
    bits = 0; p = UInt64(0x7fffffffffffffe7)
    while bits <= need
        while !is_prime_u64(p)
            p -= 2
        end
        m = Montgomery(p)
        Q = mont_pow(m, to_mont(m, 2), B)
        Qinv = mont_inv(m, Q)
        x = mont_add(m, Q, Qinv)
        ints = Vector{UInt64}(undef, N + 1)
        ints[1] = zero(UInt64); ints[2] = m.one              # [0] here is the true q-integer, and no term uses it
        for n in 2:N
            ints[n+1] = mont_sub(m, mont_mul(m, x, ints[n]), ints[n-1])
        end
        total = zero(UInt64)
        for (C, b, e) in cleared
            u = to_mont(m, mod(C, big(p)))
            u = mont_mul(m, u, b >= 0 ? mont_pow(m, Q, b) : mont_pow(m, Qinv, -b))
            for (n, a) in e
                u = mont_mul(m, u, mont_pow(m, ints[n+1], a))
            end
            total = mont_add(m, total, u)
        end
        iszero(total) || return false
        bits += ndigits(p; base = 2) - 1
        p -= 2
    end
    return true
end

"The terms of a generic form, without its root: the front power of q is `front`, in whole powers."
function _generic_terms(f::QWritten, front::Int, scale::Int = 1)
    out = Tuple{Rational{BigInt},Int,Vector{Pair{Int,Int}}}[]
    for (t, (c, p)) in enumerate(f.terms)
        e = Dict{Int,Int}(f.out)
        for (n, a) in p
            e[n] = get(e, n, 0) + a
        end
        push!(out, (Rational{BigInt}(scale * f.sign * f.coef * c), front + (isempty(f.powers) ? 0 : f.powers[t]), collect(e)))
    end
    return out
end

"""
Whether a rule with many terms vanishes at every q. Such a sum vanishes at every level whose h exceeds its factorial
arguments, so one nonzero residue at one such level, the usual screen, settles nearly every case in microseconds.
"""
function _generic_rule_vanishes(s::FactorialSum, w::Int, f::QWritten)
    K = nextpow(2, max(_rule_top(s, [s.zlo, s.zhi]), 16))
    if iszero(w)
        _level_zero_screen(s, K, level_zero_table(K)) === false && return false
    else
        _weighted_residue_vanishes(s, w, K) || return false
    end
    return _generic_sum_vanishes(_generic_terms(f, 0))
end

"Equality at every q. The [n] are multiplicatively independent, so what stays under the root has to agree."
function _generic_equal(f::QWritten, g::QWritten)
    (iszero(f) || iszero(g)) && return iszero(f) && iszero(g)
    (f.rad == g.rad && f.radcoef == g.radcoef && iseven(f.qpow - g.qpow)) || return false
    return _generic_sum_vanishes(vcat(_generic_terms(f, (f.qpow - g.qpow) ÷ 2), _generic_terms(g, 0, -1)))
end

"""
    classical(f::QForm) -> ClassicalResult

The value of a generic form at q = 1, exactly: a sign and the rational under the root. It compares with rationals and
integers, and converts with `Rational` when it is one.

```julia
classical(q6j(Exact(), 1, 1, 1, 1, 1, 1))          # 1//6
classical(q3j(Exact(), 2, 2, 0, 1, -1, 0))         # -√(1//5)
```
"""
function classical(f::QForm)
    f.k === nothing || throw(ArgumentError("this form is at level $(f.k), and has no value at q = 1"))
    f.rule === nothing && return zero(ClassicalResult)
    return classical_exact(f.rule)
end

# ---- the relations at a level ----

"""
The distribution relations at h, each once: `b` and `b + h/p` give the same relation, and so do `b` and `−b`. Sparse
rows `(largest label, label => exponent, index of p among the primes of h)`, in ascending order of the largest label.
"""
function _relation_generators(h::Int)
    primes = _prime_factors(h)
    gens = Tuple{Int,Vector{Pair{Int,Int}},Int}[]
    for (i, p) in enumerate(primes)
        step = h ÷ p
        for b in 1:step÷2
            e = Dict{Int,Int}()
            for (n, c) in Iterators.flatten(((p * b => 1,), ((j * step => 1) for j in 1:p-1),
                                             ((b + j * step => -1) for j in 0:p-1)))
                m, _ = _fold_label(n, h)
                m >= 2 && (e[m] = get(e, m, 0) + c)
            end
            filter!(kv -> !iszero(kv.second), e)
            isempty(e) || push!(gens, (maximum(keys(e)), collect(e), i))
        end
    end
    return sort!(gens; by = first)
end

"v − q·r, unless an entry passes the limit, in which case `v` stays as it was. Below the limit nothing can overflow."
function _row_sub!(v::Vector{T}, q::T, r::Vector{T}, tmp::Vector{T}, limit::T) where {T}
    @inbounds for i in eachindex(v)
        x = v[i] - q * r[i]
        abs(x) > limit && return false
        tmp[i] = x
    end
    copyto!(v, tmp)
    return true
end

"Reduce the labels below `c` of the row `v` by the rows that have those labels as pivots."
function _row_lower!(v::Vector{T}, c::Int, rows, tmp::Vector{T}, limit::T) where {T}
    for d in c-1:-1:2
        r = rows[d]
        (iszero(v[d]) || r === nothing) && continue
        q = Base.fld(v[d], r[d])
        iszero(q) || _row_sub!(v, q, r, tmp, limit)
    end
    return v
end

"Put a row at pivot `c` and reduce the rows above by it."
function _row_place!(v::Vector{T}, c::Int, rows, tmp::Vector{T}, limit::T) where {T}
    _row_lower!(v, c, rows, tmp, limit)
    rows[c] = v
    for d in c+1:length(rows)
        s = rows[d]
        (s === nothing || iszero(s[c])) && continue
        q = Base.fld(s[c], v[c])
        iszero(q) || (_row_sub!(s, q, v, tmp, limit) && _row_lower!(s, c, rows, tmp, limit))
    end
end

"""
Insert the generators, smallest labels first, into an echelon form that stays reduced: each row runs over the labels,
then over the primes of h for its constant. A generator that would pass `limit` is dropped, and every row kept is still
a true relation. Returns `(rows by pivot, number dropped)`.
"""
function _relation_echelon(gens, top::Int, np::Int, ::Type{T}, limit::T) where {T}
    rows = Vector{Union{Nothing,Vector{T}}}(nothing, top)
    tmp = zeros(T, top + np)
    dropped = 0
    for (_, e, i) in gens
        v = zeros(T, top + np)
        for (m, c) in e
            v[m] = c
        end
        v[top+i] = one(T)
        while v !== nothing
            c = findlast(!iszero, view(v, 1:top))
            c === nothing && break
            r = rows[c]
            if r === nothing
                v[c] < 0 && (v .= .-v)
                _row_place!(v, c, rows, tmp, limit)
                v = nothing
            elseif iszero(rem(v[c], r[c]))
                _row_sub!(v, div(v[c], r[c]), r, tmp, limit) || (dropped += 1; v = nothing)
            else
                g, a, b = gcdx(v[c], r[c])
                top_row = a .* v .+ b .* r                   # holds the gcd at c
                rest = div(r[c], g) .* v .- div(v[c], g) .* r # vanishes at c
                if maximum(abs, top_row) > limit || maximum(abs, rest) > limit
                    dropped += 1; v = nothing
                else
                    _row_place!(top_row, c, rows, tmp, limit)
                    v = rest
                end
            end
        end
    end
    return rows, dropped
end

const RelationRow = Tuple{Int,Vector{Int},Rational{BigInt},Vector{Int}}

"Above this h no table is built: it would be large, and the forms stay correct without one."
const RELATION_TABLE_LIMIT = 10^4
"Up to this h a table that lost generators in `Int` is built again in `Int128`."
const RELATION_WIDE_LIMIT = 2600
"A row with an exponent beyond this cannot make a form shorter, so it is not kept."
const RELATION_EXPONENT_LIMIT = 1 << 13

"""
The relations at h in echelon form, from the largest label down: rows `(pivot label, exponents, constant, labels with a
nonzero exponent)`, meaning ∏[n]^(exponents[n]) = constant, with the exponent of the pivot label positive and every larger
label absent. Labels are 2 ≤ n ≤ h/2.

The reduced table has large entries at very smooth h: 3·10⁹ at h = 2002. So generators are dropped past a limit and long
rows are left out. Rows are always true and a product reduced by them keeps its value, but at those levels reduction can
stop short. The square classes of `_class_rows` use the generators, and do not depend on this.
"""
function _level_relations(h::Int)
    out = RelationRow[]
    h <= RELATION_TABLE_LIMIT || return out
    top = h ÷ 2
    primes = _prime_factors(h)
    gens = _relation_generators(h)
    rows, dropped = _relation_echelon(gens, top, length(primes), Int, 1 << 31)
    if dropped > 0 && h <= RELATION_WIDE_LIMIT
        rows, dropped = _relation_echelon(gens, top, length(primes), Int128, Int128(1) << 62)
    end
    # What survives an incomplete table is mostly long. Short rows come from building it under the limit itself.
    dropped > 0 && ((rows, _) = _relation_echelon(gens, top, length(primes), Int, RELATION_EXPONENT_LIMIT))
    for c in top:-1:2
        r = rows[c]
        r === nothing && continue
        maximum(abs, r) <= RELATION_EXPONENT_LIMIT || continue
        con = 1 // big(1)
        for (i, p) in enumerate(primes)
            con *= (p // big(1))^Int(r[top+i])
        end
        e = Int.(r[1:top])
        push!(out, (c, e, con, findall(!iszero, e)))
    end
    return out
end

"The memory of a table of relations, in words."
_relations_weight(rows) = sum(r -> length(r[2]) + length(r[4]), rows; init = 1)

const LEVEL_RELATIONS = LRU{Int,Vector{RelationRow}}(maxsize = 8 * 10^6, by = _relations_weight)
_relations_at(h::Int) = get!(() -> _level_relations(h), LEVEL_RELATIONS, h)

"Reduce the exponents `e` (labels ≤ h/2) by the relations, returning how many times each row was taken out."
function _reduce_exponents!(e::Vector{Int}, rows)
    qs = zeros(Int, length(rows))
    for (r, (piv, row, _, labels)) in enumerate(rows)
        q = fld(e[piv], row[piv])
        iszero(q) && continue
        @inbounds for i in labels
            e[i] = Base.Checked.checked_sub(e[i], Base.Checked.checked_mul(q, row[i]))
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
    qform(rule; k = nothing, reduce = true, weight = (0, 0)) -> QForm
    qform(symbol, labels...; k = nothing, reduce = true) -> QForm
    qform(q6j(Symbolic(), labels...); k = nothing, reduce = true) -> QForm

A value written in quantum integers `[n]`, the form recoupling theory is written in. Without `k` it holds at
every generic q, and `q = 1` turns it into integers. With `k` the labels fold to `n ≤ (k+2)/2` and vanishing
factors leave their limit. With `reduce = true` the distribution relations of the level are applied too: terms
can then merge or cancel, and a zero they reveal is returned as zero. The shorter of the two forms is returned.

`qcg` and `q3j` carry powers of q, printed in front and on each term. At a level these reduce by q^h = −1 and the value
is complex, with `q` read as e^{iπ/(k+2)}. A rule with weights takes them as `weight = (w, e2)`, as `qcg_rule` returns them.

With `reduce = true`, an exact check also finds cancellations that need additive relations, generic or at a level. With
`reduce = false`, `iszero_at` can decide cancellations left in the expression. Two equal nonzero values can still print
differently, and `==` decides their equality. Numerical evaluation of a rule-backed form reuses the original factorial
rule and its stable numerical kernels.

`qform(f, labels...; k)` is what `f(Exact(k), labels...)` returns; this function also takes `reduce` and bare rules.

```julia
qform(q6j, 1, 1, 1, 1, 1, 1)          # ([5] − 1)/[4]!, which is 1/6 at q = 1
qform(q6j, 1, 1, 1, 1, 1, 1; k = 3)   # −1/[2]²
qform(q6j, 1, 1, 1, 2, 1, 1; k = 5)   # [2]/[3]²
qform(qcg, 1//2, 1//2, 1//2, -1//2, 1, 0)   # q^(−1/2) √(1/[2])
```
"""
function qform(s::FactorialSum; k = nothing, reduce::Bool = true, weight::Tuple{Int,Int} = (0, 0))
    k === nothing && return QForm(nothing, s, weight, reduce, nothing)
    classify_at_level(s, Int(k))[1] === :pole && throw(ArgumentError("the rule is singular at level $k"))
    return QForm(Int(k), s, weight, reduce, nothing)
end

function _write(s::FactorialSum, k, reduce::Bool, weight::Tuple{Int,Int})
    try
        return _qform(s, k, reduce, weight, Rational{Int})
    catch e
        (e isa OverflowError || e isa InexactError) || rethrow()
        return _qform(s, k, reduce, weight, Rational{BigInt})
    end
end

"How long a form reads: its terms first, then its exponents and the digits of its numbers."
function _qform_cost(f::QWritten)
    digits(c) = ndigits(numerator(c)) + ndigits(denominator(c)) - 2
    exps(p) = sum(abs(e) for (_, e) in p; init = 0)
    return (length(f.terms), exps(f.out) + exps(f.rad) + sum(exps(p) + digits(c) for (c, p) in f.terms; init = 0) +
                             digits(f.coef) + digits(f.radcoef))
end

"The largest factorial argument of a rule over the terms `zs`, at least 2."
function _rule_top(s::FactorialSum, zs::Vector{Int})
    N = 2
    for (n, _) in s.pre
        N = max(N, Int(n))
    end
    for f in s.fac, z in (first(zs), last(zs))
        N = max(N, _arg(f, z))
    end
    return N
end

"The prefactor squared, as doubled exponents of [1], …, [N]."
function _prefactor_exponents(s::FactorialSum, N::Int)
    pre2 = zeros(Int, N)
    for (n, c) in s.pre
        n >= 2 || continue
        pre2[2] += (s.sqrt_pre ? 1 : 2) * c
        n < N && (pre2[Int(n)+1] -= (s.sqrt_pre ? 1 : 2) * c)
    end
    for i in 3:N
        pre2[i] += pre2[i-1]
    end
    return pre2
end

"""
A rule folded to its level: `(zs, pre, prerat, es, cs)`, the contributing terms, the squared prefactor as exponents
and a rational, and each term as exponents and a rational. `nothing` when no term contributes.
"""
function _level_exponents(s::FactorialSum, level::Int, ::Type{R}) where {R}
    is_empty_sum(s) && return nothing
    st, segs = classify_at_level(s, level)
    st === :pole && throw(ArgumentError("the rule is singular at level $level"))
    st === :finite || return nothing
    zs = Int[z for seg in segs for z in seg]
    N = _rule_top(s, zs); h = level + 2
    raw = zeros(Int, N)
    lab, neg = _fold_table(N, h)
    pre, prerat, ps, nhp = _fold_exponents(_prefactor_exponents(s, N), h, lab, neg, R)
    ps < 0 && throw(ArgumentError("the radicand is negative at level $level"))
    es = Vector{Vector{Int}}(undef, length(zs)); cs = Vector{R}(undef, length(zs))
    for (t, z) in enumerate(zs)
        e, r, ts, nh = _fold_exponents(_term_exponents!(raw, s, z), h, lab, neg, R)
        nhp + 2nh == 0 || throw(ArgumentError("the terms of this rule have different valuations at level $level"))
        c = ts * r
        s.alternating && isodd(z) && (c = -c)
        es[t] = e; cs[t] = c
    end
    return zs, pre, prerat, es, cs
end

"The powers of q on the terms `zs`, in [0, h) at a level, with the signs that q^h = −1 leaves on `cs`."
function _term_powers(w::Int, zs::Vector{Int}, cs::Vector{R}, h) where {R}
    iszero(w) && return nothing
    ps = w .* zs
    h === nothing && return ps
    for t in eachindex(ps)
        n, a = fldmod(ps[t], h)
        ps[t] = a
        isodd(n) && (cs[t] = -cs[t])
    end
    return ps
end

function _qform(s::FactorialSum, k, reduce::Bool, weight::Tuple{Int,Int}, ::Type{R}) where {R}
    w, e2 = weight
    if k === nothing
        is_empty_sum(s) && return _written_zero(nothing)
        zs = collect(s.zlo:s.zhi)
        raw = zeros(Int, _rule_top(s, zs))
        es = [copy(_term_exponents!(raw, s, z)) for z in zs]
        cs = R[s.alternating && isodd(z) ? -one(R) : one(R) for z in zs]
        f = _qform_assemble(nothing, Int(s.sign0), _prefactor_exponents(s, length(raw)), one(R), es, cs,
                            _term_powers(w, zs, cs, nothing), e2)
        # a sum can vanish at every q although no two of its terms merge
        reduce && length(f.terms) > 1 && _generic_rule_vanishes(s, w, f) && return _written_zero(nothing)
        return f
    end
    level = Int(k); h = level + 2
    data = _level_exponents(s, level, R)
    data === nothing && return _written_zero(level)
    zs, pre, prerat, es, cs = data
    ps = _term_powers(w, zs, cs, h)
    plain = _qform_assemble(level, Int(s.sign0), pre, prerat, es, cs, ps, e2)
    reduce || return plain
    iszero(plain) && return plain
    # terms can cancel by additive relations that no rewriting of the products shows
    length(plain.terms) > 1 && (iszero(w) ? _level_sum_zero(s, zs, es, cs, h) : _weighted_level_zero(s, w, level)) &&
        return _written_zero(level)
    rows = _relations_at(h)
    isempty(rows) && return plain                            # h prime: the folded labels are independent
    # Exponents first. Constants are only worth computing for a form that merges terms or reads shorter.
    es2 = copy.(es); pre2r = copy(pre)
    qs, qp = try
        [_reduce_exponents!(e, rows) for e in es2], _reduce_exponents!(pre2r, rows)
    catch err
        err isa OverflowError || rethrow()
        return plain                                         # exponents this long never read shorter
    end
    _exponent_cost(pre2r, es2) < _exponent_cost(pre, es) || return plain
    cs2 = [cs[t] * _relation_constant(rows, qs[t], R) for t in eachindex(cs)]
    reduced = _qform_assemble(level, Int(s.sign0), pre2r, prerat * _relation_constant(rows, qp, R), es2, cs2, ps, e2)
    iszero(reduced) && return reduced
    return _qform_cost(reduced) < _qform_cost(plain) ? reduced : plain
end

"Merge equal terms, move what every term shares in front, and take even exponents out of the root."
function _qform_assemble(level, sign::Int, pre::Vector{Int}, prerat::R, es::Vector{Vector{Int}}, cs::Vector{R},
                         ps = nothing, e2::Int = 0) where {R}
    seen = Dict{Vector{Int},Int}(); keep = Int[]; coef = R[]
    for t in eachindex(es)
        i = get!(seen, ps === nothing ? es[t] : vcat(es[t], ps[t]), length(keep) + 1)   # powers must agree to merge
        if i > length(keep)
            push!(keep, t); push!(coef, cs[t])
        else
            coef[i] += cs[t]
        end
    end
    live = [i for i in eachindex(keep) if !iszero(coef[i])]
    isempty(live) && return _written_zero(level)
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
    sign = scale < 0 ? -sign : sign
    ps === nothing && e2 == 0 && return QWritten{R}(level, sign, c, out, under, rad, terms)
    least = ps === nothing ? 0 : minimum(ps[keep[i]] for i in ts)
    powers = ps === nothing ? Int[] : [ps[keep[i]] - least for i in ts]
    qpow = e2 + 2least
    if level !== nothing                                     # q^h = −1, nearest to 1 in front
        h = level + 2
        n = fld(qpow + h, 2h)
        qpow -= 2h * n
        isodd(n) && (sign = -sign)
    end
    return QWritten{R}(level, sign, c, out, under, rad, terms, qpow, powers)
end

qform(v::SymbolicValue; kw...) = qform(v.rule; kw...)

function qform(f::Function, args...; k = nothing, reduce::Bool = true)
    (f === qcg || f === q3j) && return _qform_weighted(f, k, reduce, args...)
    k !== nothing && !_admissible_at(f, Int(k), args...) && return _qform_zero(Int(k))
    return qform(_rule_for(f, args...); k = k, reduce = reduce)
end

function _qform_weighted(s::FactorialSum, weight, admissible, k, reduce::Bool)
    k !== nothing && !admissible(Int(k)) && return _qform_zero(Int(k))
    return qform(s; k = k, reduce = reduce, weight = weight)
end
function _qform_weighted(::typeof(qcg), k, reduce::Bool, j1::Spin, m1::Spin, j2::Spin, m2::Spin, j::Spin, m::Spin = m1 + m2)
    J1, M1, J2, M2, J, M = doubled(j1, m1, j2, m2, j, m)
    return _qform_weighted(qcg_rule(J1, M1, J2, M2, J, M)..., kk -> _qδ(J1, J2, J, kk), k, reduce)
end
_qform_weighted(::typeof(qcg), k, reduce::Bool, a::Tuple{Spin,Spin}, b::Tuple{Spin,Spin}, c::Tuple{Spin,Spin}) =
    _qform_weighted(qcg, k, reduce, a..., b..., c...)
function _qform_weighted(::typeof(q3j), k, reduce::Bool, j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1 - m2)
    J1, J2, J3, M1, M2, M3 = doubled(j1, j2, j3, m1, m2, m3)
    return _qform_weighted(q3j_rule(J1, J2, J3, M1, M2, M3)..., kk -> _qδ(J1, J2, J3, kk), k, reduce)
end

# ---- numbers ----

"The form as written, from `qn(n)` = [n] and `phase(a)` = q^(a/2)."
function _qform_value(f::QWritten, qn, phase = a -> 1)
    iszero(f) && return zero(qn(2))
    mono(p) = (v = one(qn(2)); for (n, e) in p; v *= qn(n)^e; end; v)
    total = zero(qn(2))
    for (t, (c, p)) in enumerate(f.terms)
        u = oftype(real(total), numerator(c)) / oftype(real(total), denominator(c)) * mono(p)
        total += isempty(f.powers) ? u : u * phase(2 * f.powers[t])
    end
    w = oftype(qn(2), numerator(f.coef)) / oftype(qn(2), denominator(f.coef))
    u = oftype(qn(2), numerator(f.radcoef)) / oftype(qn(2), denominator(f.radcoef))
    return f.sign * w * mono(f.out) * sqrt(u * mono(f.rad)) * total * phase(f.qpow)
end

_qform_value(f::QForm, args...) = _qform_value(_written(f), args...)

function (f::QForm)(q::Real)
    f.k === nothing || throw(ArgumentError("this form is at level $(f.k); use `numeric_value` or `float`"))
    q > 0 || throw(DomainError(q, "a QForm is evaluated at positive real q"))
    Q = float(q)
    (f.rule === nothing || is_empty_sum(f.rule)) && return zero(Q)
    f.weight == (0, 0) && return qeval(f.rule; q = Q, T = typeof(Q))
    return isone(Q) ? classical_value(f.rule, typeof(Q)) : analytic_value(f.rule, Q, typeof(Q); weight = f.weight)
end

"The form as written at its level, in the current precision."
function _qform_level(f::QForm)
    f.k === nothing && throw(ArgumentError("a generic form is evaluated with `f(q)`"))
    h = f.k + 2
    return _qform_value(_written(f), n -> sinpi(BigFloat(n) / h) / sinpi(inv(BigFloat(h))), a -> cispi(BigFloat(a) / (2h)))
end

"Round a value at a level to the current precision. Weighted values are complex."
_round_level(v::Real) = BigFloat(v)
_round_level(v::Complex) = complex(BigFloat(real(v)), BigFloat(imag(v)))

function numeric_value(f::QForm; bits::Integer = precision(BigFloat))
    bits > 0 || throw(ArgumentError("bits must be positive"))
    f.k === nothing && return numeric_value(classical(f); bits = bits)
    iszero(f) && return setprecision(() -> zero(BigFloat), BigFloat, Int(bits))
    v = setprecision(BigFloat, Int(bits) + 32) do
        f.weight == (0, 0) && return qeval(f.rule; k = f.k, T = BigFloat)
        return _weighted_level_value(f.rule, f.weight..., f.k, Complex{BigFloat})
    end
    return setprecision(() -> _round_level(v), BigFloat, Int(bits))
end
function Base.float(f::QForm)
    v = numeric_value(f; bits = 128)
    return v isa Complex ? ComplexF64(v) : Float64(v)
end
Base.Float64(f::QForm) = Float64(numeric_value(f; bits = 128))
Base.BigFloat(f::QForm) = numeric_value(f)
Base.big(f::QForm) = BigFloat(f)
Base.Rational(f::QForm) = Rational(classical(f))
Base.isapprox(f::QForm, g::QForm; kw...) = isapprox(float(f), float(g); kw...)
Base.isapprox(f::QForm, x::Number; kw...) = isapprox(float(f), x; kw...)
Base.isapprox(x::Number, f::QForm; kw...) = isapprox(x, float(f); kw...)
Base.zero(f::QForm) = _qform_zero(f.k)
Base.one(f::QForm) = qform(FactorialSum(0:0); k = f.k)

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

"q^(a/2) as text, from the doubled power `a`."
function _q_phase(a::Int)
    iszero(a) && return ""
    a == 2 && return "q"
    return iseven(a) ? "q" * to_superscript(a ÷ 2) : "q^(" * (a < 0 ? "−" : "") * string(abs(a)) * "/2)"
end

function _q_fraction(c::Rational, p::Vector{Pair{Int,Int}}, inner::String = "", phase::String = "")
    up, down = _q_product(p)
    a = numerator(c); b = denominator(c)
    top = (isone(a) && !(isempty(up) && isempty(inner) && isempty(phase)) ? "" : string(a)) * phase * up * inner
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

function _q_sum(terms; most::Int = QFORM_SHOWN_TERMS, powers::Vector{Int} = Int[])
    length(terms) == 1 && isone(first(terms)[1]) && isempty(first(terms)[2]) && return ""
    shown = first(terms, most)
    parts = String[]
    for (i, (c, p)) in enumerate(shown)
        t = _q_fraction(abs(c), p, "", isempty(powers) ? "" : _q_phase(2 * powers[i]))
        push!(parts, i == 1 ? (c < 0 ? "−" : "") * t : (c < 0 ? " − " : " + ") * t)
    end
    extra = length(terms) - length(shown)
    extra > 0 && push!(parts, " + … ($(length(terms)) terms)")
    return length(terms) == 1 ? join(parts) : "(" * join(parts) * ")"
end

"The form as text, without its level, showing at most `most` terms of its sum."
function _qform_body(f::QWritten; most::Int = QFORM_SHOWN_TERMS)
    iszero(f) && return "0"
    text = (f.sign < 0 ? "−" : "") * _q_fraction(f.coef, f.out, _q_sum(f.terms; most = most, powers = f.powers), _q_phase(f.qpow))
    isempty(f.rad) && isone(f.radcoef) && return _unwrap(text)
    root = "√(" * _q_fraction(f.radcoef, f.rad) * ")"
    text == "1" && return root
    text == "−1" && return "−" * root
    return text * " " * root
end

_qform_body(f::QForm; most::Int = QFORM_SHOWN_TERMS) = _qform_body(_written(f); most = most)

function Base.show(io::IO, f::QForm)
    w = _written(f)
    print(io, _qform_body(w))
    f.k === nothing || print(io, "   (k = ", f.k, isempty(w.powers) && iszero(w.qpow) ? ")" : ", q = e^{iπ/$(f.k + 2)})")
end

"The longest classical value that is printed beside a form."
const QFORM_SHOWN_VALUE = 40

"The form, with the number it stands for beside it: its value at q = 1, or its value at its level."
function Base.show(io::IO, ::MIME"text/plain", f::QForm)
    show(io, f)
    iszero(_written(f)) && return
    if f.k === nothing
        v = replace(sprint(show, classical(f)), "-" => "−", "//1" => "", "//" => "/")
        length(v) <= QFORM_SHOWN_VALUE && print(io, "   (", endswith(v, "/1") ? v[1:end-2] : v, " at q = 1)")
    else
        v = float(f)
        print(io, "   ≈ ", replace(string(v isa Complex ? complex(round(real(v); sigdigits = 7), round(imag(v); sigdigits = 7)) : round(v; sigdigits = 7)), "-" => "−"))
    end
end


"A sum of labelled terms over a common factor, as text."
function _q_sum_str(terms; most::Int)
    isempty(terms) && return "0"
    top = big(0); bot = big(1)
    for (c, _) in terms
        top = gcd(top, numerator(c)); bot = lcm(bot, denominator(c))
    end
    scale = top // bot
    return _q_fraction(scale, Pair{Int,Int}[], _q_sum([(c / scale, p) for (c, p) in terms]; most = most))
end
