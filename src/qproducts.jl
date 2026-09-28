# ---------------------------------------------------------------------------------
#  q-integers, q-factorials and q-binomials, on the same footing as the symbols
#
#  These are not sums: each is a bare product `Π [n]!^{c}`, which is a `FactorialSum` whose range is a
#  single point and whose factor list is empty. Saying so routes them through the same evaluator as
#  `q6j` — the level kernel with its running bound, the classical integer path, the analytic ladder at
#  generic q, and the real basis when an exact value is asked for — instead of through the monomial
#  projection they used to take. One consequence is worth stating: `qint(5)` is now a *number*, the
#  classical value, exactly as `q6j(…)` is. The cyclotomic monomial is still there, as
#  `qint(Symbolic(), 5).dcr`.
# ---------------------------------------------------------------------------------

"""
A product of q-factorials `Π [n]!^{c}` as a factorial rule: one term, a prefactor, and nothing else.

Every route the symbols use accepts this, so a q-binomial gets the level kernel's certified bound and the
generic-q ladder for free rather than having its own arithmetic.
"""
_qfact_rule(pairs) = FactorialSum(0:0; prefactor = pairs)

"The empty rule, whose value is the zero of whichever carrier is asked for."
const _ZERO_PRODUCT_RULE = FactorialSum((), false, (), false, Int8(0), 0, -1)

"""
    _product_value(pairs, k, q, exact, T; workspace) -> value

`Π [n]!^{c}` wherever the symbols evaluate. `pairs` may be empty, which is the value 1, and `nothing`
stands for the empty product, which is 0 in whatever carrier the target names.

**Numerics take the monomial route, not the rule.** Routing a bare product through the factorial-rule
evaluator was tried and is 2–30× *slower*: a product has no sum to certify, so all the rule machinery
adds is about 1.3 µs of fixed cost per call against 0.04–0.6 µs for a handful of table lookups
(`dev/results/v04_cleanup.md`). What the rule is for here is the *other* answers — `Symbolic()` hands
back the rule, and an exact level value is built from its q-factorial content in the real basis — so the
interface is the one the symbols have while the arithmetic is the cheap one.
"""
function _product_value(pairs, k, q, exact::Bool, ::Type{T}; workspace = nothing) where {T}
    if exact && !isnothing(k)
        kk = Int(k); kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
        pairs === nothing && return zero(ExactX, kk)
        return _qfact_exactx(pairs, kk)
    end
    if _is_classical(q) && !exact
        # At q = 1 the rule's kernel works in integer ratios, so `qint(5)` is 5.0 and not
        # 4.999999999999999 — which is what the monomial route returns, since it multiplies floating
        # values of Φ_d(1). Eight times the cost of that route and still a microsecond; for the one
        # place a user checks a package's arithmetic by eye, the right trade.
        return classical_value(pairs === nothing ? _ZERO_PRODUCT_RULE : _qfact_rule(pairs), T;
                               workspace = workspace)
    end
    return qeval(_pairs_mono(pairs); k = k, q = q, exact = exact, T = T)
end

# ---------------------------------------------------------------------------------
#  Numbers straight from the tables, for Float64 targets
#
#  A q-integer is one table entry, a q-factorial one split entry, a q-binomial three. Reading them is what
#  the symbols do already, so the cheap route is also the accurate one: the level table holds [n] correctly
#  rounded, and the split factorials hold [n]! as a correctly rounded mantissa and an exact exponent. The
#  monomial route this replaces built a cyclotomic monomial on every call (7–10 allocations, and 21 with
#  1.2 µs for a classical `qint(5)`), multiplied floating values of Φ_d, and near q = 1 lost up to 4e-10.
#  Anything these paths do not cover — other element types, exact targets, a level factorial that reaches
#  [h] = 0 — returns `nothing` and takes the route it always took.
# ---------------------------------------------------------------------------------

"[n] at q = e^{iπ/(k+2)}, correctly rounded from the level table, for every integer n by periodicity."
@inline function _qint_level(n::Int, k::Int)
    h = k + 2
    r = mod(n, 2h)
    (r == 0 || r == h) && return 0.0
    tab = qint_tables(Float64, k)
    return r < h ? @inbounds(tab.q[r]) : -@inbounds(tab.q[2h - r])
end

"""
[n] at a real positive or complex q, as q^{1−n}(q^{2n} − 1)/(q² − 1).

When q^{2n} is near 1 that difference is the whole value and cancels in floating point, so it is taken
from `expm1` of 2n·log q instead; otherwise the integer powers are more accurate than exp(n·log q), whose
error grows like n·|log q|·u. `nothing` where neither form applies (a negative real q).
"""
function _qint_analytic(n::Int, q::Union{Float64,ComplexF64})
    q isa Float64 && !(q > 0) && return nothing
    L = log(q)
    if abs(2n * L) < 1
        return q^(1 - n) * (expm1(2n * L) / expm1(2L))
    end
    # Written so that the large power is the one outside: for |q| > 1 the mirrored form, since
    # q^{1−n}·q^{2n} would be 0·∞ long before the value itself overflows.
    abs(q) <= 1 && return q^(1 - n) * ((q^(2n) - 1) / ((q - 1) * (q + 1)))
    return q^(n - 1) * ((1 - inv(q^(2n))) / (((q - 1) * (q + 1)) / (q * q)))
end

"Π [n]!^c from the split tables at a level; `nothing` if a factorial reaches [h] = 0."
function _qfacts_level(k::Int, pairs::NTuple{N,Tuple{Int,Int}}) where {N}
    tab = qint_tables(Float64, k)
    m, e = 1.0, 0
    for (n, c) in pairs
        n <= k + 1 || return nothing
        n >= 1 && ((m, e) = _split_mul(m, e, tab, n, c))
    end
    return ldexp(m, e)
end

"""
The Float64 value of a q-integer power, q-factorial power or q-binomial, or `nothing` to take the general
route. `what` is `:int`, `:fact` or `:binomial`; `a`, `b` are (n, p) or (n, m).
"""
function _qnumber_float(what::Symbol, a::Int, b::Int, k, q)
    # the same conventions as the pair constructors: negative n is zero, [0] = [1] = 1, [0]! = [1]! = 1
    if what === :binomial
        (b < 0 || b > a) && return 0.0
        (b == 0 || b == a) && return 1.0
    else
        a < 0 && return 0.0
        (a <= 1 || b == 0) && return 1.0
    end
    if _is_classical(q)
        if what === :int
            return Float64(a)^b
        elseif what === :fact
            tab = classical_tables(Float64, a)
            m, e = _split_mul(1.0, 0, tab, a, b)
            return ldexp(m, e)
        else
            a <= 66 && return Float64(binomial(a, b))          # exact in Int64, rounded once
            tab = classical_tables(Float64, a)
            m, e = _split_mul(1.0, 0, tab, a, 1)
            m, e = _split_mul(m, e, tab, b, -1)
            m, e = _split_mul(m, e, tab, a - b, -1)
            return ldexp(m, e)
        end
    elseif k !== nothing
        kk = Int(k)
        what === :int && return _qint_level(a, kk)^b
        what === :fact && return _qfacts_level(kk, ((a, b),))
        return _qfacts_level(kk, ((a, 1), (b, -1), (a - b, -1)))
    elseif q isa Union{Float64,ComplexF64}
        what === :int || return nothing                      # products of many [m] keep the old route
        v = _qint_analytic(a, q)
        return v === nothing ? nothing : v^b
    end
    return nothing
end

"The cyclotomic monomial of `Π [n]!^{c}` — what the numeric projections consume."
function _pairs_mono(pairs)
    pairs === nothing && return ZERO_MONOMIAL
    isempty(pairs) && return ONE_MONOMIAL
    buf = CycloBuffer(max(maximum(first, pairs), 1))
    for (m, c) in pairs
        add_qfact!(buf, Int(m), Int(c))
    end
    return snapshot(buf)
end

"`nothing` for a vanishing product, `Pair[]` for one that is identically 1, else the q-factorial content."
function _qint_pairs(n::Int, p::Int)
    n < 0 && return nothing
    (n <= 1 || p == 0) && return Pair{Int,Int}[]
    return [n => p, n - 1 => -p]
end
function _qfact_pairs(n::Int, p::Int)
    n < 0 && return nothing
    (n <= 1 || p == 0) && return Pair{Int,Int}[]
    return [n => p]
end
function _qbinomial_pairs(n::Int, m::Int)
    (m < 0 || m > n) && return nothing
    (m == 0 || m == n) && return Pair{Int,Int}[]
    return [n => 1, m => -1, n - m => -1]
end

"""
    qint(n, p = 1; k, q, exact, T)

The q-integer `[n]^p`, classical by default — `[n] = (qⁿ − q⁻ⁿ)/(q − q⁻¹)`, with the Turaev–Viro
convention `[0] = 1`. `n < 0` is zero.

The keywords are the ones every symbol takes: `k` for a level, `q` for a parameter, `exact = true` for
the exact classical rational or, with `k`, the real-basis level value. `qint(Symbolic(), n)` returns the
rule-backed symbolic x-form instead; its `.dcr` property constructs a one-term compatibility DCR.

```julia
qint(5)                 # 5.0
qint(5; k = 10)         # [5] at q = exp(iπ/12)
qint(5; q = 0.7)
qint(Exact(10), 5)      # exact, as a polynomial in x; `radical` for the surd
```
"""
Base.@constprop :aggressive @inline function qint(n::Integer, p::Integer = 1; k = nothing, q = nothing, exact::Bool = false,
              T::Type{TT} = Float64, workspace = nothing) where {TT}
    q = _evaluation_q(k, q, exact)
    if !exact && T === Float64 && !(k isa AbstractVector)
        v = _qnumber_float(:int, Int(n), Int(p), k, q)
        v === nothing || return v
        # a Float64 target lands on a Float64 or, at complex q, a ComplexF64: saying so keeps the
        # return unboxed, since the general route below is not inferable
        return _product_value(_qint_pairs(Int(n), Int(p)), k, q, exact, T; workspace = workspace)::Union{Float64,ComplexF64}
    end
    return _product_value(_qint_pairs(Int(n), Int(p)), k, q, exact, T; workspace = workspace)
end

"""
    qfact(n, p = 1; k, q, exact, T)

The q-factorial `([n]!)^p = ([1][2]⋯[n])^p`, classical by default, with the same keywords as [`qint`](@ref).
`n < 0` is zero and `[0]! = [1]! = 1`.
"""
Base.@constprop :aggressive @inline function qfact(n::Integer, p::Integer = 1; k = nothing, q = nothing, exact::Bool = false,
               T::Type{TT} = Float64, workspace = nothing) where {TT}
    q = _evaluation_q(k, q, exact)
    if !exact && T === Float64 && !(k isa AbstractVector)
        v = _qnumber_float(:fact, Int(n), Int(p), k, q)
        v === nothing || return v
        # a Float64 target lands on a Float64 or, at complex q, a ComplexF64: saying so keeps the
        # return unboxed, since the general route below is not inferable
        return _product_value(_qfact_pairs(Int(n), Int(p)), k, q, exact, T; workspace = workspace)::Union{Float64,ComplexF64}
    end
    return _product_value(_qfact_pairs(Int(n), Int(p)), k, q, exact, T; workspace = workspace)
end

"""
    qbinomial(n, m; k, q, exact, T)

The Gaussian binomial `[n]! / ([m]! [n−m]!)`, classical by default, with the same keywords as
[`qint`](@ref). Outside `0 ≤ m ≤ n` it is zero.
"""
Base.@constprop :aggressive @inline function qbinomial(n::Integer, m::Integer; k = nothing, q = nothing, exact::Bool = false,
                   T::Type{TT} = Float64, workspace = nothing) where {TT}
    q = _evaluation_q(k, q, exact)
    if !exact && T === Float64 && !(k isa AbstractVector)
        v = _qnumber_float(:binomial, Int(n), Int(m), k, q)
        v === nothing || return v
        # a Float64 target lands on a Float64 or, at complex q, a ComplexF64: saying so keeps the
        # return unboxed, since the general route below is not inferable
        return _product_value(_qbinomial_pairs(Int(n), Int(m)), k, q, exact, T; workspace = workspace)::Union{Float64,ComplexF64}
    end
    return _product_value(_qbinomial_pairs(Int(n), Int(m)), k, q, exact, T; workspace = workspace)
end

"The rule itself, which is what `Symbolic()` asks for; the `Symbolic` methods live in `targets.jl`."
_product_symbolic(pairs) = SymbolicValue(pairs === nothing ? _ZERO_PRODUCT_RULE : _qfact_rule(pairs))
