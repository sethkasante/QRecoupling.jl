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
factorial rule instead, whose `.dcr` is the cyclotomic monomial this function used to return.

```julia
qint(5)                 # 5.0
qint(5; k = 10)         # [5] at q = exp(iπ/12)
qint(5; q = 0.7)
qint(Exact(10), 5)      # exact, in radicals when they exist
```
"""
function qint(n::Integer, p::Integer = 1; k = nothing, q = nothing, exact::Bool = false,
              T::Type = Float64, workspace = nothing)
    q = _evaluation_q(k, q, exact)
    return _product_value(_qint_pairs(Int(n), Int(p)), k, q, exact, T; workspace = workspace)
end

"""
    qfact(n, p = 1; k, q, exact, T)

The q-factorial `([n]!)^p = ([1][2]⋯[n])^p`, classical by default, with the same keywords as [`qint`](@ref).
`n < 0` is zero and `[0]! = [1]! = 1`.
"""
function qfact(n::Integer, p::Integer = 1; k = nothing, q = nothing, exact::Bool = false,
               T::Type = Float64, workspace = nothing)
    q = _evaluation_q(k, q, exact)
    return _product_value(_qfact_pairs(Int(n), Int(p)), k, q, exact, T; workspace = workspace)
end

"""
    qbinomial(n, m; k, q, exact, T)

The Gaussian binomial `[n]! / ([m]! [n−m]!)`, classical by default, with the same keywords as
[`qint`](@ref). Outside `0 ≤ m ≤ n` it is zero.
"""
function qbinomial(n::Integer, m::Integer; k = nothing, q = nothing, exact::Bool = false,
                   T::Type = Float64, workspace = nothing)
    q = _evaluation_q(k, q, exact)
    return _product_value(_qbinomial_pairs(Int(n), Int(m)), k, q, exact, T; workspace = workspace)
end

"The rule itself, which is what `Symbolic()` asks for; the `Symbolic` methods live in `targets.jl`."
_product_symbolic(pairs) = SymbolicValue(pairs === nothing ? _ZERO_PRODUCT_RULE : _qfact_rule(pairs))
