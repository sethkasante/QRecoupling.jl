# --------------------------------------------
#  Structural queries: zeros, poles and level spectra
#
#  These answer questions about a symbol without computing its value. Valuations give poles and
#  all-terms-vanishing zeros in closed form; the modular test settles cancellation zeros. Microseconds
#  either way, which is what makes sweeping whole families practical.
# --------------------------------------------

_rule_for(::typeof(q6j), args...) = sixj_sum(doubled(args...)...)
_rule_for(::typeof(fsymbol), args...) = fsymbol_sum(doubled(args...)...)
_rule_for(::typeof(gsymbol), args...) = gsymbol_sum(doubled(args...)...)
_rule_for(::typeof(q3j), j1, j2, j3, m1, m2, m3 = -m1 - m2) = threej_sum(doubled(j1, j2, j3, m1, m2, m3)...)

# Labels outside the level's admissible set carry no representation, so the symbol is zero by convention
# rather than singular. The queries follow the same convention as the symbol functions.
_admissible_at(::typeof(q6j), k::Int, js...) = _qδtet(doubled(js...)..., k)
_admissible_at(::typeof(gsymbol), k::Int, js...) = _qδtet(doubled(js...)..., k)
_admissible_at(::typeof(fsymbol), k::Int, js...) = _qδtet(doubled(js...)..., k)
_admissible_at(::typeof(q3j), k::Int, j1, j2, j3, m1, m2, m3 = -m1 - m2) = _qδ(doubled(j1, j2, j3)..., k)

"""
    iszero_at(k, j1, j2, j3, j4, j5, j6) -> Bool
    iszero_at(symbol, k, labels...) -> Bool
    iszero_at(k, labels::AbstractVector) -> BitVector

Whether a symbol is *exactly* zero at level `k`, including zeros that come from cancellation between
finite terms. `symbol` is one of `q6j`, `q3j`, `fsymbol`, `gsymbol` (default `q6j`). Given a collection of
label tuples, the test runs over the batch with shared tables and threads.

The answer is exact up to a false-positive probability of about 2⁻¹²⁴ from the modular test.

```julia
iszero_at(20, 5, 5, 5, 5, 5, 5)          # true: a cancellation zero
iszero_at(q3j, 10, 1, 1, 1, 1, -1, 0)
count(iszero_at(6, all_6j(k = 6)))       # how many symbols vanish at level 6
```
"""
iszero_at(k::Integer, args::Spin...) = iszero_at(q6j, k, args...)

function iszero_at(f::Function, k::Integer, args::Spin...)
    k = Int(k)
    _admissible_at(f, k, args...) || return true
    return is_zero_at_level(_rule_for(f, args...), k)
end

iszero_at(k::Integer, labels::AbstractVector; threads = nothing) = iszero_at(q6j, k, labels; threads = threads)

function iszero_at(f::Function, k::Integer, labels::AbstractVector; threads = nothing)
    k = Int(k)
    n = length(labels)
    out = falses(n)
    n == 0 && return out
    ztab = level_zero_table(k)                  # built once, read-only afterwards
    work = function (rng)
        for i in rng
            l = labels[i]
            out[i] = !_admissible_at(f, k, l...) || is_zero_at_level(_rule_for(f, l...), k, ztab)
        end
    end
    nw = _nworkers(n, threads)
    if nw == 1
        work(1:n)
    else
        @sync for c in _chunks(n, nw)
            Threads.@spawn work(c)
        end
    end
    return out
end

"""
    issingular_at(k, j1, ..., j6) -> Bool
    issingular_at(symbol, k, labels...) -> Bool

Whether these labels are singular at level `k`: the level-k rule has a contributing term of negative
valuation, i.e. a q-factorial in a denominator that vanishes. Decided from the valuations alone, with no
arithmetic. `symbol` is one of `q6j`, `q3j`, `fsymbol`, `gsymbol` (default `q6j`).

Singular labels are exactly the ones a level cannot represent — for the 6j, those with a triangle sum above
2k. Labels that *are* admissible at the level are never singular (their prefactor and term valuations are both
non-negative), so this query is false throughout the admissible set; it answers a question about the labels,
not about the value the symbol functions return. Those follow the Turaev–Viro convention and give `0` outside
the admissible set, and `level_spectrum` reports such levels as `:inadmissible`.

```julia
issingular_at(12, 10, 10, 10, 10, 10, 10)     # true: triangle sums 30 > 2k
issingular_at(30, 10, 10, 10, 10, 10, 10)     # false: admissible from k = 30
issingular_at(4, 1, 1, 5, 1, 1, 5)            # false: no triangle, so no representation to be singular
```
"""
issingular_at(k::Integer, args::Spin...) = issingular_at(q6j, k, args...)

function issingular_at(f::Function, k::Integer, args::Spin...)
    k = Int(k)
    return classify_at_level(_rule_for(f, args...), Int(k))[1] === :pole
end

"""
    level_spectrum(j1, ..., j6; k = 2:100, cancellation = true) -> Vector{Symbol}
    level_spectrum(symbol, labels...; k, cancellation) -> Vector{Symbol}

How one symbol behaves as the level varies, one entry per level in `k`: `:inadmissible` (the labels are not
in the theory at that level, so the symbol is zero by convention), `:pole`, `:zero` (every term vanishes),
`:cancels` (finite terms summing to exactly zero) or `:finite`. Cheap because the valuations are
closed-form in the level; pass `cancellation = false` to skip the modular test and report `:finite`
wherever terms contribute.

```julia
level_spectrum(5, 5, 5, 5, 5, 5; k = 2:60)
```
"""
level_spectrum(args::Spin...; kw...) = level_spectrum(q6j, args...; kw...)

function level_spectrum(f::Function, args::Spin...; k = 2:100, cancellation::Bool = true, threads = nothing)
    K = k isa AbstractVector ? collect(k) : [k]
    out = Vector{Symbol}(undef, length(K))
    work = function (rng)
        for i in rng
            kk = Int(K[i])
            if !_admissible_at(f, kk, args...)
                out[i] = :inadmissible
                continue
            end
            s = _rule_for(f, args...)
            st, segs = classify_at_level(s, kk)
            out[i] = if st === :pole
                :pole
            elseif st === :empty || st === :zero
                :zero
            elseif cancellation && is_cancellation_zero(s, segs, kk) === true
                :cancels
            else
                :finite
            end
        end
    end
    nw = _nworkers(length(K), threads)
    if nw == 1
        work(eachindex(K))
    else
        @sync for c in _chunks(length(K), nw)
            Threads.@spawn work(c)
        end
    end
    return out
end
