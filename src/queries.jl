# --------------------------------------------
#  Structural queries: zeros, poles and level spectra
#
#  Valuations give poles and all-terms-vanishing zeros in closed form. Modular residues filter
#  cancellation candidates; definitive zero queries confirm unresolved candidates exactly.
# --------------------------------------------

# Both of these now dispatch through the symbol interface of `symbols.jl`, so a new symbol needs a
# `symbol_of` method and nothing here. Labels outside the level's admissible set carry no representation,
# so the symbol is zero by convention rather than singular, and the queries follow that convention.
_rule_for(f::Function, args...) = symbol_rule(_symbol_or_throw(f), args...)
_admissible_at(f::Function, k::Int, args...) = level_admissible(_symbol_or_throw(f), k, args...)

function _symbol_or_throw(f::Function)
    sym = symbol_of(f)
    sym === nothing && throw(ArgumentError(
        "$f is not a recoupling symbol; define `QRecoupling.symbol_of(::typeof($f))` to attach one"))
    return sym
end

"""
    iszero_at(k, j1, j2, j3, j4, j5, j6) -> Bool
    iszero_at(symbol, k, labels...) -> Bool
    iszero_at(k, labels::AbstractVector) -> BitVector

Decide whether a symbol vanishes at level `k`, including cancellation between finite terms. `symbol`
selects `q6j`, `fsymbol` or `gsymbol` (default `q6j`). Given a collection of
label tuples, the test runs over the batch with shared tables and threads.

Structural identities prove zeros directly. Other candidates pass through modular screening and exact
confirmation. In a batch, screening can use threads; exact confirmations run serially and can cost more
at high degree. Use `level_spectrum` for screening without exact confirmation.

```julia
iszero_at(20, 5, 5, 5, 5, 5, 5)          # true: a cancellation zero
iszero_at(fsymbol, 10, 1, 1, 1, 1, 1, 1)
count(iszero_at(6, all_6j(k = 6); threads=1))       # how many symbols vanish at level 6
```
"""
iszero_at(k::Integer, args::Spin...) = iszero_at(q6j, k, args...)

function iszero_at(f::Function, k::Integer, args::Spin...)
    _evaluation_q(k,nothing,false)
    k = Int(k)
    _admissible_at(f, k, args...) || return true
    return is_zero_at_level(_rule_for(f, args...), k)
end

iszero_at(k::Integer, labels::AbstractVector; threads = nothing) = iszero_at(q6j, k, labels; threads = threads)

function iszero_at(f::Function, k::Integer, labels::AbstractVector; threads = nothing)
    _evaluation_q(k,nothing,false)
    k = Int(k)
    n = length(labels)
    n == 0 && return falses(0)
    nw = _nworkers(n, threads)
    # Packed bits share storage words. Parallel workers must own independently writable bytes;
    # preserve the compact public return type by packing only after every worker has joined.
    out = nw == 1 ? falses(n) : Vector{Bool}(undef, n)
    pending = Vector{Bool}(undef, n)
    ztab = level_zero_table(k)                  # built once, read-only afterwards
    work = function (rng)
        for i in rng
            l = labels[i]
            result = !_admissible_at(f, k, l...) ? true : _level_zero_screen(_rule_for(f, l...), k, ztab)
            out[i] = result === true
            pending[i] = result === nothing
        end
    end
    if nw == 1
        work(1:n)
    else
        @sync for c in _chunks(n, nw)
            Threads.@spawn work(c)
        end
    end
    # Exact polynomial caches and field arithmetic follow the serial exact-batch policy.
    for i in eachindex(pending)
        pending[i] && (out[i] = iszero(exact_x(_rule_for(f, labels[i]...), k)))
    end
    return out isa BitVector ? out : BitVector(out)
end

"""
    issingular_at(k, j1, ..., j6) -> Bool
    issingular_at(symbol, k, labels...) -> Bool

Whether these labels are singular at level `k`: the level-k rule has a contributing term of negative
valuation, i.e. a q-factorial in a denominator that vanishes. Decided from the valuations alone, with no
arithmetic. `symbol` selects `q6j`, `fsymbol` or `gsymbol` (default `q6j`).

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
    _evaluation_q(k,nothing,false)
    k = Int(k)
    return classify_at_level(_rule_for(f, args...), Int(k))[1] === :pole
end

"""
    level_spectrum(j1, ..., j6; k = 2:100, cancellation = true, prove = false) -> Vector{Symbol}
    level_spectrum(symbol, labels...; k, cancellation, prove) -> Vector{Symbol}

How one symbol behaves as the level varies, one entry per level in `k`: `:inadmissible` (the labels are not
in the theory at that level, so the symbol is zero by convention), `:pole`, `:zero` (every term vanishes),
`:cancels` (the finite terms cancel) or `:finite`. Cheap because the valuations are closed-form in the
level; pass `cancellation = false` to skip the cancellation test and report `:finite` wherever terms
contribute.

By default `:cancels` is a modular *candidate*, not a proof. With `prove = true` every candidate is settled:
by the structural identities (pairwise or reflection cancellation) or, failing those, by exact arithmetic
at the level, so `:cancels` then means a proved zero and an unconfirmed candidate is reported as `:finite`.
Screening runs on threads; exact confirmations run serially afterwards, as in [`iszero_at`](@ref), and
can cost milliseconds each at high degree.

```julia
level_spectrum(5, 5, 5, 5, 5, 5; k = 2:60)                 # screen
level_spectrum(5, 5, 5, 5, 5, 5; k = 2:60, prove = true)   # every :cancels proved
```
"""
level_spectrum(args::Spin...; kw...) = level_spectrum(q6j, args...; kw...)

function level_spectrum(f::Function, args::Spin...; k = 2:100, cancellation::Bool = true,
                        prove::Bool = false, threads = nothing)
    K = k isa AbstractVector ? collect(k) : [k]
    _evaluation_q(K,nothing,false)
    out = Vector{Symbol}(undef, length(K))
    # entries left for exact confirmation: a candidate, or a screen that could not decide (`nothing`);
    # one byte per entry, so threaded workers never write neighbouring bits of a shared word
    pending = fill(false, length(K))
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
            elseif cancellation
                c = is_cancellation_zero(s, segs, kk)
                if prove && c !== false
                    if pairwise_zero(s) || reflection_zero(s, segs, kk)
                        :cancels                          # a structural proof, cheap enough for the worker
                    else
                        pending[i] = true; :finite        # settled exactly below
                    end
                else
                    c === true ? :cancels : :finite
                end
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
    # exact confirmations, serially (exact polynomial caches and field arithmetic follow the serial policy)
    for i in eachindex(K)
        pending[i] || continue
        iszero(exact_x(_rule_for(f, args...), Int(K[i]))) && (out[i] = :cancels)
    end
    return out
end
