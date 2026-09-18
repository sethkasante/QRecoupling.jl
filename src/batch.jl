# --------------------------------------------
#  Batched evaluation and level sweeps
#
#  A symbol function called with a collection of label tuples returns an array of values, building the
#  per-level tables once and spreading the work over threads. This is the amortisation that matters: the
#  tables are the expensive part and every symbol at a level shares them.
#
#  Worker tasks take the tables as arguments, so the common path touches no cache lock. Precision
#  escalation is task-local in Julia, so it runs on the workers too; exact and generic-q batches stay
#  serial because they go through the computer-algebra layer.
# --------------------------------------------

"Batches at least this long are threaded when more than one thread is available."
const BATCH_MIN_THREADED = 32

function _nworkers(n::Int, threads)
    threads === nothing || return clamp(Int(threads), 1, n)
    return (n >= BATCH_MIN_THREADED && Threads.nthreads() > 1) ? min(Threads.nthreads(), n) : 1
end

function _chunks(n::Int, nw::Int)
    per = cld(n, nw)
    return [((i - 1) * per + 1):min(i * per, n) for i in 1:cld(n, per)]
end

_as_vector(labels::AbstractVector) = labels
_as_vector(labels) = collect(labels)

"""
    _value_prefetched(s, k, tab, fallback) -> value

One level value from a prefetched table, with no cache lookups on the common path. Precision escalation is
task-local in Julia (>= 1.9, and the package requires >= 1.10), so this is safe to call from worker tasks.
"""
function _value_prefetched(s::FactorialSum, k::Int, tab::QIntTables{T}, fallback::F) where {T,F}
    v, st, segs = level_pass1(s, k, tab)
    st === :done && return v
    st === :fallback && return T(fallback())
    return level_escalate(s, segs, k, T, level_zero_table(k))
end

function _run(work::F, n::Int, threads) where {F}
    nw = _nworkers(n, threads)
    if nw == 1
        work(1:n)
    else
        @sync for c in _chunks(n, nw)
            Threads.@spawn work(c)
        end
    end
    return nothing
end

"""
    _level_batch(rule, fallback, labels, k, T, threads) -> Vector{T}

Values of one symbol family over many label tuples at a single level: the tables are built once and the
labels are spread over threads. `rule(label)` builds the factorial rule; `fallback(label)` is the generic
route for the rare label whose factorials fall outside the level tables.
"""
function _level_batch(rule::R, fallback::F, labels, k::Int, ::Type{T}, threads) where {R,F,T}
    L = _as_vector(labels)
    n = length(L)
    out = Vector{T}(undef, n)
    n == 0 && return out
    tab = qint_tables(T, k)                     # built once here, read-only afterwards
    _run(n, threads) do rng
        for i in rng
            l = L[i]
            out[i] = _value_prefetched(rule(l), k, tab, () -> fallback(l))
        end
    end
    return out
end

"Values for a collection of labels over a collection of levels: rows are labels, columns are levels."
function _level_grid(rule::R, fallback::F, labels, ks, ::Type{T}, threads) where {R,F,T}
    L = _as_vector(labels); K = _as_vector(ks)
    out = Matrix{T}(undef, length(L), length(K))
    for (j, k) in enumerate(K)
        kk = Int(k)
        out[:, j] = _level_batch(rule, l -> fallback(l, kk), L, kk, T, threads)
    end
    return out
end

"""
    _sweep_levels(symbol, args, ks; q, exact, T, threads) -> Vector

One symbol across a range of levels, which is what `k = 2:50` means in a symbol call. Each level needs its
own tables, so the levels are the unit of parallel work. Exact and generic-q sweeps stay serial, because
they go through the computer-algebra layer.
"""
function _sweep_levels(symbol::S, args, ks; q, exact::Bool, T::Type, threads) where {S}
    K = _as_vector(ks)
    n = length(K)
    if !isnothing(q) || exact
        return [symbol(args...; k = K[i], q = q, exact = exact, T = T) for i in 1:n]
    end
    out = Vector{T}(undef, n)
    n == 0 && return out
    _run(n, threads) do rng
        for i in rng
            out[i] = symbol(args...; k = Int(K[i]), T = T)
        end
    end
    return out
end

# ---- rules and fallbacks per symbol family ----

_pad3j(l) = length(l) == 5 ? (l[1], l[2], l[3], l[4], l[5], -l[4] - l[5]) : l

_rule_6j(l) = sixj_sum(canonical_spins(l...)...)
_rule_3j(l) = (p = _pad3j(l); threej_sum(doubled(p...)...))
_rule_f(l) = fsymbol_sum(doubled(l...)...)
_rule_g(l) = gsymbol_sum(canonical_spins(l...)...)

_fb_6j(l, k, ::Type{T}) where {T} = project_discrete(q6j_dcr(canonical_spins(l...)...), k, T)
_fb_3j(l, k, ::Type{T}) where {T} = (p = _pad3j(l); project_discrete(q3j_dcr(doubled(p...)...), k, T))
_fb_f(l, k, ::Type{T}) where {T} = project_discrete(fsymbol_dcr(doubled(l...)...), k, T)
_fb_g(l, k, ::Type{T}) where {T} = project_discrete(gsymbol_dcr(canonical_spins(l...)...), k, T)

"""
    _normalize_labels(labels, nlab, fname) -> Vector of label tuples

Accepts a collection of label tuples, or one label tuple of spins, and returns a vector of label tuples.
"""
function _normalize_labels(labels, nlab, fname)
    L = _as_vector(labels)
    isempty(L) && return L
    l = first(L)
    (l isa Union{Tuple,AbstractVector} && length(l) in nlab) && return L
    (l isa Spin && length(L) in nlab) && return [Tuple(L)]      # one label passed as a single tuple
    throw(ArgumentError("$fname(labels; ...) expects a collection of label tuples with " *
                        "$(join(nlab, " or ")) entries, or one such tuple; " *
                        "got an element of type $(typeof(l))"))
end

"""
    _batched(rule, fallback, labels, nlab, fname; k, q, exact, T, threads)

Shared implementation of the collection-argument methods: the fast level path when `k` is given and `q` is
not, otherwise an element-by-element fall back to the generic single-symbol route.
"""
function _batched(rule::R, fallback::F, symbol::S, labels, nlab, fname;
                  k, q, exact::Bool, T::Type, threads) where {R,F,S}
    labels = _normalize_labels(labels, nlab, fname)
    if !isnothing(k) && isnothing(q) && !exact
        k isa AbstractVector && return _level_grid(rule, (l, kk) -> fallback(l, kk, T), labels, k, T, threads)
        return _level_batch(rule, l -> fallback(l, Int(k), T), labels, Int(k), T, threads)
    end
    L = _as_vector(labels)
    if k isa AbstractVector
        return [symbol(l...; k = kk, q = q, exact = exact, T = T) for l in L, kk in k]
    end
    return [symbol(l...; k = k, q = q, exact = exact, T = T) for l in L]
end

for (f, rule, fb, nlab) in ((:q6j, :_rule_6j, :_fb_6j, :((6,))),
                            (:q3j, :_rule_3j, :_fb_3j, :((5, 6))),
                            (:fsymbol, :_rule_f, :_fb_f, :((6,))),
                            (:gsymbol, :_rule_g, :_fb_g, :((6,))))
    @eval function $f(labels::Union{AbstractVector,Base.Generator,Base.Iterators.Filter,Tuple{Any,Vararg{Any}}};
                      k = nothing, q = nothing, exact::Bool = false, T::Type = Float64, threads = nothing)
        _batched($rule, $fb, $f, labels, $nlab, $(string(f)); k = k, q = q, exact = exact, T = T, threads = threads)
    end
end

# ---- label generators ----

"""
    all_6j(; k, jmax = k, canonical = false)

Every 6j label set that is admissible at level `k`, as spin tuples (multiples of 1/2 as `Rational`), with
each spin at most `jmax`. With `canonical = true` only one representative per Regge symmetry class is kept,
which is what a sweep over a level wants.

```julia
labels = all_6j(k = 6)
vals = q6j(labels; k = 6)        # tables built once, work threaded
```
"""
function all_6j(; k::Integer, jmax::Real = k, canonical::Bool = false)
    k = Int(k)
    Jmax = min(2 * jmax, 2k) |> x -> Int(floor(x))
    out = NTuple{6,Rational{Int}}[]
    seen = canonical ? Set{NTuple{6,Int}}() : nothing
    half(J) = J // 2
    for J1 in 0:Jmax, J2 in 0:Jmax
        for J3 in abs(J1 - J2):2:min(J1 + J2, Jmax)
            J1 + J2 + J3 <= 2k || continue
            for J4 in 0:Jmax, J5 in 0:Jmax
                _qδ(J3, J4, J5, k) || continue
                for J6 in max(abs(J1 - J5), abs(J2 - J4)):2:min(J1 + J5, J2 + J4, Jmax)
                    (_qδ(J1, J5, J6, k) && _qδ(J2, J4, J6, k)) || continue
                    if canonical
                        c = canonical_spins(half(J1), half(J2), half(J3), half(J4), half(J5), half(J6))
                        c in seen && continue
                        push!(seen, c)
                    end
                    push!(out, (half(J1), half(J2), half(J3), half(J4), half(J5), half(J6)))
                end
            end
        end
    end
    return out
end
