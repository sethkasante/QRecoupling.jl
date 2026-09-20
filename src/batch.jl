# --------------------------------------------
#  Batched evaluation and level sweeps
#
#  A symbol function called with a collection of label tuples returns an array of values, building the
#  per-level tables once and spreading the work over threads. This is the amortisation that matters: the
#  tables are the expensive part and every symbol at a level shares them.
#
#  Worker tasks take the tables as arguments, so the common path touches no cache lock. Precision
#  scopes are task-local on Julia 1.12+. Older runtimes use serial batches because BigFloat precision
#  is shared there. Exact and generic-q batches also stay serial.
# --------------------------------------------

"Batches at least this long are threaded when more than one thread is available."
const BATCH_MIN_THREADED = 32

function _nworkers(n::Int, threads)
    VERSION < v"1.12" && return 1
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

One level value from a prefetched table. The caller validates level admissibility before
selecting a standard-symbol family shortcut and owns its workspace for the whole worker task.
"""
function _value_prefetched(s::FactorialSum, k::Int, tab::QIntTables{T}, fallback::F;
                           labels=nothing, family=nothing, workspace=nothing) where {T,F}
    v, st, segs = level_pass1(s,k,tab; family=family,workspace=workspace)
    st === :done && return v
    st === :fallback && return T(fallback())
    return level_escalate(s,segs,k,T,level_zero_table(k); labels=labels,workspace=workspace)
end

function _run(work::F, n::Int, threads) where {F}
    n == 0 && return nothing
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
function _level_batch(rule::R, fallback::F, labels, k::Int, ::Type{T}, threads;
                      family = nothing) where {R,F,T}
    L = _as_vector(labels)
    n = length(L)
    out = Vector{T}(undef, n)
    n == 0 && return out
    tab = qint_tables(T, k)                     # built once here, read-only afterwards
    J = T === Float64 ? _doubled_labels(L, family) : nothing
    if J === nothing
        _run(n, threads) do rng
            work = EvaluationWorkspace()
            for i in rng
                l = L[i]
                if !_rule_admissible(rule,l,k)
                    out[i] = zero(T)
                    continue
                end
                out[i] = _value_prefetched(rule(l),k,tab,()->fallback(l);
                                          family=_rule_family(rule),workspace=work)
            end
        end
    else
        rest = _family_pass!(out, J, family, LevelQ(tab, k), threads)
        _run_workspace(rest,n,threads) do i,work
            if !_qδtet(J[i]...,k)
                out[i] = zero(T)
            else
                out[i] = _value_prefetched(_family_rule(family,J[i]),k,tab,
                                           ()->_family_fallback(family,J[i],k,T);
                                           labels=family === Val(:sixj) ? J[i] : nothing,
                                           family=family,workspace=work)
            end
        end
    end
    return out
end

"Classical (q = 1) values over many labels: families by recurrence, the rest per symbol, threaded."
function _classical_batch(rule::R, labels, ::Type{T}, threads; family = nothing) where {R,T}
    L = _as_vector(labels)
    n = length(L)
    out = Vector{T}(undef, n)
    n == 0 && return out
    J = T === Float64 ? _doubled_labels(L, family) : nothing
    if J === nothing
        _run(n, threads) do rng
            work = EvaluationWorkspace()
            for i in rng
                out[i] = classical_value(rule(L[i]),T; workspace=work)
            end
        end
    else
        rest = _family_pass!(out, J, family, ClassicalQ(), threads)
        _run_workspace(rest,n,threads) do i,work
            out[i] = classical_value(_family_rule(family, J[i]), T;
                                     labels=family === Val(:sixj) ? J[i] : nothing,workspace=work)
        end
    end
    return out
end

"One reusable scratch per task, including when processing scattered leftover indices."
function _run_workspace(f::F, idx, n::Int, threads) where {F}
    _run(idx === nothing ? n : length(idx),threads) do rng
        work = EvaluationWorkspace()
        for j in rng
            f(idx === nothing ? j : idx[j],work)
        end
    end
    return nothing
end

"Run `f(i)` for i in `idx` (or 1:n), spread over threads."
function _run_over(f::F, idx, n::Int, threads) where {F}
    idx === nothing && return _run(rng -> foreach(f, rng), n, threads)
    _run(length(idx), threads) do rng
        for j in rng
            f(idx[j])
        end
    end
    return nothing
end

# ---- families: runs of labels along one spin go to the three-term recurrence ----
#
# A batch that contains whole columns {x j2 j3; l1 l2 l3} (x varying, the rest fixed) — a whole level from
# `all_6j`, an F-matrix built by a comprehension, a sweep over one spin — is computed column by column by the
# recurrence of `families.jl`: 40–85 ns per entry at any spin, every entry to ~1e-16, where a single symbol
# costs 0.1–1 µs and much more once its Racah sum cancels. Columns are found as runs of consecutive labels
# in which one spin steps by one (O(1) per label); any of the six positions is moved to the first by a
# tetrahedral symmetry. Entries within 10⁻¹² of the column's largest (near a node, or exact zeros) and
# labels outside runs take the certified single-symbol path, so exact zeros stay exact.

"Shortest run of labels worth a recurrence (below it the column setup costs more than it saves)."
const FAMILY_MIN_RUN = 4

"Symbol families that the recurrence serves: the value is the 6j times a per-entry factor."
_family(::typeof(q6j)) = Val(:sixj)
_family(::typeof(fsymbol)) = Val(:f)
_family(::typeof(gsymbol)) = Val(:g)
_family(_) = nothing

"The one position in which b steps up from a by one spin (2 in doubled labels), or 0."
@inline function _step_position(a::NTuple{6,Int}, b::NTuple{6,Int})
    p = 0
    for t in 1:6
        d = b[t] - a[t]
        d == 0 && continue
        (d == 2 && p == 0) || return 0
        p = t
    end
    return p
end

"Maximal runs (first, last, position) of consecutive labels stepping in one position."
function _find_runs(J::Vector{NTuple{6,Int}}, minrun::Int)
    runs = NTuple{3,Int}[]
    n = length(J)
    i = 1
    while i < n
        p = _step_position(J[i], J[i+1])
        if p == 0
            i += 1
            continue
        end
        j = i + 1
        while j < n && _step_position(J[j], J[j+1]) == p
            j += 1
        end
        if j - i + 1 >= minrun
            push!(runs, (i, j, p))
            i = j + 1
        else
            i = j
        end
    end
    return runs
end

"""
    DoubledLabels

Label tuples kept as doubled integers (J = 2j). It behaves as a vector of spin tuples, so every function
takes it like any other collection of labels, while the batch path reads the integers directly and skips the
conversion (~15 ns per label for `Rational` spins). `all_6j(; doubled = true)` returns one.
"""
struct DoubledLabels <: AbstractVector{NTuple{6,Rational{Int}}}
    J::Vector{NTuple{6,Int}}
end

Base.size(L::DoubledLabels) = size(L.J)
Base.@propagate_inbounds Base.getindex(L::DoubledLabels, i::Int) = L.J[i] .// 2
Base.IndexStyle(::Type{DoubledLabels}) = IndexLinear()

"The factorial rule of a family member from its doubled labels."
_family_rule(::Val{:sixj}, J::NTuple{6,Int}) = sixj_sum(J...)
_family_rule(::Val{:f}, J::NTuple{6,Int}) = fsymbol_sum(J...)
_family_rule(::Val{:g}, J::NTuple{6,Int}) = gsymbol_sum(J...)

"Doubled label of one spin tuple, without the rational arithmetic of the generic `doubled`."
@inline _doubled_spin(j::Integer) = 2 * Int(j)
@inline function _doubled_spin(j::Rational)
    d = denominator(j)
    d == 1 && return 2 * Int(numerator(j))
    d == 2 && return Int(numerator(j))
    throw(ArgumentError("spin labels must be multiples of 1/2, got $j"))
end
@inline function _doubled_spin(j::AbstractFloat)
    J = round(Int, 2j)
    J == 2j || throw(ArgumentError("spin labels must be multiples of 1/2, got $j"))
    return J
end

"Doubled labels of a family batch, or nothing (not a family symbol, too short, or not six labels)."
function _doubled_labels(L, family)
    (family === nothing || length(L) < FAMILY_MIN_RUN) && return nothing
    J = Vector{NTuple{6,Int}}(undef, length(L))
    for i in eachindex(J)
        l = L[i]
        length(l) == 6 || return nothing
        J[i] = (_doubled_spin(l[1]), _doubled_spin(l[2]), _doubled_spin(l[3]),
                _doubled_spin(l[4]), _doubled_spin(l[5]), _doubled_spin(l[6]))
    end
    return J
end

_doubled_labels(L::DoubledLabels, family) = family === nothing ? nothing : L.J

"The generic (DCR) route for one family member, from doubled labels."
_family_fallback(::Val{:sixj}, J::NTuple{6,Int}, k::Int, ::Type{T}) where {T} =
    project_discrete(q6j_dcr(J...), k, T)
_family_fallback(::Val{:f}, J::NTuple{6,Int}, k::Int, ::Type{T}) where {T} =
    project_discrete(fsymbol_dcr(J...), k, T)
_family_fallback(::Val{:g}, J::NTuple{6,Int}, k::Int, ::Type{T}) where {T} =
    project_discrete(gsymbol_dcr(J...), k, T)

"Factor from the 6j to the symbol, as a double word: 1, (−1)^{j1+j2+j4+j5}√([2j3+1][2j6+1]), or √Π[2j+1]."
@inline _family_factor(::Val{:sixj}, Q, J) = (1.0, 0.0)
@inline function _family_factor(::Val{:f}, Q, J)
    d = _dwsqrt(_dwm(_qint_dw(Q, J[3] + 1), _qint_dw(Q, J[6] + 1)))
    return iseven((J[1] + J[2] + J[4] + J[5]) ÷ 2) ? d : _dwn(d)
end
@inline function _family_factor(::Val{:g}, Q, J)
    p = (1.0, 0.0)
    for t in 1:6
        p = _dwm(p, _qint_dw(Q, J[t] + 1))
    end
    return _dwsqrt(p)
end

"""
    _family_pass!(out, J, family, Q, threads) -> indices still to compute, or nothing

Fills `out` for the labels (doubled, `J`) that sit in runs along one spin, from column recurrences; returns
the indices left for the single-symbol path (`nothing` means all of them).
"""
function _family_pass!(out::Vector{Float64}, J::Vector{NTuple{6,Int}}, family, Q, threads)
    n = length(J)
    runs = _find_runs(J, FAMILY_MIN_RUN)
    isempty(runs) && return nothing
    done = zeros(Bool, n)                        # one byte per entry: safe to write from several threads
    _run(length(runs), threads) do rr
        work = ColumnWork()
        for r in rr
            _family_run!(out, done, J, runs[r], family, Q, work)
        end
    end
    return [i for i in 1:n if !done[i]]
end

function _family_run!(out, done, J, run, family, Q, work)
    i0, i1, p = run
    σ = _TO_FIRST[p]
    c = J[i0]
    J2, J3, L1, L2, L3 = c[σ[2]], c[σ[3]], c[σ[4]], c[σ[5]], c[σ[6]]
    X = _column_range(Q, J2, J3, L1, L2, L3)
    nx = length(X)
    (nx >= FAMILY_MIN_RUN && nx <= 4 * (i1 - i0 + 1) + 8) || return nothing   # a small part of a long column
    sixj_column!(nothing, J2, J3, L1, L2, L3, Q, work)
    w = work.w
    mx = 0.0
    @inbounds for t in 1:nx
        mx = max(mx, abs(w[t][1]))
    end
    @inbounds for i in i0:i1
        d = J[i][p] - first(X)
        (iseven(d) && 0 <= d <= last(X) - first(X)) || continue
        v = w[d ÷ 2 + 1]
        abs(v[1]) >= 1e-12 * mx || continue      # near a node or an exact zero: the certified path decides
        vv = _dwm(v, _family_factor(family, Q, J[i]))
        out[i] = vv[1] + vv[2]
        done[i] = true
    end
    return nothing
end

"Values for a collection of labels over a collection of levels: rows are labels, columns are levels."
function _level_grid(rule::R, fallback::F, labels, ks, ::Type{T}, threads; family = nothing) where {R,F,T}
    L = _as_vector(labels); K = _as_vector(ks)
    out = Matrix{T}(undef, length(L), length(K))
    for (j, k) in enumerate(K)
        kk = Int(k)
        out[:, j] = _level_batch(rule, l -> fallback(l, kk), L, kk, T, threads; family = family)
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

_rule_6j(l) = sixj_sum(doubled(l...)...)          # symmetric by construction
_rule_3j(l) = (p = _pad3j(l); threej_sum(doubled(p...)...))
_rule_f(l) = fsymbol_sum(doubled(l...)...)
_rule_g(l) = gsymbol_sum(doubled(l...)...)        # symmetric by construction

_rule_family(::typeof(_rule_6j)) = Val(:sixj)
_rule_family(::typeof(_rule_3j)) = Val(:threej)
_rule_family(::typeof(_rule_f)) = Val(:f)
_rule_family(::typeof(_rule_g)) = Val(:g)
_rule_admissible(rule, l, k) = _qδtet(doubled(l...)...,k)
_rule_admissible(::typeof(_rule_3j),l,k) = _qδ(doubled(l[1],l[2],l[3])...,k)

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
    fam = _family(symbol)
    if !isnothing(k) && isnothing(q) && !exact
        k isa AbstractVector &&
            return _level_grid(rule, (l, kk) -> fallback(l, kk, T), labels, k, T, threads; family = fam)
        return _level_batch(rule, l -> fallback(l, Int(k), T), labels, Int(k), T, threads; family = fam)
    end
    _is_classical(q) && !exact && !(k isa AbstractVector) &&
        return _classical_batch(rule, labels, T, threads; family = fam)
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
    all_6j(; k, jmax = k, canonical = false, doubled = false)

Every 6j label set that is admissible at level `k`, as spin tuples (multiples of 1/2 as `Rational`), with
each spin at most `jmax`. With `canonical = true` only one representative per Regge symmetry class is kept,
which is what a sweep over a level wants. With `doubled = true` the labels come back as `DoubledLabels`,
which behaves the same but lets batches skip the conversion to doubled integers (~10% on a whole level).

```julia
labels = all_6j(k = 6)
vals = q6j(labels; k = 6)        # tables built once, work threaded
```
"""
function all_6j(; k::Integer, jmax::Real = k, canonical::Bool = false, doubled::Bool = false)
    k = Int(k)
    Jmax = min(2 * jmax, 2k) |> x -> Int(floor(x))
    out = NTuple{6,Int}[]
    seen = canonical ? Set{NTuple{6,Int}}() : nothing
    for J1 in 0:Jmax, J2 in 0:Jmax
        for J3 in abs(J1 - J2):2:min(J1 + J2, Jmax)
            J1 + J2 + J3 <= 2k || continue
            for J4 in 0:Jmax, J5 in 0:Jmax
                _qδ(J3, J4, J5, k) || continue
                for J6 in max(abs(J1 - J5), abs(J2 - J4)):2:min(J1 + J5, J2 + J4, Jmax)
                    (_qδ(J1, J5, J6, k) && _qδ(J2, J4, J6, k)) || continue
                    if canonical
                        c = canonical_spins(J1 // 2, J2 // 2, J3 // 2, J4 // 2, J5 // 2, J6 // 2)
                        c in seen && continue
                        push!(seen, c)
                    end
                    push!(out, (J1, J2, J3, J4, J5, J6))
                end
            end
        end
    end
    L = DoubledLabels(out)
    return doubled ? L : NTuple{6,Rational{Int}}[l for l in L]
end
