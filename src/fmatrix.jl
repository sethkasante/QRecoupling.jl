# ---------------------------------------------------------------------------------
#  The F-matrix (associativity / recoupling matrix)
#
#  For four spins a, b, c, d the two ways of fusing to d,
#
#      ((a b)_e c)_d     and     (a (b c)_f)_d ,
#
#  are related by a change of basis. Its matrix is the F-matrix, and its entries are exactly the package's
#  F-symbol:
#
#      [F^{abc}_d]_{ef} = (-1)^{a+b+c+d} √([2e+1][2f+1]) {a b e; c d f} = fsymbol(a, b, e, c, d, f).
#
#  In that normalisation F is **orthogonal** — Σ_e [F]_{ef}[F]_{ef'} = δ_{ff'} is the orthogonality relation
#  the package already verifies — so the inverse is the transpose. At a level the rows and columns run over
#  the level-admissible intermediates only, which is why both label vectors are returned with the matrix.
#
#  Computing it column by column is what makes it cheap. With f fixed, e runs over a *column* of 6j symbols,
#  which the three-term recurrence of `families.jl` fills in O(n) instead of n independent Racah sums — and
#  the row and column ranges do not depend on f at all, because e's triangle conditions are (e,a,b) and
#  (e,c,d). Better still, the recurrence coefficients are shared: `Bn`, `Dn` and hence `up(x−1)²` depend only
#  on (b, a, d, c), which are fixed for the whole matrix, and of `di` only the single scalar
#  C = [ (a-f-d)/2 ][ (a+f-d)/2 + 1 ] changes from column to column. So the four-factor products are formed
#  once for the entire matrix.
# ---------------------------------------------------------------------------------

"""
Evaluation context for a family, from the same keywords the symbol functions take: a level `k`, a real `q`,
or the classical limit.
"""
function _family_context(k, q, J2::Int, J3::Int, L1::Int, L2::Int, L3::Int)
    if k !== nothing
        q === nothing || throw(ArgumentError("give a level `k` or a parameter `q`, not both"))
        kk = Int(k)
        kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
        return LevelQ(qint_tables(Float64, kk), kk), kk
    end
    (q === nothing || _is_classical(q)) && return ClassicalQ(), nothing
    q isa Real && return real_q_tables(float(q), J2, J3, L1, L2, L3), nothing
    return complex_q_tables(ComplexF64(q), J2, J3, L1, L2, L3), nothing
end

"""
Element type of the matrix: whatever the caller asked for, else Float64/ComplexF64,
preserving BigFloat for a parameter supplied at that precision.
"""
_fmatrix_bigq(q) = q !== nothing && real(typeof(float(q))) === BigFloat
function _fmatrix_eltype(::Nothing, q)
    _fmatrix_bigq(q) && return q isa Real ? BigFloat : Complex{BigFloat}
    return (q !== nothing && !(q isa Real) && !_is_classical(q)) ? ComplexF64 : Float64
end
_fmatrix_eltype(T::Type, _) = T

"""
Validate the target keywords up front, so that an inadmissible or nonsensical level is reported even when
the label ranges come out empty and the evaluation context is never built.
"""
function _fmatrix_level(k, q)
    k === nothing && return nothing
    q === nothing || throw(ArgumentError("give a level `k` or a parameter `q`, not both"))
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
    return kk
end

"The intermediates `f` in b ⊗ c that also satisfy a ⊗ f ∋ d, in doubled labels."
function _fmatrix_cols(A::Int, B::Int, C::Int, D::Int, k)
    lo = max(abs(B - C), abs(A - D))
    hi = min(B + C, A + D)
    k === nothing && return lo:2:hi
    return lo:2:min(hi, 2k - B - C, 2k - A - D)
end

"""
    fmatrix(a, b, c, d; k=nothing, q=nothing, T=nothing) -> (F, e, f)

The F-matrix (associativity matrix) for four spins, relating the two fusion bases

    ((a b)ₑ c)_d   ⟷   (a (b c)_f)_d ,

together with the intermediate spins labelling its rows (`e`) and columns (`f`). Entries are

    F[i, j] = fsymbol(a, b, e[i], c, d, f[j]) = (-1)^{a+b+c+d} √([2e+1][2f+1]) {a b e; c d f} ,

which is the **orthogonal** normalisation: `F' * F ≈ I`, so the inverse change of basis is `transpose(F)`.

Use `k` for a unitary level, `q` for any deformation — positive real, or complex — or neither for the
classical limit; `k` and `q` are mutually exclusive, as elsewhere in the package. At a level the ranges are
the level-admissible intermediates, so `F` can be smaller than the classical matrix, or empty.

`T` defaults to the natural element type: `Float64` at a level, classically and at real `q`, `ComplexF64`
off the real axis; BigFloat inputs preserve their real or complex element type. Requests for `BigFloat`
or `Complex{BigFloat}` entries, and BigFloat inputs, use scalar `fsymbol` evaluations at the requested
precision instead of the machine-precision recurrence. This path is slower but retains the extra digits.
The orthogonality relation `Σ_e [F]_{ef} [F]_{ef'} = δ_{ff'}` is an algebraic identity,
so at complex `q` it makes `F` **complex orthogonal — `transpose(F) * F ≈ I`, not unitary**. `F' * F` is
not the identity there, and nothing says it should be: the sum is bilinear, not sesquilinear.

At machine precision the whole matrix is built from the column recurrence, at `O(1)` work per entry rather than one Racah sum
each, and the recurrence coefficients are shared across all columns. Off the real axis the column runs in
complex double words (`CDWord`), at the same `u²` per component; `q` at a root of unity is refused
there, because that is a level and has its own tables.

```julia
using LinearAlgebra
F, e, f = fmatrix(1, 1, 1, 1; k = 6)
F' * F ≈ I                      # orthogonal
F[1, 1] ≈ fsymbol(1, 1, e[1], 1, 1, f[1]; k=6)

G, e, f = fmatrix(1, 1, 1, 1; q = 0.8 + 0.3im)
transpose(G) * G ≈ I            # complex orthogonal off the real axis, not unitary
```
"""
function fmatrix(a::Spin, b::Spin, c::Spin, d::Spin;
                 k = nothing, q = nothing, T::Union{Type,Nothing} = nothing)
    A, B, C, D = doubled(a, b, c, d)
    kk = _fmatrix_level(k, q)
    fs = _fmatrix_cols(A, B, C, D, kk)
    E = _fmatrix_eltype(T, q)
    # e runs over the column of {x b a; f d c}; its range does not involve f
    es = isempty(fs) ? (0:2:-2) : sixj_column_range(B, A, first(fs), D, C, kk)
    F = zeros(E, length(es), length(fs))
    (isempty(es) || isempty(fs)) && return F, collect(es) .// 2, collect(fs) .// 2

    if real(E) === BigFloat || _fmatrix_bigq(q)
        for (j, f) in enumerate(fs), (i, e) in enumerate(es)
            F[i, j] = fsymbol(a, b, e // 2, c, d, f // 2; k=kk, q=q, T=BigFloat)
        end
        return F, collect(es) .// 2, collect(fs) .// 2
    end

    Q, _ = _family_context(k, q, B, A, first(fs), D, C)
    W = _wordtype(Q)
    work = ColumnWork(W)
    _resize!(work, length(es))
    sgn = iseven((A + B + C + D) ÷ 2) ? 1 : -1
    shared = Vector{W}(undef, length(es))              # the f-independent part of di
    _column_shared!(shared, work, Q, es, B, A, D, C)   # and work.e, for the whole matrix
    col = Vector{W === CDWord ? ComplexF64 : Float64}(undef, length(es))
    for (j, f) in enumerate(fs)
        _column_di!(work, Q, es, A, f, D, shared)
        # √[2f+1] in the same branch convention the column folds √[2e+1] in — the Ψ basis at complex q,
        # the ordinary positive root everywhere else. A principal root of the product here negated whole
        # columns against `fsymbol` at four of twelve sampled complex q.
        cf = _dws(_xdim_root(Q, f + 1), Float64(sgn))
        _sixj_column_core!(col, B, A, f, D, C, Q, work, es; xdim = true, cfac = cf)
        @inbounds for i in eachindex(es)
            F[i, j] = E(col[i])
        end
    end
    return F, collect(es) .// 2, collect(fs) .// 2
end

"""
    fmatrix_labels(a, b, c, d; k=nothing) -> (e, f)

The row and column intermediates of [`fmatrix`](@ref) without computing it — the admissible `e ∈ a⊗b` with
`e⊗c ∋ d`, and `f ∈ b⊗c` with `a⊗f ∋ d`.
"""
function fmatrix_labels(a::Spin, b::Spin, c::Spin, d::Spin; k = nothing)
    A, B, C, D = doubled(a, b, c, d)
    kk = _fmatrix_level(k, nothing)
    fs = _fmatrix_cols(A, B, C, D, kk)
    es = isempty(fs) ? (0:2:-2) : sixj_column_range(B, A, first(fs), D, C, kk)
    return collect(es) .// 2, collect(fs) .// 2
end
