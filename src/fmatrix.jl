# ---------------------------------------------------------------------------------
#  The F-matrix, [F^{abc}_d]_{ef} = fsymbol(a, b, e, c, d, f), orthogonal in this normalisation. With f fixed,
#  e runs over a 6j column, filled by the recurrence of `families.jl`. The row and column ranges and the
#  coefficients `Bn`, `Dn` are fixed for the whole matrix; only C = [(a−f−d)/2][(a+f−d)/2 + 1] changes per column.
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
    qq = _analytic_q(q)
    qq isa Real && return real_q_tables(qq, J2, J3, L1, L2, L3), nothing
    return complex_q_tables(ComplexF64(qq), J2, J3, L1, L2, L3), nothing
end

"""
Element type of the matrix: T is a precision floor, preserving complex values
and a higher-precision parameter, as in the scalar API.
"""
_fmatrix_bigq(q) = q !== nothing && real(typeof(float(q))) === BigFloat
function _fmatrix_eltype(T, q)
    E = T === nothing ? Float64 : T
    q === nothing || (E = promote_type(E,typeof(float(q))))
    return q isa Real && q < 0 ? Complex{real(E)} : E
end

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

with the intermediate spins labelling its rows (`e`) and columns (`f`). Entries are

    F[i, j] = fsymbol(a, b, e[i], c, d, f[j]) = (-1)^{a+b+c+d} √([2e+1][2f+1]) {a b e; c d f} ,

so `transpose(F) * F ≈ I`. At complex `q` the matrix is complex orthogonal, not unitary.

Use `k` for a level, `q` for a real or complex parameter, or neither for the classical limit (`k` and `q`
are exclusive). At a level the ranges hold only admissible intermediates, so `F` can be smaller or empty.
Negative real q uses the branch of `complex(q)`. `T` sets a precision floor: machine precision uses the
column recurrence (O(1) per entry); `BigFloat` requests and inputs use scalar `fsymbol` calls instead.

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
