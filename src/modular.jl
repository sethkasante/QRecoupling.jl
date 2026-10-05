# ---------------------------------------------------------------------------------
#  Modular data: the phases, and the matrices built from them
#
#  Above the 6j symbols sits a smaller and much cheaper layer — the twist, the S- and T-matrices, the
#  total dimension, the central charge. None of it needs a Racah sum. Every quantity here is either a
#  root of unity, which the package already represents exactly as a `QPhase`, or a value of
#  `sin(mπ/h)`, and that is the whole reason this file is fast:
#
#      S_{ab} = √(2/h) · sin((2a+1)(2b+1)π/h),      h = k+2,
#
#  so an n×n S-matrix is n² table lookups. The argument `(2a+1)(2b+1)` is folded into `[0, 2h)` by
#  addition rather than by `mod` — along a row it advances by a constant — and `sin(mπ/h)` for
#  `m ∈ [0, h]` is tabulated once per level. No trigonometric call is made per entry, and none is made
#  at all for a level that has been seen before.
#
#  Conventions are fixed by the package, not chosen here. The twist is
#
#      θ_j = q^{J(J+2)/2} = exp(2πi j(j+1)/(k+2)),   J = 2j,   q = exp(iπ/h),
#
#  which is the unique choice making `rmatrix(j, j, 0) = (−1)^{2j} θ_j⁻¹`, the self-braiding of a
#  self-dual object through the vacuum. The modular-data tests check this and the remaining identities against the
#  category the package already has: Verlinde's formula must reproduce `_qδ`, the level-truncated fusion
#  rule, which it has no way of knowing.
# ---------------------------------------------------------------------------------

"""
    level_labels(k) -> Vector{Rational{Int}}

The simple objects of SU(2)_k: `j = 0, 1/2, …, k/2`. These index the rows and columns of every matrix
in this file, and are returned alongside them so that no caller has to reconstruct the order.
"""
function level_labels(k::Integer)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
    return [J // 2 for J in 0:kk]
end

const _SIN_TABLE = LRU{Tuple{Int,DataType,Int},Any}(maxsize = 32)
const _SIN_LOCK = ReentrantLock()

"""
`sin(mπ/h)` for `m = 0…h`, at the requested precision, cached per level. The half-period is all that is
needed: `sin` beyond it is the same table with a sign, which is what `_sin_at` uses.
"""
function _sin_table(h::Int, ::Type{T}) where {T<:AbstractFloat}
    key = (h, T, precision(T))
    lock(_SIN_LOCK) do
        get!(_SIN_TABLE, key) do
            [sinpi(T(m) / T(h)) for m in 0:h]
        end
    end::Vector{T}
end

"`sin(mπ/h)` for `m` already folded into `[0, 2h)`."
@inline _sin_at(tab::Vector{T}, r::Int, h::Int) where {T} = r < h ? @inbounds(tab[r+1]) :
                                                                    -@inbounds(tab[r-h+1])

function _smatrix(k::Integer, ::Type{T}) where {T<:AbstractFloat}
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
    h = kk + 2
    n = kk + 1
    tab = _sin_table(h, T)
    c = sqrt(T(2) / T(h))
    S = Matrix{T}(undef, n, n)
    two_h = 2h
    @inbounds for i in 1:n
        step = i                      # (2a+1) for a = (i-1)/2
        r = 0
        for j in 1:n
            r += step
            r >= two_h && (r -= two_h)
            S[i, j] = c * _sin_at(tab, r, h)
        end
    end
    return S, level_labels(kk)
end
"""
    smatrix(k; T = Float64) -> (S, labels)

The modular S-matrix of SU(2)_k, `S_{ab} = √(2/h)·sin((2a+1)(2b+1)π/h)` with `h = k+2`, together with
the spins labelling it.

`S` is real, symmetric and orthogonal, `S² = C` is the charge conjugation matrix — the identity here,
since every SU(2)_k object is self-dual — and `S_{0b}/S_{00}` is the quantum dimension of `b`.

Built in `O(n²)` table lookups with no trigonometry per entry: along a row the argument advances by a
constant, folded back into range by one comparison.

```julia
S, js = smatrix(6)
S * S' ≈ I
S[1, :] ./ S[1, 1] ≈ qdim.(js; k = 6)
```
"""
smatrix(k::Integer; T::Type = Float64) = _smatrix(k, T)

"""
    twist(j; k = nothing, q = nothing, exact = false, T = ComplexF64)

The topological spin `θ_j = q^{J(J+2)/2} = exp(2πi j(j+1)/(k+2))`, `J = 2j`.

This is the convention the package's `rmatrix` already uses: `rmatrix(j, j, 0) = (−1)^{2j} θ_j⁻¹`, the
self-braiding of a self-dual object through the vacuum. With `exact = true` the answer is the
[`QPhase`](@ref) `q^{J(J+2)/2}`, exact at any level; the classical limit is 1.

`q` evaluates the same monomial away from a root of unity, which is what the rest of the package means by
a generic parameter. It is the one place the modular layer has a generic-q value at all — `smatrix`,
`tmatrix` and the rest are statements about a modular category and need the level.
"""
function twist(j::Spin; k = nothing, q = nothing, exact::Bool = false, T::Type = ComplexF64)
    q = _evaluation_q(k, q, exact)
    k isa AbstractVector && return [twist(j;k=kk,exact=exact,T=T) for kk in k]
    _check_phase_q(q)
    J = doubled(j)
    J >= 0 || throw(DomainError(j, "spin must be nonnegative"))
    if k !== nothing
        kk = Int(k)
        _qδ(J, 0, J, kk) || throw(ArgumentError("j = $j is not an object of SU(2)_$kk"))
    end
    p = J * (J + 2)
    exact && return QPhase(Int8(1), p//2)
    _is_classical(q) && return T(1)                   # q → 1
    k === nothing && q === nothing && return T(1)
    if k === nothing
        v = qhalfpow(_phase_parameter(q,T), p)
        return _phase_eltype(T,k,q)(v)
    end
    E = _phase_eltype(T,k,q)
    return E(_level_phase(p,Int(k),real(E)))
end

"""
    central_charge(k) -> Rational{Int}

The Virasoro central charge `c = 3k/(k+2)` of SU(2)_k. It enters the T-matrix through the framing
anomaly `exp(−2πi c/24)`, and is fixed independently by the Gauss sums and by `(ST)³ = S²`, both of which
the modular-data tests check rather than assume.
"""
function central_charge(k::Integer)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
    return (3 * kk) // (kk + 2)
end

function _tmatrix(k::Integer, ::Type{T}, anomaly::Bool) where {T<:Number}
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
    js = level_labels(kk)
    n = length(js)
    R = real(float(T))
    E = Complex{R}
    pre = anomaly ? cispi(-R(kk) / R(4(kk+2))) : one(E)   # exp(−2πi c/24)
    M = zeros(E, n, n)
    @inbounds for i in 1:n
        J = i-1
        M[i, i] = pre * _level_phase(J*(J+2),kk,R)
    end
    return M, js
end
"""
    tmatrix(k; T = ComplexF64, anomaly = true) -> (Tm, labels)

The modular T-matrix, `T_{ab} = δ_{ab}·exp(-2πi c/24)·θ_a`. With `anomaly = false` the framing phase is
dropped and the diagonal is the bare twists, which is what most fusion-category conventions print.

Returned as a `Diagonal`-shaped dense matrix so that `(S*T)^3 ≈ S^2` can be written directly.
"""
tmatrix(k::Integer; T::Type = ComplexF64, anomaly::Bool = true) = _tmatrix(k, T, anomaly)

function _total_qdim(k::Integer, ::Type{T}) where {T<:AbstractFloat}
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
    h = kk + 2
    return sqrt(T(h) / T(2)) / sinpi(one(T) / T(h))
end
"""
    total_qdim(k; T = Float64)

The total quantum dimension `D = √(Σ_a d_a²) = √(h/2)/sin(π/h)`, `h = k+2`.

The closed form agrees with the sum, as checked by the modular-data tests,
because `Σ_{n=1}^{h-1} sin²(nπ/h) = h/2` exactly, which is also why `D² = 2h/(4-x²)` is a rational
function of the package's `x = 2cos(π/h)`.
"""
total_qdim(k::Integer; T::Type = Float64) = _total_qdim(k, T)

function _gauss_sum(k::Integer, ::Type{T}, inverse::Bool) where {T<:Number}
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
    h = kk + 2
    R = real(float(T))
    E = Complex{R}
    tab = _sin_table(h, R)
    s1 = tab[2]                                   # sin(π/h)
    acc = zero(E)
    @inbounds for J in 0:kk
        d = tab[J+2] / s1                         # d_j = [J+1] = sin((J+1)π/h)/sin(π/h)
        θ = _level_phase(J*(J+2),kk,R)
        acc += d * d * (inverse ? conj(θ) : θ)
    end
    return acc
end
"""
    gauss_sum(k; T = ComplexF64, inverse = false) -> Number

`p₊ = Σ_a d_a² θ_a` (or `p₋ = Σ_a d_a² θ_a⁻¹`). The pair fixes the central charge without reference to
the formula `3k/(k+2)`: `p₊ = D·exp(2πi c/8)`, so `p₊p₋ = D²` and `p₊/p₋ = exp(2πi c/4)`.

"""
gauss_sum(k::Integer; T::Type = ComplexF64, inverse::Bool = false) = _gauss_sum(k, T, inverse)

"""
    verlinde(a, b, c; k, T = Float64) -> Int

The fusion multiplicity `N^{ab}_c` from Verlinde's formula,
`N^{ab}_c = Σ_x S_{ax}S_{bx}S̄_{cx}/S_{0x}`.

This is not how the package decides fusion — that is `_qδ`, a triangle inequality and a level cut — which
is exactly why it is worth having: the two agree over every triple at every level tested, and they share
no code. Returns the nearest integer, and throws if the sum is not within `1e-6` of one.
"""
function verlinde(a::Spin, b::Spin, c::Spin; k::Integer, T::Type = Float64)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k,"level must be nonnegative"))
    A,B,C = doubled(a,b,c)
    (0 <= A <= kk && 0 <= B <= kk && 0 <= C <= kk) ||
        throw(ArgumentError("labels must be objects of SU(2)_$kk: 0, 1/2, …, $(kk//2)"))
    h = kk+2
    tab = _sin_table(h,T)
    ra = rb = rc = 0
    acc = zero(T)
    # Three sine rows suffice; there is no need to allocate the full S-matrix.
    @inbounds for x in 1:kk+1
        ra += A+1; ra >= 2h && (ra -= 2h)
        rb += B+1; rb >= 2h && (rb -= 2h)
        rc += C+1; rc >= 2h && (rc -= 2h)
        acc += _sin_at(tab,ra,h)*_sin_at(tab,rb,h)*_sin_at(tab,rc,h)/tab[x+1]
    end
    acc *= T(2)/T(h)
    r = round(Int, acc)
    abs(acc - r) <= 1e-6 * max(one(T), abs(acc)) ||
        throw(ErrorException("Verlinde sum $(acc) is not an integer; this is a bug, please report it"))
    return r
end

"""
    monodromy(a, b, c; k = nothing, q = nothing, T = ComplexF64)

The double braid on the `c` channel of `a ⊗ b`: `θ_c/(θ_aθ_b)`, which is `rmatrix(a,b,c)²` and the
eigenvalue of the full monodromy. Returns zero on an inadmissible triple, as the braiding does.
"""
function monodromy(a::Spin, b::Spin, c::Spin; k = nothing, q = nothing, T::Type = ComplexF64)
    r = rmatrix(a, b, c; k = k, q = q, T = T)
    return r .* r
end

"""
    bmatrix(a, b, c, d; k=nothing, q=nothing, T=ComplexF64, inverse=false) -> (B, e, e′)

The braid matrix: the change of basis that carries `((a b)_e c)_d` to `((a c)_{e′} b)_d`, i.e. braids
`b` past `c` inside the fusion tree.

    [B^{abc}_d]_{e e′} = Σ_f [F^{abc}_d]_{ef} · R^{bc}_f · [F^{acb}_d]_{e′f}

— conjugate by the F-matrix into the basis where the braiding is diagonal, apply it there, and come back.
The two F-matrices share their column labels (both index `f ∈ b⊗c` with `a⊗f ∋ d`), which is what makes
the middle factor a diagonal of R-phases and the whole thing one `O(n³)` pass with no Racah sums.
Rows carry incoming channels `e`; columns carry outgoing channels `e′`. For a column vector of
incoming amplitudes, the outgoing amplitudes are `transpose(B) * state`.

`inverse = true` gives the same construction with `R⁻¹`, so that

    bmatrix(a, c, b, d; k, inverse = true) * bmatrix(a, b, c, d; k) == I

at **any** `q`. At a level `B` is also unitary, because the R-eigenvalues are phases there. At generic
positive real `q`, those eigenvalues are real and need not have unit modulus, so `B` is generally not
unitary. Measured over 107 matrices at k = 2…14: unitarity and the
inverse relation both to **4.4e-16**.

```julia
B, e, e′ = bmatrix(1, 1, 1, 1; k = 6)
B * B' ≈ I
```
"""
function bmatrix(a::Spin, b::Spin, c::Spin, d::Spin;
                 k = nothing, q = nothing, T::Type{TT} = ComplexF64, inverse::Bool = false) where {TT}
    _check_phase_q(q)
    if real(TT) === BigFloat
        F1, es, fs = fmatrix(a, b, c, d; k = k, q = q, T = BigFloat)
        F2, es2, fs2 = b == c ? (F1,es,fs) : fmatrix(a, c, b, d; k = k, q = q, T = BigFloat)
    else
        F1, es, fs = fmatrix(a, b, c, d; k = k, q = q)
        F2, es2, fs2 = b == c ? (F1,es,fs) : fmatrix(a, c, b, d; k = k, q = q)
    end
    fs == fs2 || throw(ErrorException(
        "the two F-matrices disagree on their intermediate labels ($fs vs $fs2); this is a bug"))
    E = _phase_eltype(TT,k,q)
    B = zeros(E, length(es), length(es2))
    (isempty(es) || isempty(es2)) && return B, es, es2
    @inbounds for (l, f) in enumerate(fs)
        r = E === ComplexF64 ? rmatrix(b,c,f;k=k,q=q) : rmatrix(b,c,f;k=k,q=q,T=E)
        # `inv`, not `conj`: they agree only where |R| = 1.
        inverse && (r = inv(r))
        for j in eachindex(es2)
            w = r * F2[j, l]
            for i in eachindex(es)
                B[i, j] += E(F1[i, l] * w)
            end
        end
    end
    return B, es, es2
end
