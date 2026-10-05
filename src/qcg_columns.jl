# ---------------------------------------------------------------------------------
#  Coupling columns by the Casimir recurrence: whole coupling matrices, and a near-edge tier for single
#  coefficients
#
#  Δ(Ω), Ω = FE + [m+½]², acts on |j m⟩ by [j+½]². At fixed m it is symmetric tridiagonal in the product basis,
#  so the column C(m1) = ⟨j1 m1; j2 m−m1|j m⟩, m1 = lo … hi, solves
#
#    A(m1) C(m1+1) + (D(m1) − [j+½]²) C(m1) + A(m1−1) C(m1−1) = 0,
#    A(m1) = q^{m2−m1−1} √([j1+m1+1][j1−m1][j2−m2+1][j2+m2]),                     m2 = m − m1
#    D(m1) = [j1−m1][j1+m1+1] q^{2m2} + q^{−2m1}[j2−m2][j2+m2+1] + [m+½]²
#
#  (from Δ(F)Δ(E) = FE⊗K² + FK⁻¹⊗KE + K⁻¹E⊗FK + K⁻²⊗FE and [x−y][x+y] = [x]² − [y]²). In r = q^{1/2},
#  [x]_q = [2x]_r/[2]_r, and the common [2]_r² drops out, so every coefficient is a product of integer
#  r-numbers at doubled labels and a power of r — no half-integer q-numbers and no fractional powers.
#
#  Run as in `families.jl`, two-sided in double words. Scale: at real q, Σ C² = 1 (positive terms) with
#  C(hi) > 0; at a level, one certified single-term endpoint. Error: relative n·u² on the forward stretch,
#  absolute beyond the meeting point. Measured ≤ 4e−16 against 256-bit references.
# ---------------------------------------------------------------------------------

"""
Tables for coupling columns: `[n]_r` for n = 0…I and `r^n` for n = −P…P, r = q^{1/2}, in double words.
`k` is the level, or −1 off a level.
"""
struct CGTables{P}
    ints::Vector{DWNum}
    pows::Vector{P}
    I::Int
    P::Int
    k::Int
    exact::Bool          # q = 1: every coefficient is an exact integer combination, so cancellation rounds nothing
end
CGTables{V}(ints, pows, I, P, k) where {V} = CGTables{V}(ints, pows, I, P, k, false)

@inline _cgint(T::CGTables, n::Int) = @inbounds T.ints[n+1]
@inline _cgpow(T::CGTables, n::Int) = @inbounds T.pows[n+T.P+1]

_dw_of(x::BigFloat) = DWNum(x)

"Real q > 0: built in 128 bits from r = √q itself, never from a rounded r (a rounded root would be amplified
by the powers, as a rounded e^{iπ/h} was at a level)."
function _cg_tables_real(q::Float64, I::Int, P::Int)
    setprecision(BigFloat, 128) do
        r = sqrt(BigFloat(q))
        x = r + inv(r)
        ints = Vector{DWNum}(undef, I + 1)
        prev = zero(BigFloat); cur = one(BigFloat)
        ints[1] = zero(DWNum)
        for n in 1:I                                   # [n+1] = x[n] − [n−1], the dominant solution
            ints[n+1] = DWNum(cur)
            prev, cur = cur, x * cur - prev
        end
        pows = Vector{DWNum}(undef, 2P + 1)
        p = one(BigFloat); ri = inv(r); pm = one(BigFloat)
        pows[P+1] = one(DWNum)
        for n in 1:P
            p *= r; pm *= ri
            pows[P+1+n] = DWNum(p); pows[P+1-n] = DWNum(pm)
        end
        return CGTables{DWNum}(ints, pows, I, P, -1)
    end
end

_cg_tables_classical(I::Int, P::Int) =
    CGTables{DWNum}([DWNum(Float64(n)) for n in 0:I], fill(one(DWNum), 2P + 1), I, P, -1, true)

"""
A level: r = e^{iπ/(2h)} is the root of level 2k + 2, whose tables already hold [n]_r as double words. The
powers follow from them without transcendental functions: r^n + r^{−n} = [n+1]_r − [n−1]_r and
r^n − r^{−n} = [n]_r (r − r^{−1}), with sin(π/2h) the one number taken from 128 bits.
"""
function _cg_tables_level(k::Int, I::Int, P::Int)
    Q = LevelQ(qint_tables(Float64, 2k + 2), 2k + 2)
    d(n) = (t = _qint_dw(Q, n); DWNum(t[1], t[2]))
    ints = [d(n) for n in 0:I]
    s = setprecision(() -> DWNum(sinpi(inv(BigFloat(2(k + 2))))), BigFloat, 128)
    half = DWNum(0.5)
    pows = [Complex(half * (d(n + 1) - d(n - 1)), d(n) * s) for n in -P:P]
    return CGTables{Complex{DWNum}}(ints, pows, I, P, k)
end

# Real-q tables are cached on q and grown on demand; level tables live in a `LevelCache`.
const CG_REAL_TABLES = LRU{Float64,Any}(maxsize = 16)
const CG_REAL_LOCK = ReentrantLock()
const CG_LEVEL_TABLES = LevelCache{CGTables{Complex{DWNum}}}()
const CG_CLASSICAL_TABLES = Ref{Any}(nothing)

"Table sizes a block (j1, j2) can reach, doubled labels."
_cg_sizes(J1::Int, J2::Int) = (2(J1 + J2) + 4, 2(J1 + J2) + 4)

function _cg_tables(q, k, J1::Int, J2::Int)
    I, P = _cg_sizes(J1, J2)
    if k !== nothing
        kk = Int(k)
        # one table per level, at the largest size the level admits (labels are bounded by k)
        Ik, Pk = _cg_sizes(kk, kk)
        return get_level!(() -> _cg_tables_level(kk, Ik, Pk), CG_LEVEL_TABLES, kk)
    elseif q === nothing || _is_classical(q)
        t = CG_CLASSICAL_TABLES[]
        (t isa CGTables{DWNum} && t.I >= I && t.P >= P) && return t
        t = _cg_tables_classical(max(I, 64), max(P, 64)); CG_CLASSICAL_TABLES[] = t
        return t
    end
    qq = Float64(q)
    t = @lock CG_REAL_LOCK get(CG_REAL_TABLES, qq, nothing)
    (t isa CGTables{DWNum} && t.I >= I && t.P >= P) && return t
    I2 = t isa CGTables{DWNum} ? max(I, t.I) : I; P2 = t isa CGTables{DWNum} ? max(P, t.P) : P
    t = _cg_tables_real(qq, max(I2, 32), max(P2, 32))
    @lock CG_REAL_LOCK (CG_REAL_TABLES[qq] = t)
    return t
end

clear_cg_caches!() = (@lock CG_REAL_LOCK empty!(CG_REAL_TABLES); empty!(CG_LEVEL_TABLES);
                      CG_CLASSICAL_TABLES[] = nothing; nothing)

"""
Scratch for the columns of one sector. `a`, `ia` and `d` are the sector's coefficients — A(m1), 1/A(m1) and
D(m1) without λ — which do not depend on j: every column of a sector shares them and differs only in
λ = [j+½]². So the square roots, divisions and table reads are paid once per sector, not once per column.
"""
struct CGWork{V,O}
    a::Vector{V}
    ia::Vector{V}
    d::Vector{V}
    dmag::Vector{Float64}     # Σ|terms| of d[k], to measure how much the diagonal cancels
    w::Vector{V}
    out::Vector{O}
    est::Vector{Float64}
end
function CGWork(::CGTables{P}, n::Int) where {P}
    O = P === DWNum ? Float64 : ComplexF64
    return CGWork{P,O}(Vector{P}(undef, n), Vector{P}(undef, n), Vector{P}(undef, n), Vector{Float64}(undef, n),
                       Vector{P}(undef, n), Vector{O}(undef, n), Vector{Float64}(undef, n))
end

# magnitudes for growth tests and rescaling, from the high words: no square root in the loop
@inline _mag(x::DWNum) = abs(x.hi)
@inline _mag(z::Complex{DWNum}) = abs(real(z).hi) + abs(imag(z).hi)

const CG_EST = 64.0     # the n·u² constant of the error estimate; measured errors are ≤ 1/20 of it

"""
    _cg_sector!(work, T, J1, J2, M) -> (lo, n)

The coefficients of the sector m (doubled labels): m1 = lo, lo+1, … over n values.
"""
function _cg_sector!(work::CGWork{V}, T::CGTables, J1::Int, J2::Int, M::Int) where {V}
    lo = max(-J1, M - J2); hi = min(J1, M + J2)
    n = (hi - lo) ÷ 2 + 1
    mh = _cgint(T, abs(M + 1))^2
    # in r = q^{1/2}: [x]_q ∝ [2x]_r and q^x = r^{2x}; the common [2]_r² is dropped
    @inbounds for k in 1:n
        M1 = lo + 2(k - 1); M2 = M - M1
        t1 = _cgint(T, J1 - M1) * _cgint(T, J1 + M1 + 2) * _cgpow(T, 2M2)
        t2 = _cgpow(T, -2M1) * _cgint(T, J2 - M2) * _cgint(T, J2 + M2 + 2)
        work.d[k] = t1 + t2 + mh
        work.dmag[k] = _mag(t1) + _mag(t2) + _mag(mh)
        if k < n
            a = _cgpow(T, M2 - M1 - 2) *
                sqrt(_cgint(T, J1 + M1 + 2) * _cgint(T, J1 - M1) * _cgint(T, J2 - M2 + 2) * _cgint(T, J2 + M2))
            work.a[k] = a; work.ia[k] = inv(a)
        end
    end
    return lo, n
end

"""
    _cg_twosided!(w, a, ia, d, λ, n) -> meet

The symmetric three-term recurrence a[k] w[k+1] + (d[k] − λ) w[k] + a[k−1] w[k−1] = 0 with w[0] = w[n+1] = 0,
up to scale: forward from 1 to the first local maximum, backward from n to it (each in the direction in which
the solution grows), the left part scaled to match. Returns the meeting index. Both the Casimir columns (in m1)
and the dual rows (in j) are this recurrence.
"""
function _cg_twosided!(w::Vector{V}, a::Vector{V}, ia::Vector{V}, d::Vector{V}, λ, n::Int) where {V}
    big = 2.0^500
    meet = 1
    w[1] = one(V)
    n == 1 && return meet
    # forward from 1, to the first local maximum
    w[2] = -((d[1] - λ) * w[1]) * ia[1]
    meet = n
    if _mag(w[2]) < _mag(w[1])
        meet = 1
    else
        @inbounds for k in 2:n-1
            w[k+1] = -((d[k] - λ) * w[k] + a[k-1] * w[k-1]) * ia[k]
            if _mag(w[k+1]) > big
                for i in 1:k+1; w[i] = _aldexp(w[i], -500); end
            end
            _mag(w[k+1]) < _mag(w[k]) && (meet = k; break)
        end
    end
    # backward from n down to the meeting point, then the left part scaled to it
    if meet < n
        fm = w[meet]
        w[n] = one(V)
        w[n-1] = -((d[n] - λ) * w[n]) * ia[n-1]
        @inbounds for k in n-1:-1:meet+1
            w[k-1] = -((d[k] - λ) * w[k] + a[k] * w[k+1]) * ia[k-1]
            if _mag(w[k-1]) > big
                for i in k-1:n; w[i] = _aldexp(w[i], -500); end
            end
        end
        sc = w[meet] / fm
        @inbounds for i in 1:meet-1; w[i] *= sc; end
    end
    return meet
end

"Scale `w` to `out` by a certified value `(m, e)` at index `i`; for a level, where the entries are complex."
function _cg_scale_to!(out, w::Vector{Complex{DWNum}}, n::Int, i::Int, m::ComplexF64, e::Int)
    sc = Complex{DWNum}(DWNum(real(m)), DWNum(imag(m))) / w[i]
    @inbounds for k in 1:n
        z = w[k] * sc
        out[k] = _aldexp(ComplexF64(Float64(real(z)), Float64(imag(z))), e)
    end
    return out
end

"Scale `w` to unit norm with the entry at `i` positive (real q > 0 and q = 1)."
function _cg_normalise!(out, w::Vector{DWNum}, n::Int, i::Int)
    ss = zero(DWNum)
    @inbounds for k in 1:n; ss += w[k] * w[k]; end
    sc = inv(sqrt(ss))
    signbit(w[i].hi) && (sc = -sc)
    @inbounds for k in 1:n; out[k] = Float64(w[k] * sc); end
    return out
end

"""
Per-entry error estimates: relative n·u² on the forward stretch (up to `meet`), absolute against the backward
branch's running maximum beyond it, plus the final rounding and the relative error of the scale; the u² term is
multiplied by the cancellation of the recurrence's diagonal.
"""
function _cg_estimates!(est, out, n::Int, meet::Int, seedrel::Float64, cancel::Float64 = 1.0)
    u2 = eps(Float64)^2
    run = 0.0
    @inbounds for k in n:-1:1
        x = abs(out[k])
        run = max(run, x)
        # `cancel`: how far the recurrence's diagonal cancelled when it was formed — its rounding, u²·cancel,
        # perturbs every step; at strong deformation it decides whether the recurrence can be used at all
        e = CG_EST * n * u2 * (1 + cancel) * (k <= meet ? x : run) / x + 2eps(Float64) + seedrel
        est[k] = isfinite(e) ? e : Inf
    end
    return est
end

"The certified value of a single-term coefficient at a level, from the level pass: (mantissa, exponent, error)."
function _cg_level_seed(k::Int, J1::Int, M1::Int, J2::Int, M2::Int, J::Int)
    s, (wt, e2) = qcg_rule(J1, M1, J2, M2, J, M1 + M2)
    v, relb, _ = _weighted_level_pass(s, wt, e2, k, Float64)
    return v.m, v.e, Float64(relb)
end

"""
    _cg_column!(work, T, J1, J2, J, M, lo, n)

The column C(j1 m1, j2 m−m1 | j m), m1 = lo, lo+1, …, in `work.out[1:n]`, with relative error estimates in
`work.est[1:n]`; `_cg_sector!` must have filled the sector's coefficients. Labels admissible, and inside the
fusion rule at a level.
"""
function _cg_column!(work::CGWork{V}, T::CGTables, J1::Int, J2::Int, J::Int, M::Int, lo::Int, n::Int) where {V}
    λ = _cgint(T, J + 1)^2
    meet = _cg_twosided!(work.w, work.a, work.ia, work.d, λ, n)
    cancel = 1.0
    T.exact || @inbounds for k in 1:n
        cancel = max(cancel, (work.dmag[k] + _mag(λ)) / _mag(work.d[k] - λ))
    end
    hi = lo + 2(n - 1)
    seedrel = 0.0
    if T.k >= 0
        # a level: scale from the certified endpoint C(hi), a single-term sum
        m, e, seedrel = _cg_level_seed(T.k, J1, hi, J2, M - hi, J)
        _cg_scale_to!(work.out, work.w, n, n, m, e)
    else
        _cg_normalise!(work.out, work.w, n, n)                # C(hi) > 0
    end
    _cg_estimates!(work.est, work.out, n, meet, seedrel, isfinite(cancel) ? cancel : Inf)
    return nothing
end

# ---------------------------------------------------------------------------------
#  The dual recurrence: rows in j
#
#  Cᵀ diag(q^{2m1}) C is tridiagonal in j (q^{2m1} plays the part of J1z), so at fixed (m1, m2) the coefficients
#  C(j) = ⟨j1 m1; j2 m2|j m⟩ solve
#
#    s_{j+1} C(j+1) − e_j C(j) + s_j C(j−1) = 0,
#    s_j = q^m √([j−m][j+m][j1+j2+1−j][j1+j2+1+j][j−j1+j2][j+j1−j2]) / ([2j] √([2j−1][2j+1])),
#    e_j = −q^{m1+j1}[j1−m1] + q^{m+j}[j−m][j+j1−j2][j1+j2+j+1]/([2j][2j+1])
#          + q^{m−j−1}[j1+j2−j][j+m+1][j−j1+j2+1]/([2j+1][2j+2]).
#
#  The q-Hahn recurrence (Koekoek–Swarttouw 14.6.3: p = q², n = j1+j2−j, N = j1+j2−m, α = p^{−2j1−1},
#  β = p^{−2j2−1}) in q-numbers with q − q⁻¹ divided out, so q → 1 is the classical recurrence.
# ---------------------------------------------------------------------------------

"""
    _cg_row_coefficients!(work, T, J1, M1, J2, M2, jlo, n)

Coefficients of the row recurrence for j = jlo, jlo+1, … over n values (doubled labels), in the form of
`_cg_twosided!` with λ = 0: a[k] = s_{j_k + 1}, d[k] = −e_{j_k}. All in r = q^{1/2}, [x]_q ∝ [2x]_r; the common
[2]_r drops out because every term has one more bracket upstairs than down.
"""
function _cg_row_coefficients!(work::CGWork{V}, T::CGTables, J1::Int, M1::Int, J2::Int, M2::Int, jlo::Int, n::Int) where {V}
    M = M1 + M2
    I(x) = _cgint(T, x)
    @inbounds for k in 1:n
        J = jlo + 2(k - 1)
        t1 = -_cgpow(T, M1 + J1) * I(J1 - M1)
        t2 = (J > 0 && J != M) ? _cgpow(T, M + J) * I(J - M) * I(J + J1 - J2) * I(J1 + J2 + J + 2) / (I(2J) * I(2J + 2)) : zero(V)
        t3 = J1 + J2 != J ? _cgpow(T, M - J - 2) * I(J1 + J2 - J) * I(J + M + 2) * I(J - J1 + J2 + 2) / (I(2J + 2) * I(2J + 4)) : zero(V)
        work.d[k] = -(t1 + t2 + t3)
        work.dmag[k] = _mag(t1) + _mag(t2) + _mag(t3)
        if k < n
            Jn = J + 2
            sv = _cgpow(T, M) * sqrt(I(Jn - M) * I(Jn + M) * I(J1 + J2 + 2 - Jn) * I(J1 + J2 + 2 + Jn) *
                                     I(Jn - J1 + J2) * I(Jn + J1 - J2)) / (I(2Jn) * sqrt(I(2Jn - 2) * I(2Jn + 2)))
            work.a[k] = sv; work.ia[k] = inv(sv)
        end
    end
    return nothing
end

"""
    _cg_row!(work, T, J1, M1, J2, M2, jlo, n) -> Bool

The row C(j) = ⟨j1 m1; j2 m2|j m⟩, j = jlo, jlo+1, …, in `work.out[1:n]` with error estimates in `work.est`.
At real q the row is complete and is normalised (Σ_j C² = 1, C(j1+j2) > 0); at a level the row is used only
when the fusion rule keeps every j (else `false`), and the scale comes from a single-term end — j = |j1−j2|
when that is the lowest j, else j = j1+j2.
"""
function _cg_row!(work::CGWork{V}, T::CGTables, J1::Int, M1::Int, J2::Int, M2::Int, jlo::Int, n::Int) where {V}
    jhi = jlo + 2(n - 1)
    # At a level the row must be complete. Where the fusion rule cuts it, the first coefficient beyond the cut
    # is not zero — for half-integer spins its [h] upstairs is cancelled by [2j+1] = [h] downstairs — so the
    # backward branch would start from a false boundary condition (measured: wrong top entries at k = 10).
    T.k >= 0 && jhi != J1 + J2 && return false
    _cg_row_coefficients!(work, T, J1, M1, J2, M2, jlo, n)
    meet = _cg_twosided!(work.w, work.a, work.ia, work.d, zero(DWNum), n)
    cancel = 1.0
    T.exact || @inbounds for k in 1:n
        cancel = max(cancel, work.dmag[k] / _mag(work.d[k]))
    end
    seedrel = 0.0
    if T.k >= 0
        if jlo == abs(J1 - J2)
            m, e, seedrel = _cg_level_seed(T.k, J1, M1, J2, M2, jlo)
            _cg_scale_to!(work.out, work.w, n, 1, m, e)
        else
            m, e, seedrel = _cg_level_seed(T.k, J1, M1, J2, M2, jhi)
            _cg_scale_to!(work.out, work.w, n, n, m, e)
        end
    else
        _cg_normalise!(work.out, work.w, n, n)                # the stretched coefficient is positive
    end
    _cg_estimates!(work.est, work.out, n, meet, seedrel, isfinite(cancel) ? cancel : Inf)
    return true
end

"One caller-owned CG sector; the numeric type is fixed, but the target and labels may change."
mutable struct CGSectorScratch{P,O}
    tables::CGTables{P}
    work::CGWork{P,O}
    J1::Int
    J2::Int
    M::Int
    lo::Int
    n::Int
end

function _cg_sector_workspace(T::CGTables, ::Nothing, J1::Int, J2::Int, M::Int)
    lo = max(-J1,M-J2); hi = min(J1,M+J2)
    work = CGWork(T,(hi-lo)÷2+1)
    lo,n = _cg_sector!(work,T,J1,J2,M)
    return work,lo,n
end

function _cg_sector_workspace(T::CGTables{P}, workspace::EvaluationWorkspace, J1::Int, J2::Int, M::Int) where {P}
    O = P === DWNum ? Float64 : ComplexF64
    saved = workspace.cg[]
    if saved isa CGSectorScratch{P,O}
        if saved.tables === T && saved.J1 == J1 && saved.J2 == J2 && saved.M == M
            return saved.work,saved.lo,saved.n
        end
        n = (min(J1,M+J2)-max(-J1,M-J2))÷2+1
        work = saved.work
        if length(work.w) < n
            capacity = max(n,2length(work.w))
            for v in (work.a,work.ia,work.d,work.dmag,work.w,work.out,work.est)
                resize!(v,capacity)
            end
        end
        lo,n = _cg_sector!(work,T,J1,J2,M)
        saved.tables=T; saved.J1=J1; saved.J2=J2; saved.M=M; saved.lo=lo; saved.n=n
        return work,lo,n
    end
    work,lo,n = _cg_sector_workspace(T,nothing,J1,J2,M)
    workspace.cg[] = CGSectorScratch(T,work,J1,J2,M,lo,n)
    return work,lo,n
end

"""
    _cg_entry(q, k, J1, M1, J2, M2, J; workspace=nothing) -> value or nothing

One coefficient from its column, for the near-edge tier, returned only if its estimate meets
`RTOL_PLAIN/8`. `q` real and positive, or a level. The row through the entry is the same length (a sector is
square), and the column can reuse workspace coefficients and is never cut by a level.
"""
function _cg_entry(q, k, J1::Int, M1::Int, J2::Int, M2::Int, J::Int; workspace=nothing)
    T = _cg_tables(q, k, J1, J2)
    work,lo,n = _cg_sector_workspace(T,workspace,J1,J2,M1+M2)
    _cg_column!(work, T, J1, J2, J, M1 + M2, lo, n)
    i = (M1 - lo) ÷ 2 + 1
    return work.est[i] <= RTOL_PLAIN / 8 ? work.out[i] : nothing
end

"Whether the column fast path applies: classical, a level, or real q > 0 (complex q stays entry by entry)."
_cg_column_target(q, k) = k !== nothing || q === nothing || _is_classical(q) ||
                          (q isa Real && !(q isa BigFloat) && q > 0)

"""
    qcg_matrix(j1, j2, m; k=nothing, q=nothing, T=nothing) -> (C, m1, j)

One sector of the coupling matrix: `C[a, b] = qcg(j1, m1[a], j2, m − m1[a], j[b], m)`, with `m1` running
down from its largest value and `j` up over the spins of j1 ⊗ j2 that reach m (inside the fusion rule at a
level). At positive real `q` and classically `C` is square and orthogonal, `transpose(C) * C ≈ I`.
At a level, the retained columns satisfy the same bilinear identity, but fusion truncation can make
`C` rectangular; then `C * transpose(C)` need not be the identity.

Built column by column from the Casimir recurrence in double words, `O(1)` work per entry, at real `q > 0`,
classically and at a level; every entry whose own error estimate misses the package's accuracy promise is
recomputed by [`qcg`](@ref). Complex or negative `q` and `BigFloat` entries are evaluated entry by entry.
`T` sets a precision floor, as in [`qcg`](@ref): level and negative-q values retain their complex type,
and a higher-precision `q` is never narrowed. For example, `k = 10, T = BigFloat` gives complex BigFloat entries.
"""
function qcg_matrix(j1::Spin, j2::Spin, m::Spin; k = nothing, q = nothing, T::Union{Type,Nothing} = nothing)
    J1, J2, M = doubled(j1, j2, m)
    kk = _fmatrix_level(k, q)
    iseven(J1 + J2 + M) && abs(M) <= J1 + J2 || throw(ArgumentError("m = $m is not a magnetic label of $j1 ⊗ $j2"))
    js = _qcg_js(J1, J2, M, kk)
    lo = max(-J1, M - J2); hi = min(J1, M + J2)
    E = _qcg_eltype(T, q, kk)
    C = zeros(E, (hi - lo) ÷ 2 + 1, length(js))
    if !isempty(js)
        jidx = Dict(J => b for (b, J) in enumerate(js))
        tw = _qcg_fast(E, q, kk) ? _qcg_tables_work(q, kk, J1, J2) : nothing
        _qcg_sector_columns!((M1, J, v) -> (C[(hi - M1) ÷ 2 + 1, jidx[J]] = v), E, q, kk, J1, J2, M, js, tw)
    end
    return C, collect(hi:-2:lo) .// 2, js .// 2
end

"""
    qcg_matrix(j1, j2; k=nothing, q=nothing, T=nothing) -> (C, product, coupled)

The coupling matrix of j1 ⊗ j2: `C[a, b] = qcg(j1, m1, j2, m2, j, m)` with `product[a] = (m1, m2)` and
`coupled[b] = (j, m)`. Rows follow `kron` order — m1 from j1 down to −j1 and, for each, m2 from j2 down — so a
column is the coupled vector |j m⟩ in the basis `kron(e_{m1}, e_{m2})`; columns run over j upwards and, within
each j, m from j down. At positive real `q` and classically `C` is orthogonal; at a level only the `j` inside
the fusion rule are columns, and `transpose(C) * C ≈ I` is a bilinear identity, without conjugation.

The matrix is dense, with (2j1+1)(2j2+1) rows and as many columns before fusion truncation (about 110 MB
for Float64 at j1 = j2 = 30). `qcg_matrix(j1, j2, m)` holds one sector.

```julia
using LinearAlgebra
C, prod, coup = qcg_matrix(1, 1//2; q = 0.8)
transpose(C) * C ≈ I
b = findfirst(==((3//2, 3//2)), coup)
C[1, b] ≈ qcg(1, 1, 1//2, 1//2, 3//2, 3//2; q = 0.8)
```
"""
function qcg_matrix(j1::Spin, j2::Spin; k = nothing, q = nothing, T::Union{Type,Nothing} = nothing)
    J1, J2 = doubled(j1, j2)
    kk = _fmatrix_level(k, q)
    prod = [(M1 // 2, M2 // 2) for M1 in J1:-2:-J1 for M2 in J2:-2:-J2]
    jall = [J for J in abs(J1 - J2):2:(J1 + J2) if kk === nothing || (J1 + J2 + J) ÷ 2 <= kk]
    coupled = [(J // 2, M // 2) for J in jall for M in J:-2:-J]
    E = _qcg_eltype(T, q, kk)
    C = zeros(E, length(prod), length(coupled))
    isempty(jall) && return C, prod, coupled
    # column of (J, M): the J block starts after Σ_{J' < J} (J' + 1) columns, and m runs down within it
    offset = Dict{Int,Int}(); acc = 0
    for J in jall; offset[J] = acc; acc += J + 1; end
    tw = _qcg_fast(E, q, kk) ? _qcg_tables_work(q, kk, J1, J2) : nothing
    for M in -(J1 + J2):2:(J1 + J2)
        js = _qcg_js(J1, J2, M, kk)
        isempty(js) && continue
        _qcg_sector_columns!(E, q, kk, J1, J2, M, js, tw) do M1, J, v
            C[((J1 - M1) ÷ 2) * (J2 + 1) + (J2 - (M - M1)) ÷ 2 + 1, offset[J] + (J - M) ÷ 2 + 1] = v
        end
    end
    return C, prod, coupled
end

"The spins of the sector m (doubled), inside the fusion rule at a level."
_qcg_js(J1, J2, M, kk) = [J for J in abs(J1 - J2):2:(J1 + J2) if J >= abs(M) && (kk === nothing || (J1 + J2 + J) ÷ 2 <= kk)]

function _qcg_eltype(T, q, kk)
    # As in the scalar API, T is a precision floor: it cannot discard the
    # imaginary part or narrow a higher-precision parameter supplied by the caller.
    E = T === nothing ? Float64 : T
    q === nothing || (E = promote_type(E, typeof(float(q))))
    return kk !== nothing || (q isa Real && q < 0) ? Complex{real(E)} : E
end

_qcg_fast(::Type{E}, q, kk) where {E} = real(E) === Float64 && _cg_column_target(q, kk)

function _qcg_tables_work(q, kk, J1, J2)
    T = _cg_tables(q, kk, J1, J2)
    return T, CGWork(T, min(J1,J2) + 1)
end

"""
Every column of the sector m, handed to `put(M1, J, value)` (doubled labels): by the recurrence when `tw`
holds tables and scratch, entry by entry through `qcg` otherwise or for an entry whose estimate misses the
promise.
"""
function _qcg_sector_columns!(put::F, ::Type{E}, q, kk, J1::Int, J2::Int, M::Int, js, tw) where {F,E}
    j1 = J1 // 2; j2 = J2 // 2
    # `::E`: without it the fallback's type is not inferred, the value handed to `put` is boxed, and every
    # entry allocates — although no entry falls back (5,856 allocations per (30,30) sector, 1.3× the time)
    direct(M1, J) = E(kk === nothing ? qcg(j1, M1 // 2, j2, (M - M1) // 2, J // 2; q = q, T = real(E)) :
                                       qcg(j1, M1 // 2, j2, (M - M1) // 2, J // 2; k = kk, T = real(E)))::E
    lo = max(-J1, M - J2); hi = min(J1, M + J2)
    if tw === nothing
        for J in js, M1 in hi:-2:lo
            put(M1, J, direct(M1, J))
        end
        return nothing
    end
    T, work = tw
    lo, n = _cg_sector!(work, T, J1, J2, M)
    for J in js
        _cg_column!(work, T, J1, J2, J, M, lo, n)
        @inbounds for i in 1:n
            M1 = lo + 2(i - 1)
            put(M1, J, work.est[i] <= RTOL_PLAIN ? E(work.out[i]) : direct(M1, J))
        end
    end
    return nothing
end

"""
    qcg_row(j1, m1, j2, m2; k=nothing, q=nothing, T=nothing) -> (c, j)

The product state |j1 m1⟩⊗|j2 m2⟩ in the coupled basis: `c[i] = qcg(j1, m1, j2, m2, j[i])` for every `j` of
j1 ⊗ j2 that reaches m = m1 + m2 (inside the fusion rule at a level), `j` ascending. It is one row of
[`qcg_matrix`](@ref), restricted to its allowed coupling channels (some coefficients can vanish).
At positive real `q` and classically `sum(abs2, c) ≈ 1`.

At real `q > 0`, classically and at a level the row comes from the three-term recurrence in j that
q^{2m1} satisfies in the coupled basis (the q-analogue of J1z, from the q-Hahn polynomials), in double words:
`O(1)` work per entry. Entries whose error estimate misses the package's promise are recomputed by
[`qcg`](@ref); complex or negative `q` and `BigFloat` are evaluated entry by entry.
`T` is a precision floor; `k = 10, T = BigFloat`, for example, returns complex BigFloat entries.
"""
function qcg_row(j1::Spin, m1::Spin, j2::Spin, m2::Spin; k = nothing, q = nothing, T::Union{Type,Nothing} = nothing)
    J1, M1, J2, M2 = doubled(j1, m1, j2, m2)
    kk = _fmatrix_level(k, q)
    (abs(M1) <= J1 && abs(M2) <= J2 && iseven(J1 + M1) && iseven(J2 + M2)) ||
        throw(ArgumentError("($m1, $m2) are not magnetic labels of $j1 ⊗ $j2"))
    M = M1 + M2
    js = _qcg_js(J1, J2, M, kk)
    E = _qcg_eltype(T, q, kk)
    c = Vector{E}(undef, length(js))
    isempty(js) && return c, js .// 2
    # `::E` keeps the fallback inferred, so the row does not box every entry (as in `_qcg_sector_columns!`)
    direct(J) = E(kk === nothing ? qcg(j1, m1, j2, m2, J // 2; q = q, T = real(E)) :
                                   qcg(j1, m1, j2, m2, J // 2; k = kk, T = real(E)))::E
    if _qcg_fast(E, q, kk)
        Tb = _cg_tables(q, kk, J1, J2)
        work = CGWork(Tb, length(js))
        if _cg_row!(work, Tb, J1, M1, J2, M2, first(js), length(js))
            for i in eachindex(js)
                c[i] = work.est[i] <= RTOL_PLAIN ? E(work.out[i]) : direct(js[i])
            end
            return c, js .// 2
        end
    end
    for (i, J) in enumerate(js)
        c[i] = direct(J)
    end
    return c, js .// 2
end
