# Exact level values, continued from exact_display.jl: arithmetic, `ExactXSum`, scalars and exact equality.

# ---------------------------------------------------------------------------------
#  Arithmetic
#
#  A value is `P(x)·√(R(x))` with `R = ∏_{e ∈ rad} ψ_e` reduced modulo `Ψ_h`, and the convention that
#  `√` is the **positive** root — which is forced, since the value is real and `P` is.
#
#  Multiplication is where the square classes do their work. With `A = ∏_S ψ` and `B = ∏_T ψ`,
#
#      A·B = G²·D,     G = ∏_{S ∩ T} ψ,     D = ∏_{S △ T} ψ,
#
#  so `√A·√B = √(G²D) = |G|·√D`: the shared factors leave the root, the symmetric difference stays, and
#  the class arithmetic is 𝔽₂ over the ψ basis exactly as in `generic_x.jl`. The one thing that is *not*
#  inherited from the generic case is the absolute value — over ℚ(x) there is no sign to take, but at a
#  level `G` is a number and `|G|` is `±G` according to it. Getting that wrong is a sign error in every
#  product, so the sign is taken from a certified evaluation rather than assumed.
#
#  Addition cannot fuse classes, so a sum is kept keyed by class, like `XSum` upstream and
#  `CompositeExactResult` in the cyclotomic layer. Unlike `XSum`, the keys here are **not** guaranteed
#  independent: a class that is squarefree in ℚ(x) can become a square modulo `Ψ_h`, so an empty term
#  list proves zero but a surviving coefficient proves nothing. `iszero` therefore asks for a proof —
#  the norm over all sign choices of the roots, which lands in the base field — exactly as
#  `types_projection.jl` does, and at half the degree.
# ---------------------------------------------------------------------------------

const _PSIQ_RED = Dict{Tuple{Int,Int},Any}()

"`ψ_e` reduced modulo `Ψ_h`, over ℚ and cached: the radicals of a level are built from these repeatedly."
function _psiq_red(e::Int, h::Int)
    lock(X_LOCK) do
        get!(_PSIQ_RED, (h, e)) do
            _redq(_toqq(psi_x(e)), _psiq(h))
        end
    end
end

"`∏_{e ∈ idx} ψ_e` reduced modulo `Ψ_h`."
function _psi_prod(idx, h::Int)
    S, _ = _qqx()
    isempty(idx) && return S(1)
    Ψ = _psiq(h)
    acc = S(1)
    for e in idx
        acc = _redq(acc * _psiq_red(e, h), Ψ)
    end
    return acc
end

"""
Sign of a nonzero field element at `x = 2cos(π/h)`, from a certified evaluation. A product of two values
needs `|G|`, and `G` is an exact nonzero element, so this is a decision and not an estimate — it throws
rather than guess if the evaluation cannot be certified.
"""
function _field_sign(f, h::Int)
    iszero(f) && return 0
    v = _eval_at_x(f, h; rtol = 1e-8)
    v === nothing && throw(DomainError(h,
        "could not certify the sign of a radical's shared factor at level $(h - 2); " *
        "the exact values are unaffected, but their product cannot be formed"))
    return v < 0 ? -1 : 1
end

const _PSI_SIGN = Dict{Tuple{Int,Int},Int}()

"""
Sign of `ψ_e(2cos(π/h))`, cached. `sign(G) = ∏ sign(ψ_e)` for a product, so the per-multiplication cost
of `|G|` is a dictionary lookup rather than a certified evaluation.

Measured, this is insurance and not a fix: over every radical class of every admissible symbol at
k = 6, 8, 10, 12 — 18 distinct ψ indices and 191,514 pairs sharing a class — **no** `ψ_e` was negative and
no shared factor was. That is not an accident (the indices that occur divide `2n` for factorial arguments
`n < h`), but it is also not a theorem anyone has written down here, so the absolute value stays.
"""
function _psi_sign(e::Int, h::Int)
    lock(X_LOCK) do
        get!(_PSI_SIGN, (e, h)) do
            _field_sign(_psi_prod([e], h), h)
        end
    end
end

Base.:-(v::ExactX) = ExactX(v.k, -v.p, v.sqclass, v.r)
function _one_exactx(k::Integer)
    S, _ = _qqx()
    return ExactX(Int(k), S(1), Int[], S(1))
end
Base.one(::Type{ExactX}, k::Integer) = _one_exactx(k)

"""
    _qfact_exactx(pairs, k; sgn = 1) -> ExactX

`sgn · ∏ₙ [n]!^{cₙ}` at level `k`, from a collection of `n => c`. This is the whole content of every
value in the package that has no sum — a quantum dimension, a theta net, any `CyclotomicMonomial` whose
q-power cancels — so those reach the real basis by the same route the prefactor of a Racah sum does,
through the ψ exponents, with no cyclotomic field anywhere.
"""
function _qfact_exactx(pairs, k::Integer; sgn::Int = 1)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
    h = kk + 2
    Ψ = _psiq(h)
    S, _ = _qqx()
    iszero(sgn) && return zero(ExactX, kk)
    num, den = _psi_monomial(psi_exponents(pairs); modulus = psi_level(h), level = h)
    d = _redq(_toqq(den), Ψ)
    iszero(d) && throw(DomainError(kk, "the value is singular at level $kk: a q-integer vanishes"))
    n = _redq(_toqq(num) * QQ(sgn), Ψ)
    iszero(n) && return zero(ExactX, kk)
    p = _divmod_psi(n, d, h)
    p === nothing && (p = _redq(n * _invmodq(d, Ψ), Ψ))
    return ExactX(kk, p, Int[], S(1))
end

"""
    _qint_exactx(k, n) -> ExactX

The q-integer `[n]` at level `k` as an exact value. Identities carry these as coefficients — the
`[x+1]` of Biedenharn–Elliott, say — and they belong in the same arithmetic as the symbols they multiply.
"""
function _qint_exactx(k::Integer, n::Integer)
    S, _ = _qqx()
    kk = Int(k)
    p = _redq(_toqq(qint_x(Int(n))), _psiq(kk + 2))
    return ExactX(kk, p, Int[], S(1))
end
Base.one(v::ExactX) = _one_exactx(v.k)

_same_level(a, b) = a.k == b.k ||
    throw(ArgumentError("level mismatch: $(a.k) and $(b.k) are values at different roots of unity"))

"""
    a * b

The product, with the shared square classes leaving the root. Closed on single values: two classes always
fuse to one.
"""
function Base.:*(a::ExactX, b::ExactX)
    _same_level(a, b)
    (iszero(a.p) || iszero(b.p)) && return zero(ExactX, a.k)
    h = a.k + 2
    Ψ = _psiq(h)
    inter = intersect(a.sqclass, b.sqclass)
    diff = sort!(symdiff(a.sqclass, b.sqclass))
    G = _psi_prod(inter, h)
    iszero(G) && return zero(ExactX, a.k)          # a radical that vanishes at this level
    sg = prod(e -> _psi_sign(e, h), inter; init = 1)
    p = _redq(_redq(a.p * b.p, Ψ) * (sg < 0 ? -G : G), Ψ)
    iszero(p) && return zero(ExactX, a.k)
    return ExactX(a.k, p, diff, _psi_prod(diff, h))
end

Base.:*(c::Union{Integer,Rational}, v::ExactX) =
    iszero(c) ? zero(ExactX, v.k) : ExactX(v.k, v.p * QQ(c), v.sqclass, v.r)
Base.:*(v::ExactX, c::Union{Integer,Rational}) = c * v
function Base.:/(v::ExactX, c::Union{Integer,Rational})
    iszero(c) && throw(DivideError())
    return (1 // c) * v
end

"""
    inv(v)

`1/(P√R) = (1/(P·R))·√R`: the root moves to the numerator against itself, so the class is unchanged and
no sign decision is needed.
"""
function Base.inv(v::ExactX)
    iszero(v.p) && throw(DivideError())
    h = v.k + 2
    Ψ = _psiq(h)
    d = _redq(v.p * v.r, Ψ)
    iszero(d) && throw(DivideError())
    S, _ = _qqx()
    q = _divmod_psi(S(1), d, h)
    q === nothing && (q = _invmodq(d, Ψ))
    return ExactX(v.k, q, copy(v.sqclass), v.r)
end
Base.:/(a::ExactX, b::ExactX) = a * inv(b)
Base.:/(c::Union{Integer,Rational}, v::ExactX) = c * inv(v)

"""
    is_provably_nonzero(v::ExactX) -> Bool

Whether the value is nonzero by an exact argument. A single value carries one square class, so its
coefficient decides — there is nothing here for the norm of [`is_provably_nonzero(::ExactXSum)`](@ref)
to do, and that is the whole difference between a value and a sum of them.
"""
is_provably_nonzero(v::ExactX) =
    !iszero(v.p) && (isempty(v.sqclass) || !iszero(v.r))

function Base.:^(v::ExactX, n::Integer)
    n < 0 && return inv(v)^(-n)
    n == 0 && return _one_exactx(v.k)
    acc = v
    for _ in 2:n
        acc = acc * v
    end
    return acc
end

# ---------------------------------------------------------------------------------
#  Sums, keyed by square class
# ---------------------------------------------------------------------------------

"""
    ExactXSum

A sum `Σ_S c_S(x)·√(∏_{e ∈ S} ψ_e(x))` at one level, keyed by the square class `S`. Products fuse classes
and stay exact; sums cannot, so they are kept apart.

The keys are the classes of `generic_x.jl`, squarefree over ℚ(x) — but a class that is squarefree there
can become a **square** modulo `Ψ_h`, so two keys may denote the same root and `isempty(terms)` is the only
immediate proof of vanishing. `iszero` first uses the structural and norm tests, then resolves
dependent specialized classes with exact algebraic numbers when necessary.
"""
struct ExactXSum
    k::Int
    terms::Dict{Vector{Int},QQPolyRingElem}
end

ExactXSum(k::Integer) = ExactXSum(Int(k), Dict{Vector{Int},QQPolyRingElem}())
function ExactXSum(v::ExactX)
    iszero(v.p) && return ExactXSum(v.k)
    return ExactXSum(v.k, Dict(copy(v.sqclass) => v.p))
end
Base.zero(::Type{ExactXSum}, k::Integer) = ExactXSum(k)
Base.zero(s::ExactXSum) = ExactXSum(s.k)
Base.length(s::ExactXSum) = length(s.terms)
Base.isempty(s::ExactXSum) = isempty(s.terms)
level(s::ExactXSum) = s.k

"The value of one class as an `ExactX`."
_term(s::ExactXSum, S) = ExactX(s.k, s.terms[S], copy(S), _psi_prod(S, s.k + 2))

function Base.:+(a::ExactXSum, b::ExactXSum)
    _same_level(a, b)
    out = copy(a.terms)
    for (S, c) in b.terms
        if haskey(out, S)
            w = out[S] + c
            iszero(w) ? delete!(out, S) : (out[S] = w)
        else
            out[S] = c
        end
    end
    return ExactXSum(a.k, out)
end
Base.:-(a::ExactXSum) = ExactXSum(a.k, Dict(S => -c for (S, c) in a.terms))
Base.:-(a::ExactXSum, b::ExactXSum) = a + (-b)

function Base.:*(a::ExactXSum, b::ExactXSum)
    _same_level(a, b)
    out = ExactXSum(a.k)
    for S in keys(a.terms), T in keys(b.terms)
        out = out + ExactXSum(_term(a, S) * _term(b, T))
    end
    return out
end
Base.:*(c::Union{Integer,Rational}, s::ExactXSum) =
    iszero(c) ? zero(s) : ExactXSum(s.k, Dict(S => v * QQ(c) for (S, v) in s.terms))
Base.:*(s::ExactXSum, c::Union{Integer,Rational}) = c * s

function Base.:^(s::ExactXSum, n::Integer)
    n >= 0 || throw(ArgumentError("negative powers of a sum are not supported"))
    acc = ExactXSum(_one_exactx(s.k))
    for _ in 1:n
        acc = acc * s
    end
    return acc
end

for op in (:+, :-, :*)
    @eval Base.$op(a::ExactX, b::ExactXSum) = $op(ExactXSum(a), b)
    @eval Base.$op(a::ExactXSum, b::ExactX) = $op(a, ExactXSum(b))
end
"""
A sum that carries one square class **is** a value, and is handed back as one.

Without this, `v * 2` returns an `ExactX` and `v + 2` an `ExactXSum` — the same number in two types,
differing by which method the user happened to call. Adding two values whose classes genuinely differ
still gives an [`ExactXSum`](@ref), because that is a thing the value type cannot hold; and arithmetic
that starts in `ExactXSum` stays there, so the sum type is never taken away from code that asked for it.
"""
_collapse(s::ExactXSum) = isempty(s.terms) ? zero(ExactX, s.k) :
                          length(s.terms) == 1 ? _term(s, first(keys(s.terms))) : s

Base.:+(a::ExactX, b::ExactX) = _collapse(ExactXSum(a) + ExactXSum(b))
Base.:-(a::ExactX, b::ExactX) = _collapse(ExactXSum(a) - ExactXSum(b))

# ---------------------------------------------------------------------------------
#  Scalars, on the same footing as values
# ---------------------------------------------------------------------------------

"""
A rational as an exact value at a level: the constant polynomial, no class, no root. This is what lets
`v + 1` mean what it says — a scalar is a value like any other, and refusing to add one was an accident
of which methods happened to be written.
"""
function _const_exactx(c::Union{Integer,Rational}, k::Integer)
    S, _ = _qqx()
    return ExactX(Int(k), S(QQ(c)), Int[], S(1))
end

for T in (:ExactX, :ExactXSum)
    # `v + c` goes through the value/value method for an `ExactX`, so it collapses like any other sum
    @eval Base.:+(v::$T, c::Union{Integer,Rational}) = v + _const_exactx(c, v.k)
    @eval Base.:+(c::Union{Integer,Rational}, v::$T) = v + _const_exactx(c, v.k)
    @eval Base.:-(v::$T, c::Union{Integer,Rational}) = v + _const_exactx(-c, v.k)
    @eval Base.:-(c::Union{Integer,Rational}, v::$T) = (-v) + _const_exactx(c, v.k)
    @eval Base.:(==)(v::$T, c::Union{Integer,Rational}) = v == _const_exactx(c, v.k)
    @eval Base.:(==)(c::Union{Integer,Rational}, v::$T) = v == _const_exactx(c, v.k)
end

"""
Floats do not mix with exact values, and the error says what to do instead.

Accepting one would be worse than refusing: `Float64` is exactly a dyadic rational, so `v + 0.1` would
silently commit to `3602879701896397//36028797018963968` and print as though it were the value the user
meant. Either the scalar is exact, and `1//10` says so, or the computation is numerical, and `Float64(v)`
says that.
"""
function _inexact_mix(x)
    throw(ArgumentError(
        "an exact value does not combine with $(typeof(x)): pass an exact scalar (for example `1//10` " *
        "rather than `0.1`), or leave the exact field first with `Float64(v)`"))
end
for op in (:+, :-, :*, :/)
    @eval Base.$op(::Union{ExactX,ExactXSum}, x::AbstractFloat) = _inexact_mix(x)
    @eval Base.$op(x::AbstractFloat, ::Union{ExactX,ExactXSum}) = _inexact_mix(x)
    @eval Base.$op(::Union{ExactX,ExactXSum}, x::Complex) = _inexact_mix(x)
    @eval Base.$op(x::Complex, ::Union{ExactX,ExactXSum}) = _inexact_mix(x)
end

Base.one(::Type{ExactXSum}, k::Integer) = ExactXSum(_one_exactx(k))
Base.one(s::ExactXSum) = one(ExactXSum, s.k)
Base.isone(v::ExactX) = isempty(v.sqclass) && isone(v.p)
Base.isone(s::ExactXSum) = length(s.terms) == 1 && isone(_term(s, first(keys(s.terms))))

"""
    radical_norm(s::ExactXSum) -> field element or nothing

`∏` over all sign choices of the roots, which lies in ℚ[x]/Ψ_h because every root then appears to an even
total power. `nothing` when the elimination has not closed after a few steps, which callers must treat as
undecided rather than as an answer.
"""
function radical_norm(s::ExactXSum)
    isempty(s.terms) && return nothing
    cur = s
    for _ in 1:6
        ks = collect(keys(cur.terms))
        if length(ks) == 1
            isempty(ks[1]) || return nothing        # still under a root: not in the base field
            return cur.terms[ks[1]]
        end
        pick = ks[end]
        conj = ExactXSum(cur.k, Dict(S => (S == pick ? -c : c) for (S, c) in cur.terms))
        cur = cur * conj
        isempty(cur.terms) && return nothing
    end
    return nothing
end

"""
    is_provably_nonzero(s::ExactXSum) -> Bool

Whether the value is nonzero **by an exact argument**, never by assuming that distinct class keys are
independent. One term is nonzero when its coefficient and its radical are; otherwise a nonzero
[`radical_norm`](@ref) is the proof. `false` means undecided as well as zero, so it is safe to branch on.
"""
function is_provably_nonzero(s::ExactXSum)
    isempty(s.terms) && return false
    if length(s.terms) == 1
        S, c = first(s.terms)
        iszero(c) && return false
        return isempty(S) || !iszero(_psi_prod(S, s.k + 2))
    end
    N = radical_norm(s)
    return N !== nothing && !iszero(N)
end

"""
Whether the sum is exactly zero. Structural zeros and a nonzero radical norm settle the common cases.
When specialized classes are dependent, Nemo's algebraic numbers select the real embedding
`x = 2cos(π/(k+2))` and the nonnegative square roots exactly. No numerical tolerance decides equality.
This fallback can cost more for large algebraic degrees; numerical symbol evaluation never uses it.
"""
function Base.iszero(s::ExactXSum)
    isempty(s.terms) && return true
    if length(s.terms) == 2
        # Different generic classes often reduce to the very same radicand at a level.
        # Merge that pair in the base field before either the norm or the algebraic fallback.
        (S, a), state = iterate(s.terms)
        (T, b), _ = iterate(s.terms, state)
        r = _psi_prod(S, s.k + 2)
        if r == _psi_prod(T, s.k + 2)
            return iszero(r) || iszero(a + b)
        end
    end
    is_provably_nonzero(s) && return false
    length(s.terms) == 1 && return true              # its coefficient or its radical vanishes
    return iszero(_algebraic_value(s))
end

"Resolve an undecided exact sum in its real embedding, including dependencies between radical classes."
function _algebraic_value(s::ExactXSum)
    h = s.k + 2
    x = 2 * cospi(QQBar(1 // h))
    acc = QQBar(0)
    for (S, c) in s.terms
        p = evaluate(c, x)
        iszero(p) && continue
        if !isempty(S)
            r = evaluate(_psi_prod(S, h), x)
            r < 0 && throw(DomainError(s.k, "the exact sum has a negative real radicand"))
            p *= sqrt(r)
        end
        acc += p
    end
    return acc
end

Base.:(==)(a::ExactXSum, b::ExactXSum) = a.k == b.k && iszero(a - b)
Base.:(==)(a::ExactX, b::ExactX) = ExactXSum(a) == ExactXSum(b)
Base.:(==)(a::ExactX, b::ExactXSum) = ExactXSum(a) == b
Base.:(==)(a::ExactXSum, b::ExactX) = a == ExactXSum(b)

"""
    numeric_value(s::ExactXSum) -> BigFloat or nothing

The number the sum denotes, each term certified as in [`numeric_value(::ExactX)`](@ref).
"""
function numeric_value(s::ExactXSum)
    isempty(s.terms) && return BigFloat(0)
    acc = BigFloat(0)
    for S in keys(s.terms)
        t = numeric_value(_term(s, S))
        t === nothing && return nothing
        acc += t
    end
    return acc
end
function Base.float(s::ExactXSum, ::Type{T} = Float64) where {T<:AbstractFloat}
    v = numeric_value(s)
    v === nothing && throw(DomainError(s.k,
        "the sum's value could not be certified numerically; see `numeric_value`"))
    return T(v)
end
Base.Float64(s::ExactXSum) = float(s, Float64)
evaluate_exact(s::ExactXSum, ::Type{T} = ComplexF64) where {T} = T(float(s, Float64))

function Base.show(io::IO, s::ExactXSum)
    isempty(s.terms) && return print(io, "ExactXSum(k = ", s.k, ", 0)")
    parts = [_xform_str(_term(s, S)) for S in sort!(collect(keys(s.terms)))]
    print(io, "ExactXSum(k = ", s.k, ", ", join(parts, " + "), ")")
end

"""
The same display an [`ExactX`](@ref) gets, one line per square class: no arithmetic, the shape of a term
that is too long to write out, and the `≈` only when `IOContext(io, :approximate => true)` asks for it.
A sum and a value differ in how many classes they carry, not in how they are read — and a sum of one
class is read exactly like the value it is.
"""
function Base.show(io::IO, ::MIME"text/plain", s::ExactXSum)
    h = s.k + 2
    d = euler_phi(2h) ÷ 2
    # one class reads exactly like a single value, because that is what it is; the count appears only
    # when there is something to count
    head = "Exact value at level k = " * string(s.k) *
           (length(s.terms) > 1 ? ", " * string(length(s.terms)) * " square classes" : "") *
           "   (x = 2cos(π/" * string(h) * "), degree " * string(d) * " over ℚ)"
    lines = [head]
    if isempty(s.terms)
        print(io, lines[1], "\n  = 0")
        return
    end
    long = false
    for (i, S) in enumerate(sort!(collect(keys(s.terms))))
        t = _term(s, S)
        push!(lines, (i == 1 ? "  = " : "  + ") * _xform_str(t))
        if _xform_long(t)
            long = true
            for l in _xform_dims(t)
                push!(lines, "      " * l)
            end
        end
    end
    long && push!(lines, "    `v.x_value` gives each (P, R) in full")
    if get(io, :approximate, false)
        v = numeric_value(s)
        push!(lines, v === nothing ? "  (the value could not be certified numerically)" :
                     "  ≈ " * string(Float64(v)))
    end
    print(io, join(lines, "\n"))
end

"""
Properties beyond the stored fields, the same two an [`ExactX`](@ref) has.

`s.x_value` is the vector of `(P, R)` pairs, one per square class, in the order the display lists them;
`s.rad` is [`radical(s)`](@ref radical).
"""
function Base.getproperty(s::ExactXSum, f::Symbol)
    if f === :x_value
        t = getfield(s, :terms)
        return [(P = t[S], R = _psi_prod(S, getfield(s, :k) + 2)) for S in sort!(collect(keys(t)))]
    end
    f === :rad && return radical(s)
    return getfield(s, f)
end
Base.propertynames(::ExactXSum, private::Bool = false) =
    private ? (:k, :terms, :x_value, :rad) : (:k, :terms, :x_value, :rad)

"""
    radical(s::ExactXSum; maxlen = 0, degree_limit = 32) -> RadExpr or NoRadical

Term by term, added up. A sum of square classes has a radical form exactly when each of its terms does,
and the first term that does not is the answer — with the *sum's* number under it, not the term's.
"""
function radical(s::ExactXSum; maxlen::Int = 0, degree_limit::Int = RADICAL_DESCENT_MAX_DEGREE)
    isempty(s.terms) && return RadExpr(0)
    acc = RadExpr(Rational{BigInt}(0))
    for S in sort!(collect(keys(s.terms)))
        e = radical(_term(s, S); maxlen = maxlen, degree_limit = degree_limit)
        if e isa NoRadical
            nv = numeric_value(s)
            return NoRadical(s.k, e.kind, e.degree, e.field, e.limit,
                             nv === nothing ? nothing : Float64(nv))
        end
        acc = RadExpr(acc.rat + e.rat, vcat(acc.terms, e.terms))
    end
    return simplify(acc)
end
