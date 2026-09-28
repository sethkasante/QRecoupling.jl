# ---------------------------------------------------------------------------------
#  Exact level values in x = q + q⁻¹, and the display ladder
#
#  A level-k value is real, but ℚ(ζ₂ₕ) does not show that: `−2/3 ζ⁶ + 4/3 ζ² − 1` is a real number whose
#  display contains an `i` (ζ₂₄⁶ = i). The real subfield does show it. Writing `x = q + q⁻¹ = 2cos(π/h)`,
#  `h = k+2`, every value of the package is
#
#      v = P(x) · √(R(x)),      P, R ∈ ℚ[x] reduced modulo Ψ_h,
#
#  with `Ψ_h` the minimal polynomial of `2cos(π/h)` — degree `φ(2h)/2`, **half** the cyclotomic degree.
#  This is the same shape as the generic `XValue` of `generic_x.jl`, so a level is not a second
#  architecture: it is *reduction modulo Ψ_h*, applied while the Racah sum is formed rather than after.
#
#  On top of that sits a display ladder, tried in order and stopping at the first form that is short:
#
#      rational  →  single surd  →  nested radical  →  polynomial in x  →  (the generic Φ-form in q)
#
#  How far the radicals go is not a matter of effort. `v²` lies in the real cyclotomic field, which is
#  abelian of degree `φ(2h)/2`; an abelian field is a tower of quadratic extensions **iff its degree is a
#  power of two**. So a nested-square-root form exists exactly when `φ(2k+4)/2` is a power of 2 —
#  k = 0,1,2,3,4,6,8,10,13,14,15,18,22,… — and for k = 5,7,9,11,12,16,17,19,… the degree carries an odd
#  prime factor and no real radical form exists at all (casus irreducibilis). `has_radical_form` decides
#  this in one line, and the display says so rather than failing quietly.
#
#  The descent itself is Lagrange's, not a general solver: with σ an automorphism whose square fixes u,
#
#      u = (u + σu)/2 + √( ((u − σu)/2)² ),
#
#  and both parts lie in the fixed field of ⟨Stab(u), σ⟩, which is strictly larger — so the recursion
#  terminates. The Galois action needs no number field: σ_j is the substitution x ↦ C_j(x), because
#  `C_j(2cos(π/h)) = 2cos(jπ/h)` runs over the conjugates as j runs over the residues coprime to 2h.
#
#  Measured coverage (`dev/results/exact_display_ladder.md`): the radical form is the whole story at
#  k ≤ 4, about half of it at k = 6–10, and 0.01% at k = 14. That is why the ladder falls through on
#  length rather than promising radicals everywhere.
# ---------------------------------------------------------------------------------

const _QQX = Ref{Any}()

"The polynomial ring ℚ[x]; the level values live in its quotient by Ψ_h."
function _qqx()
    isassigned(_QQX) && return _QQX[]::Tuple{QQPolyRing,QQPolyRingElem}
    lock(X_LOCK) do
        isassigned(_QQX) || (_QQX[] = polynomial_ring(QQ, "x"))
    end
    return _QQX[]::Tuple{QQPolyRing,QQPolyRingElem}
end

function _toqq(f::ZZPolyRingElem)
    S, _ = _qqx()
    return S([QQ(coeff(f, i)) for i in 0:max(degree(f), 0)])
end

"Ψ_h over ℚ, cached alongside its degree."
const _PSIQ = Dict{Int,Any}()
function _psiq(h::Int)
    lock(X_LOCK) do
        get!(_PSIQ, h) do
            _toqq(psi_level(h))
        end
    end
end

_redq(f, Ψ) = mod(f, Ψ)

"""
Word-sized primes for the modular division below. Sixty bits each, so two of them already carry more than
the answers ever need; the list is long only so that a run of unlucky primes cannot exhaust it.
"""
const _MM_PRIMES = let out = UInt[], p = ZZ(2)^60
    for _ in 1:40
        p = next_prime(p)
        push!(out, UInt(p))
    end
    out
end

"`f` as an integer polynomial together with the common denominator cleared out of it."
function _clear_denoms(f)
    d = Int(degree(f))
    R, _ = xring()
    d < 0 && return R(0), ZZ(1)
    den = ZZ(1)
    for i in 0:d
        den = lcm(den, denominator(coeff(f, i)))
    end
    return R([numerator(coeff(f, i) * den) for i in 0:d]), den
end

"""
    _divmod_psi(num, den, h) -> QQPoly or nothing

`num · den⁻¹ mod Ψ_h`, from images modulo word-sized primes, or `nothing` if the primes run out.

Why not simply invert. The inputs are enormous and the answer is not: at j = 20, k = 60 the denominator
carries 273-bit coefficients and the result 23-bit ones, and at j = 30, k = 100 it is 501 against 42. A
Euclidean inversion over ℚ[x] does all its work at the input's size — 1.4 ms there, and 3.2 ms at
j = 25, k = 80 — while the modular route works at 60 bits per prime and needs only enough primes to
reconstruct the *output*. Measured, that is **13–23× faster** at those sizes and a wash at small ones.

It is not a heuristic. Each pass adds a prime, CRTs, attempts rational reconstruction, and then **checks
`p·den ≡ num (mod Ψ)` exactly**; nothing is returned until that holds, so a wrong reconstruction cannot
escape and no coefficient bound has to be proved in advance.
"""
function _divmod_psi(num, den, h::Int)
    Ψ = _psiq(h)
    iszero(den) && return nothing
    # Rational denominators need no modular images, CRT or reconstruction.
    degree(den) == 0 && return _redq(num / coeff(den, 0), Ψ)
    S, _ = _qqx()
    d = Int(degree(Ψ))
    d <= 0 && return nothing
    Z, _ = _clear_denoms(Ψ)                    # Ψ is monic, so this is Ψ itself over ℤ
    N, aN = _clear_denoms(num)
    D, aD = _clear_denoms(den)
    iszero(D) && return nothing
    scal = QQ(aD, aN)
    M = ZZ(1)
    acc = ZZRingElem[]
    for pp in _MM_PRIMES
        F = Nemo.Native.GF(pp)
        Fx, _ = polynomial_ring(F, "X")
        Ψp = Fx([F(coeff(Z, i)) for i in 0:Int(degree(Z))])
        Dp = Fx([F(coeff(D, i)) for i in 0:max(Int(degree(D)), 0)])
        isone(gcd(Dp, Ψp)) || continue         # den is not invertible here: try another prime
        Np = Fx([F(coeff(N, i)) for i in 0:max(Int(degree(N)), 0)])
        r = mulmod(Np, invmod(Dp, Ψp), Ψp)
        v = ZZRingElem[lift(ZZ, coeff(r, i)) for i in 0:d-1]
        if isempty(acc)
            acc = v; M = ZZ(pp)
        else
            acc = ZZRingElem[crt(acc[i], M, v[i], ZZ(pp)) for i in 1:d]
            M *= ZZ(pp)
        end
        out = try
            S([QQ(Nemo.reconstruct(acc[i], M)) for i in 1:d]) * scal
        catch
            nothing
        end
        out === nothing && continue
        _redq(out * den, Ψ) == num && return out
    end
    return nothing
end

"Inverse modulo the irreducible Ψ, by the extended Euclidean algorithm."
function _invmodq(f, Ψ)
    g, s, _ = gcdx(_redq(f, Ψ), Ψ)
    iszero(g) && throw(DivideError())
    return _redq(s * inv(coeff(g, 0)), Ψ)
end

# ---------------------------------------------------------------------------------
#  The value
# ---------------------------------------------------------------------------------

"""
    ExactX

An exact value at level `k`, written in the real variable `x = q + q⁻¹ = 2cos(π/h)`, `h = k + 2`:

    v = P(x) · √(R(x)),     P, R ∈ ℚ[x] reduced modulo Ψ_h.

`sqclass` lists the ψ indices under the root — an 𝔽₂ square class over the ψ basis of `generic_x.jl`,
empty when the value needs no root at all, in which case `v = P(x)` outright. It is normalised at
construction: a class whose product collapses to a rational square at this level is folded into `P` and
cleared. Degrees are below `φ(2h)/2`, half the degree of the cyclotomic field the same value occupies in
`Exact(k; form = :canonical)`.

**Printing shows `P(x)·√(R(x))`**, the form the value is stored in, which costs nothing to produce and
is available at every level. [`radical`](@ref) rewrites it in nested square roots on request, and says
why there is none when there is none. Two properties reach past the display:

| | |
|---|---|
| `v.x_value` | the pair `(P, R)`, the value being `P(x)·√(R(x))` |
| `v.rad` | shorthand for `radical(v)` |

[`radical_form`](@ref) is [`radical`](@ref) with a length budget and `nothing` in place of the
explanation, [`xpolynomial`](@ref) and [`radicand`](@ref) give `P` and `R` as coefficient vectors, and
[`has_radical_form`](@ref) answers the existence question without computing the expression.
"""
struct ExactX
    k::Int
    p::QQPolyRingElem          # P, reduced mod Ψ_h
    sqclass::Vector{Int}       # ψ indices under the root
    r::QQPolyRingElem          # R = ∏_{e ∈ sqclass} ψ_e, reduced mod Ψ_h

    function ExactX(k::Integer, p::QQPolyRingElem, cls::Vector{Int}, r::QQPolyRingElem)
        kk = Int(k)
        S, _ = _qqx()
        # A square class is symbolic over the ψ basis of ℚ(x); *at a level* the product can collapse to a
        # rational square, and then there is no root left to carry. Measured on 885 exact 6j values to
        # k = 16, **63** had a non-empty class whose radicand reduced to a constant — every one of them
        # printed `· √(1)`, compared unequal to the same number with an empty class, and propagated the
        # phantom class through every product it entered. Normalising here is the only place that cannot
        # be bypassed.
        iszero(p) && return new(kk, p, Int[], S(1))
        if !isempty(cls) && degree(r) <= 0
            s = _rat_sqrt(Rational{BigInt}(coeff(r, 0)))
            s === nothing || return new(kk, p * QQ(numerator(s), denominator(s)), Int[], S(1))
        end
        return new(kk, p, cls, r)
    end
end

"The exact rational square root of `c`, or `nothing` when `c` is negative or not a square."
function _rat_sqrt(c::Rational{BigInt})
    c < 0 && return nothing
    iszero(c) && return c
    n, d = numerator(c), denominator(c)
    sn, sd = isqrt(n), isqrt(d)
    (sn * sn == n && sd * sd == d) || return nothing
    return Rational{BigInt}(sn, sd)
end

"""
Properties beyond the stored fields.

`v.x_value` is the whole stored form as the pair `(P, R)`, so that `P, R = v.x_value` and
`v.x_value.P`, `v.x_value.R` both work: the value is `P(x)·√(R(x))` at `x = 2cos(π/(k+2))`, with `R = 1`
exactly when `v.sqclass` is empty. It used to be `P` alone, which silently dropped the root from every
value that had one.

`v.rad` is shorthand for [`radical(v)`](@ref radical): the value in **nested square roots**, computed on
demand and with no length limit, or a [`NoRadical`](@ref) saying why there is none.
"""
function Base.getproperty(v::ExactX, s::Symbol)
    s === :rad && return radical(v)
    s === :x_value && return (P = getfield(v, :p), R = getfield(v, :r))
    return getfield(v, s)
end
Base.propertynames(::ExactX, private::Bool = false) =
    private ? (:k, :p, :sqclass, :r, :rad, :x_value) : (:k, :x_value, :rad)

level(v::ExactX) = v.k
Base.iszero(v::ExactX) = iszero(v.p)

"The exact zero at a level — what an inadmissible label set evaluates to, by the package's convention."
function Base.zero(::Type{ExactX}, k::Integer)
    S, _ = _qqx()
    return ExactX(Int(k), S(0), Int[], S(1))
end
Base.zero(v::ExactX) = zero(ExactX, v.k)

"""
    xpolynomial(v::ExactX) -> Vector{Rational{BigInt}}

Coefficients of `P` in `1, x, x², …`, so that `v = P(x)·√(R(x))` with `x = 2cos(π/(k+2))`.
"""
xpolynomial(v::ExactX) = _ratvec(v.p)

"""
    radicand(v::ExactX) -> Vector{Rational{BigInt}}

Coefficients of `R`, the polynomial under the square root. Empty radical class gives `R = 1`.
"""
radicand(v::ExactX) = _ratvec(v.r)

_ratvec(f) = degree(f) < 0 ? Rational{BigInt}[] :
             [Rational{BigInt}(coeff(f, i)) for i in 0:degree(f)]

"""
    exact_x(s::FactorialSum, k) -> ExactX

The exact value of a factorial rule at level `k`, in the real basis. The Racah sum is formed over ℤ[x] by
the same cleared ratio recursion as `generic_value`, reduced modulo `Ψ_h` at every step so the
degree never exceeds `φ(2h)/2`, and the numerator and denominator are then divided in the field.
"""
function exact_x(s::FactorialSum, k::Integer)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
    h = kk + 2
    Ψ = _psiq(h)
    S, _ = _qqx()
    # Reduce as the sum is formed. Cancellation in ℤ[x] can hide behind the reduction — a numerator and
    # denominator sharing a factor that vanishes at this level reduce to 0/0 — so a vanishing denominator
    # is not a verdict: fall back to the unreduced generic value, where the factor cancels exactly, and
    # only then call the rule singular.
    v = nothing
    try
        v = generic_value(s; modulus = psi_level(h), level = h)
        iszero(v.den) && (v = nothing)
    catch e
        e isa DivideError || rethrow()
        v = nothing
    end
    reduced = v !== nothing
    if !reduced
        v = generic_value(s)
    end
    iszero(v.num) && return ExactX(kk, S(0), Int[], S(1))
    den = _redq(_toqq(v.den), Ψ)
    iszero(den) && throw(DomainError(kk, "the rule is singular at level $kk: its denominator vanishes"))
    nred = _redq(_toqq(v.num), Ψ)
    p = _divmod_psi(nred, den, h)
    p === nothing && (p = _redq(nred * _invmodq(den, Ψ), Ψ))
    iszero(p) && return ExactX(kk, S(0), Int[], S(1))
    rp = S(1)
    for e in v.rad                                   # `v` is an `XValue` here: its class field is `rad`
        rp = _redq(rp * _toqq(psi_x(e)), Ψ)
    end
    iszero(rp) && throw(DomainError(kk, "the radical of the rule vanishes at level $kk"))
    return ExactX(kk, p, copy(v.rad), rp)
end

# ---------------------------------------------------------------------------------
#  Numbers
# ---------------------------------------------------------------------------------

"x = 2cos(π/h) at the requested precision."
_x0(h::Int, ::Type{T}) where {T} = 2 * cos(T(pi) / T(h))

function _evalq(f, x0::T) where {T}
    degree(f) < 0 && return zero(T)
    acc = zero(T)
    for i in degree(f):-1:0
        acc = acc * x0 + T(Rational{BigInt}(coeff(f, i)))
    end
    return acc
end

"""
Horner, carrying the running mass `Σ|cᵢ||x|ⁱ` alongside the value. The mass bounds the relative error of
the result — `(2d+1)·u·mass/|value|` — which is the same shape of certificate the numeric kernels use, and
here it is essential rather than decorative: an exact value in `x` is evaluated near `x = 2` where the
polynomial cancels, and the cancellation grows with the level.
"""
function _horner_mass(f, x0::T) where {T}
    acc = zero(T); mass = zero(T); ax = abs(x0)
    for i in degree(f):-1:0
        c = T(Rational{BigInt}(coeff(f, i)))
        acc = acc * x0 + c
        mass = mass * ax + abs(c)
    end
    return acc, mass
end

"Precision the certified evaluation will climb to before giving up."
const EVAL_MAX_BITS = 1 << 17

"""
And the much lower cap the Lagrange descent uses. It needs a *sign*, not a value, so `rtol = 1e-8` is
generous; and it asks for one at every node of a recursion that runs thousands of times inside a sweep, so
letting it climb to `EVAL_MAX_BITS` on a badly cancelling intermediate would cost seconds for an answer
worth a bit. A node that cannot be signed at 4096 bits abandons the descent, which reports no radical form
— the same conservative outcome as a field that has none.
"""
const DESCENT_MAX_BITS = 4096

"""
Largest degree `radical_form` will descend on when it has a length budget.

The descent produces one rational leaf per degree, and the leaves grow as it squares its intermediates,
so past a handful of them the expression cannot fit in `maxlen` however it is written — and running the
descent to find that out is pure cost. Degree 8 is where the line sits in practice: the deepest form the
package has ever displayed is the four-radical nest at k = 14, which is degree 8 and 44 characters, while
at k = 100 a degree-16 value took **8.5 s** to build and render and was then thrown away for being far
past 80. `radical_form(v; maxlen = 0)` ignores the cap and computes it anyway.
"""
const RADICAL_MAX_DEGREE = 8

"""
Largest field degree at which [`radical_levels`](@ref) examines a level the sufficient condition rejects.
Deciding a value's own degree is cheap per level and ruinous per sweep — 11.7 s to `kmax = 400` without a
bound — and what lies past it is, in every case measured, a rational value, which is recognised for free.
"""
const REFINE_MAX_DEGREE = 64

"""
Largest degree of the *value* that [`radical`](@ref) will run the descent on unasked.

The descent is exponential in that degree — one rational leaf per degree, and the leaves grow as it
squares its intermediates. Measured: degree 8 is microseconds, degree 32 is `0.13` s at k = 94, and
`q6j(Exact(254), 45, 45, 45, 30, 30, 30)` — field degree 128, *value* degree 64 — exhausted memory and
took the session down with it. A default that can do that is not a default. Past the cap `radical`
returns a [`NoRadical`](@ref) naming the degree, and `radical(v; degree_limit = …)` is the way to say
that the wait is wanted.
"""
const RADICAL_DESCENT_MAX_DEGREE = 32


"""
Largest field degree at which the *exact* degree of a value is computed.

Deciding whether a radical exists is cheap at any size (a modular conjugation), but the degree itself is
a minimal polynomial over ℚ: 1 ms at field degree 32, 12 ms at 64, 168 ms at 128 measured. Past this
the question is declined rather than paid for, since a value in a field that large is past the descent
budget anyway and the exact number would only be for the sentence.
"""
const RADICAL_MINPOLY_MAX_DEGREE = 64

"""
    _eval_at_x(f, h; rtol) -> BigFloat or nothing

`f(2cos(π/h))`, at whatever precision it takes to certify `rtol`, or `nothing` when `EVAL_MAX_BITS` is not
enough. The precision is not guessed twice: the first pass measures the actual cancellation and the second
is made wide enough for it.

This exists because a fixed precision is wrong. At `k = 420` the value of `{1 1 1; 1 1 1}` is a degree-208
polynomial in `x = 2cos(π/422) = 1.99989…`, where 256 bits print **−64.011** for a symbol whose value is
**0.16666**. The expression was exact; only the number under it was not. An exact value that cannot say
what its number is should say nothing, not something.
"""
function _eval_at_x(f, h::Int; rtol::Real = 1e-16, maxbits::Int = EVAL_MAX_BITS,
                    minbits::Int = 256)
    degree(f) < 0 && return BigFloat(0)
    bits = min(maxbits, max(32, minbits))
    while true
        acc, need = setprecision(BigFloat, bits) do
            a, m = _horner_mass(f, _x0(h, BigFloat))
            d = max(degree(f), 1)
            if iszero(a)
                return a, bits * 2
            end
            rel = (2d + 1) * ldexp(BigFloat(1), -bits) * m / abs(a)
            rel <= rtol && return a, 0
            # how many more bits the measured cancellation asks for
            return a, bits + max(1, Int(ceil(Float64(log2(rel / rtol))))) + 16
        end
        need == 0 && return acc
        bits >= maxbits && return nothing
        bits = min(max(need, bits + 1), maxbits)
    end
end

"""
    float(v::ExactX, T = Float64)

The number, evaluated by [`_eval_at_x`](@ref) at whatever precision certifies it and rounded once. Throws
when no reachable precision certifies it — see [`numeric_value`](@ref) for the form that returns `nothing`
instead, which is what the display uses.
"""
function Base.float(v::ExactX, ::Type{T} = Float64) where {T<:AbstractFloat}
    r = numeric_value(v; bits=precision(T))
    r === nothing && throw(DomainError(v.k,
        "the exact value at level $(v.k) cancels beyond $(EVAL_MAX_BITS) bits when evaluated at " *
        "x = 2cos(π/$(v.k + 2)); the expression is still exact. For a number use the certified numeric " *
        "path, `q6j(labels...; k = $(v.k))`."))
    return T(r)
end
Base.Float64(v::ExactX) = float(v, Float64)
Base.BigFloat(v::ExactX) = float(v, BigFloat)

"""
    evaluate_exact(v, T = ComplexF64)

The number an exact value denotes, in `T`. Named for the cyclotomic carrier this replaces so that code
asking an exact value for its number does not have to know which carrier produced it — the real basis
answers to the same call, and `real(evaluate_exact(v))` keeps meaning what it meant.
"""
evaluate_exact(v::ExactX, ::Type{T} = ComplexF64) where {T} =
    T(float(v, typeof(real(zero(T)))))

"""
    numeric_value(v::ExactX; bits=precision(BigFloat)) -> BigFloat or nothing

Evaluate the selected real embedding with adaptive guard precision, then round to `bits` bits.
Returns `nothing` if the working-precision budget is exhausted. The Horner mass is an error estimate,
not an outward-rounded interval certificate. Negative radicands are rejected rather than clamped to zero.
"""
function numeric_value(v::ExactX; bits::Integer=precision(BigFloat))
    0 < bits <= EVAL_MAX_BITS || throw(ArgumentError("bits must be between 1 and $EVAL_MAX_BITS"))
    target=Int(bits)
    iszero(v.p) && return setprecision(()->BigFloat(0),BigFloat,target)
    h = v.k + 2
    tol=setprecision(()->ldexp(BigFloat(1),-target-8),BigFloat,target+16)
    p = _eval_at_x(v.p, h;rtol=tol,minbits=max(256,target+32))
    p === nothing && return nothing
    isempty(v.sqclass) && return setprecision(()->BigFloat(p),BigFloat,target)
    r = _eval_at_x(v.r, h;rtol=tol,minbits=max(256,target+32))
    r === nothing && return nothing
    r < 0 && throw(DomainError(v.k,"the x-form radicand is negative at this real embedding"))
    result=setprecision(() -> p * sqrt(r), BigFloat, max(precision(p), precision(r)))
    return setprecision(()->BigFloat(result),BigFloat,target)
end

# ---------------------------------------------------------------------------------
#  Radical expressions
# ---------------------------------------------------------------------------------

"""
    RadExpr

A real number written as `rat + Σᵢ cᵢ √(eᵢ)`, with each `eᵢ` again a `RadExpr`. This is what the Lagrange
descent produces and what [`radical_form`](@ref) returns; `float` evaluates it, so a rendered formula can
be checked against the value it claims to be.
"""
struct RadExpr
    rat::Rational{BigInt}
    terms::Vector{Tuple{Rational{BigInt},RadExpr}}
end
RadExpr(r::Rational{BigInt}) = RadExpr(r, Tuple{Rational{BigInt},RadExpr}[])
RadExpr(r::Integer) = RadExpr(Rational{BigInt}(r))

is_rational(e::RadExpr) = isempty(e.terms)

# Structural equality, which `simplify` needs to merge two terms with the same radicand: the default for a
# struct holding a vector is identity, so without this the merge silently never fired.
Base.:(==)(a::RadExpr, b::RadExpr) = a.rat == b.rat && a.terms == b.terms
Base.hash(e::RadExpr, h::UInt) = hash(e.rat, hash(e.terms, h))

function Base.float(e::RadExpr, ::Type{T} = Float64) where {T<:AbstractFloat}
    # A denested radical can cancel savagely — `√((20905 − 9349√5)/40)` loses ten digits, because
    # 9349√5 = 20905.00004… — and how badly depends on the expression, not on a constant. Evaluate at
    # doubling precision until two agree, then round once.
    if T === Float64
        prev = setprecision(() -> _rad_float(e), BigFloat, 256)
        bits = 512
        while bits <= EVAL_MAX_BITS
            cur = setprecision(() -> _rad_float(e), BigFloat, bits)
            isapprox(cur, prev; rtol = 1e-25, atol = 0) && return T(cur)
            prev = cur; bits *= 2
        end
        return T(prev)
    end
    return T(_rad_float(e))
end

function _rad_float(e::RadExpr)
    acc = BigFloat(e.rat)
    for (c, inner) in e.terms
        acc += BigFloat(c) * sqrt(max(zero(BigFloat), _rad_float(inner)))
    end
    return acc
end
Base.Float64(e::RadExpr) = float(e, Float64)

"Nesting depth: 0 for a rational, 1 for a plain surd."
depth(e::RadExpr) = isempty(e.terms) ? 0 : 1 + maximum(depth(t[2]) for t in e.terms)
nleaves(e::RadExpr) = isempty(e.terms) ? 1 : 1 + sum(nleaves(t[2]) for t in e.terms)

# ---------------------------------------------------------------------------------
#  Simplification: rational square parts, and the one denesting that matters
# ---------------------------------------------------------------------------------

_issq(r::Rational{BigInt}) = r >= 0 && isqrt(numerator(r))^2 == numerator(r) &&
                             isqrt(denominator(r))^2 == denominator(r)
_rsqrt(r::Rational{BigInt}) = Rational{BigInt}(isqrt(numerator(r)), isqrt(denominator(r)))

"""
How far `_split_surd` trial-divides before leaving the rest under the root.

Without a bound this is a hang, not a slowdown. The Lagrange descent squares its intermediates, so the
leaf integers grow with the level: 29 digits at k = 94, **41** at k = 100 for the same labels — and
trial division runs to `√n`, i.e. 10¹⁴ steps against 10²⁰. That is exactly the difference the levels
showed, one printing "too long to be useful" in a couple of seconds and the other never returning.

The bound costs nothing mathematically. A factor left inside the root makes the surd less tidy, never
wrong, and the one case worth catching past small primes — the whole cofactor being a perfect square — is
one integer square root away.
"""
const SURD_TRIAL_LIMIT = 10_000

"Primes below `SURD_TRIAL_LIMIT`, sieved once: trial dividing by composites is four fifths wasted work."
const _SURD_PRIMES = let n = SURD_TRIAL_LIMIT
    sieve = trues(n)
    sieve[1] = false
    for i in 2:isqrt(n)
        sieve[i] || continue
        for j in i*i:i:n
            sieve[j] = false
        end
    end
    BigInt[i for i in 2:n if sieve[i]]
end

"`√(a/b)` as `c·√m`: rationalise, pull out the square part of the small primes, then test the rest."
function _split_surd(r::Rational{BigInt})
    r >= 0 || return nothing
    n = numerator(r) * denominator(r)        # √(n/d) = √(nd)/d
    d = denominator(r)
    c = Rational{BigInt}(1, d)
    m = BigInt(1)
    for f in _SURD_PRIMES
        f * f > n && break
        e = 0
        while n % f == 0
            n ÷= f; e += 1
        end
        iszero(e) && continue
        c *= Rational{BigInt}(f)^(e ÷ 2)
        isodd(e) && (m *= f)
    end
    if n > 1
        sq = isqrt(n)
        sq * sq == n ? (c *= Rational{BigInt}(sq)) : (m *= n)
    end
    return c, m
end

"""
Collapse what can be collapsed: a rational radicand becomes `c√m`, a perfect square leaves the root, and
`√(t + c√d)` denests to `√X ± √Y` whenever `t² − c²d` is a rational square. Everything else is left alone.
"""
function simplify(e::RadExpr)
    acc = e.rat
    out = Tuple{Rational{BigInt},RadExpr}[]
    for (c, inner0) in e.terms
        iszero(c) && continue
        inner = simplify(inner0)
        if is_rational(inner)
            inner.rat < 0 && return e            # not a real number in this shape; leave it
            s = _split_surd(inner.rat)
            s === nothing && return e
            cc, m = s
            if m == 1
                acc += c * cc
            else
                push!(out, (c * cc, RadExpr(Rational{BigInt}(m))))
            end
        else
            dn = _denest(inner)
            if dn === nothing
                push!(out, (c, inner))
            else
                for (c2, in2) in dn
                    s = _split_surd(in2)
                    if s === nothing
                        push!(out, (c * c2, RadExpr(in2)))
                    else
                        cc, m = s
                        m == 1 ? (acc += c * c2 * cc) : push!(out, (c * c2 * cc, RadExpr(Rational{BigInt}(m))))
                    end
                end
            end
        end
    end
    # merge equal radicands
    merged = Tuple{Rational{BigInt},RadExpr}[]
    for (c, inner) in out
        i = findfirst(t -> t[2] == inner, merged)
        i === nothing ? push!(merged, (c, inner)) : (merged[i] = (merged[i][1] + c, inner))
    end
    filter!(t -> !iszero(t[1]), merged)
    sort!(merged, by = t -> (is_rational(t[2]) ? float(t[2]) : Inf))
    return RadExpr(acc, merged)
end

"√(t + c√d) = √X ± √Y with X + Y = t and 4XY = c²d, when t² − c²d is a rational square."
function _denest(inner::RadExpr)
    length(inner.terms) == 1 || return nothing
    c, sub = inner.terms[1]
    is_rational(sub) || return nothing
    t = inner.rat; d = sub.rat
    disc = t * t - c * c * d
    (disc < 0 || !_issq(disc)) && return nothing
    g = _rsqrt(disc)
    X = (t + g) / 2; Y = (t - g) / 2
    (X < 0 || Y < 0) && return nothing
    return [(Rational{BigInt}(1), X), (c < 0 ? Rational{BigInt}(-1) : Rational{BigInt}(1), Y)]
end

# ---------------------------------------------------------------------------------
#  The Lagrange descent
# ---------------------------------------------------------------------------------

"Residues `j` with `gcd(j, 2h) = 1` and `1 ≤ j < h`: one per conjugate `2cos(jπ/h)` of `x`."
_conj_reps(h::Int) = [j for j in 1:max(h - 1, 1) if gcd(j, 2h) == 1]

function _rep(j::Int, h::Int)
    r = mod(j, 2h)
    return r > h ? 2h - r : r
end

const _CHEB_RED = Dict{Tuple{Int,Int},Any}()

"`C_j` reduced modulo `Ψ_h`: the image of `x` under the j-th conjugate, built once per `(j, h)`."
function _cheb_red(j::Int, h::Int)
    lock(X_LOCK) do
        get!(_CHEB_RED, (h, j)) do
            _redq(_toqq(cheb_x(j)), _psiq(h))
        end
    end
end

"σ_j applied to `f`: substitute x ↦ C_j(x), the j-th conjugate, and reduce."
_sigma(f, j::Int, h::Int, Ψ) = j == 1 ? f : _redq(evaluate(f, _cheb_red(j, h)), Ψ)

const _ODD_GENS = Dict{Int,Vector{Int}}()

"`a^e` as a conjugate index, folded back into `1 ≤ · ≤ h`."
_powrep(a::Int, e::Int, h::Int) = _rep(Int(powermod(a, e, 2h)), h)

"""
A generating set for the odd part of the Galois group `G = (ℤ/2h)*/±1`, empty when `G` is a 2-group.

`G` is abelian, so `G ≅ G₂ × G_odd` and raising to the 2-part of `|G|` kills `G₂` and permutes `G_odd`:
the image of that map *is* `G_odd`. A generating set of it is then a couple of elements, found by
closure over the integers — no field arithmetic anywhere in here.
"""
function _odd_gens(h::Int)
    lock(X_LOCK) do
        get!(_ODD_GENS, h) do
            reps = _conj_reps(h)
            d = length(reps)
            m = trailing_zeros(d)
            (d >> m) == 1 && return Int[]
            odd = unique!([_powrep(a, 1 << m, h) for a in reps])
            gens = Int[]
            cur = Set{Int}(1)
            for g in odd
                g in cur && continue
                push!(gens, g)
                cur = _close_group(gens, h)
                length(cur) >= length(odd) && break
            end
            return gens
        end
    end
end

"""
Is `[ℚ(u):ℚ]` a power of two — the whole radical question — without computing the degree?

`[G : Stab(u)]` is a power of two exactly when `G_odd ⊆ Stab(u)`: one direction because the index then
divides `|G/G_odd|`, the other because an odd-order group has no nontrivial image in a 2-group. So the
answer is a couple of conjugations, against a minimal polynomial of degree `φ(2h)/2`. Measured at
k = 420: **12 ms against 164 ms**, and the gap grows with the level.
"""
function _degree_is_2power(u, h::Int)
    degree(u) <= 0 && return true
    gens = _odd_gens(h)
    isempty(gens) && return true                 # a 2-group has only 2-power indices
    Ψ = _psiq(h)
    # A difference seen modulo one word-sized prime is a difference, full stop, and that is the answer
    # in the overwhelming majority of cases — 0.9 ms against 121 ms at k = 420 measured. Agreement
    # modulo a prime proves nothing, so the exact conjugation is still run when no prime separates them.
    for g in gens
        _sigma_differs_mod(u, g, h, Ψ) && return false
    end
    for g in gens
        _sigma(u, g, h, Ψ) == u || return false
    end
    return true
end

"""
Do `σ_j(u)` and `u` differ? `true` is a proof; `false` means only that this prime saw no difference.

Composition modulo `Ψ_h` over `𝔽_p` is machine-word arithmetic where the exact route carries the
value's own coefficients — 144 bits at k = 420 — through 208 polynomial multiplications.
"""
function _sigma_differs_mod(u, j::Int, h::Int, Ψ)
    j == 1 && return false
    U, _ = _clear_denoms(u)                      # the common denominator cancels: σ is ℚ-linear
    C, aC = _clear_denoms(_cheb_red(j, h))
    Z, _ = _clear_denoms(Ψ)
    du = Int(degree(U))
    du < 0 && return false
    for pp in _MM_PRIMES[1:2]
        F = Nemo.Native.GF(pp)
        Fx, _ = polynomial_ring(F, "X")
        iszero(F(aC)) && continue
        Ψp = Fx([F(coeff(Z, i)) for i in 0:Int(degree(Z))])
        Cp = Fx([F(coeff(C, i)) for i in 0:max(Int(degree(C)), 0)]) * inv(F(aC))
        Up = Fx([F(coeff(U, i)) for i in 0:du])
        acc = zero(Fx)
        for i in du:-1:0
            acc = mulmod(acc, Cp, Ψp) + Fx(F(coeff(U, i)))
        end
        acc == Up || return true
    end
    return false
end


"""
    has_radical_form(k) -> Bool

Whether every exact value at level `k` admits a closed form in real nested square roots. True exactly when
`φ(2k+4)/2` is a power of two, because the real cyclotomic field is then a tower of quadratic extensions;
false when the degree carries an odd prime factor, where no real radical expression exists at all.

Equivalently, and more memorably: **the levels with radicals are the ones where the regular `(k+2)`-gon is
constructible with ruler and compass** — `k+2` a power of two times distinct Fermat primes. Up to 30 that
is `k = 0,1,2,3,4,6,8,10,13,14,15,18,22,28,30`; Gauss's 17-gon is `k = 15`.

A *particular* symbol can still be a radical at a level this rejects, if its value falls into a 2-power
subfield — see [`has_radical_form(::ExactX)`](@ref) and [`radical_levels`](@ref).

```julia
has_radical_form.(0:15)   # false only at k = 5, 7, 9, 11, 12
```
"""
function has_radical_form(k::Integer)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "level must be nonnegative"))
    d = euler_phi(2 * (kk + 2)) ÷ 2
    return d >= 1 && count_ones(d) == 1
end

const _NFIELD = Dict{Int,Any}()

"The number field `ℚ(2cos(π/h))` as Nemo knows it, cached; used only for minimal polynomials."
function _nfield(h::Int)
    lock(X_LOCK) do
        get!(_NFIELD, h) do
            number_field(_psiq(h), "a")
        end
    end
end

"""
    _value_degree(u, h) -> Int

Degree of a field element of `ℚ(2cos(π/h))` over ℚ, from its minimal polynomial.

The obvious route — count the conjugates that fix it, orbit–stabiliser — is one polynomial composition
modulo `Ψ_h` per conjugate, and that is quadratic in a degree that is itself the field's: **6 s** at
k = 420, where the field has degree 210, against **0.12 s** for Nemo's `minpoly` in the number field. The
predicate is called from `show`, so the difference is the difference between a display and a hang.
"""
function _value_degree(u, h::Int)
    degree(u) <= 0 && return 1
    K, _ = _nfield(h)
    return Int(degree(minpoly(K(u))))
end

"""
    has_radical_form(v::ExactX) -> Bool

Whether **this** value has a closed form in real nested square roots — a weaker question than
[`has_radical_form(k)`](@ref), which asks it of every value at the level.

`v²` lies in the real cyclotomic field, so `ℚ(v²)` is abelian whatever the level, and an abelian field is
a tower of quadratic extensions exactly when its degree is a power of two. The value can therefore land in
a 2-power *subfield* of a field that has none: at `k = 5` the whole field has degree 3, but a symbol whose
value happens to be rational there is still a radical expression.

```julia
has_radical_form(q6j(Exact(5; form = :x), 1, 1, 1, 1, 1, 1))
```
"""
function has_radical_form(v::ExactX)
    iszero(v.p) && return true
    has_radical_form(v.k) && return true
    return _degree_is_2power(_square(v), v.k + 2)
end

"`v²` as an element of the field: it is real and abelian whatever the level, which is what the radical
criterion rests on."
function _square(v::ExactX)
    Ψ = _psiq(v.k + 2)
    return _redq(_redq(v.p * v.p, Ψ) * v.r, Ψ)
end

"""
    radical_levels(f, labels...; kmax = 64) -> Vector{Int}

The levels `k ≤ kmax` at which this symbol is admissible **and** its exact value can be written in real
nested square roots.

Two things put a level in the list. Most come for free: when `φ(2k+4)/2` is a power of two *every* value
at that level is a radical expression, and no symbol needs to be computed. The rest are levels where the
field itself has an odd prime in its degree but this particular value falls into a 2-power subfield, and
those are found by computing the value.

Inadmissible levels are left out; a level where the symbol vanishes is included, since `0` is as closed a
form as there is. This answers "for which k can I see this symbol in radicals?", which is not the same
question as [`has_radical_form(k)`](@ref).

**Cost.** The second kind of level is the expensive one, and it gets more expensive with `k`: an unbounded
sweep took 0.56 s to `kmax = 200`, 3.3 s to 300 and 11.7 s to 400. So the refinement is only attempted
while the field degree is at most `refine_max_degree` (default `REFINE_MAX_DEGREE`), and
`refine = false` turns it off entirely, leaving the sufficient condition — instant, and every level it
lists is certain. What a bounded sweep can miss is a level where the value lies in a 2-power *subfield*
of a field that has none; measured across k = 5, 7, 9, 11, 12, 16, 17, every such value was **rational**,
and rational values are recognised for free at any degree.

```julia
radical_levels(q6j, 1, 1, 1, 1, 1, 1; kmax = 24)
radical_levels(q6j, 1, 1, 1, 1, 1, 1; kmax = 500, refine = false)   # instant, sufficient only
```
"""
function radical_levels(f, labels...; kmax::Integer = 64, refine::Bool = true,
                        refine_max_degree::Integer = REFINE_MAX_DEGREE)
    kk = Int(kmax)
    kk >= 0 || throw(DomainError(kmax, "kmax must be nonnegative"))
    sym = symbol_of(f)
    sym === nothing && throw(ArgumentError(
        "radical_levels needs a function attached to a symbol rule: q6j, q3j, fsymbol, gsymbol"))
    out = Int[]
    for k in 0:kk
        level_admissible(sym, k, labels...) || continue
        if has_radical_form(k)                       # free: the level answers for every value
            push!(out, k)
            continue
        end
        refine || continue
        euler_phi(2 * (k + 2)) ÷ 2 <= refine_max_degree || continue
        v = exact_x(symbol_rule(sym, labels...), k)
        # A rational value is a radical form and costs nothing to recognise; measured, it is also the
        # *only* way a value has ever beaten the level's criterion, so this shortcut is the usual path.
        if iszero(v.p) || (degree(v.p) <= 0 && isempty(v.sqclass)) || has_radical_form(v)
            push!(out, k)
        end
    end
    return out
end

"The subgroup of `G = (ℤ/2h)*/±1` generated by `gens`, as a set of representatives."
function _close_group(gens, h::Int)
    S = Set{Int}(1)
    frontier = Int[1]
    while !isempty(frontier)
        a = pop!(frontier)
        for g in gens
            c = _rep(a * g, h)
            c in S && continue
            push!(S, c); push!(frontier, c)
        end
    end
    return S
end

"""
One representative per coset of the subgroup `fixed`. Everything in a coset acts the same way on an
element `fixed` already fixes — `G` is abelian, so `σ_{af}(u) = σ_a(σ_f(u)) = σ_a(u)` — so conjugating by
more than one of them is pure repetition.
"""
function _coset_reps(reps::Vector{Int}, fixed::Set{Int}, h::Int)
    length(fixed) <= 1 && return reps
    seen = Set{Int}()
    out = Int[]
    for a in reps
        a in seen && continue
        push!(out, a)
        for b in fixed
            push!(seen, _rep(a * b, h))
        end
    end
    return out
end

"""
Write a field element as a `RadExpr`, or `nothing` when the field is not a 2-group tower, or when a
sign could not be certified. `u` is reduced modulo `Ψ_h`.
"""
function _descend(u, h::Int, Ψ, reps::Vector{Int}, fixed::Set{Int} = Set{Int}(1))
    degree(u) <= 0 && return RadExpr(Rational{BigInt}(degree(u) < 0 ? 0 : coeff(u, 0)))
    # Conjugation is 94% of this recursion — 112 calls and 4.18 ms of 4.46 at j = 8, k = 30 — and most of
    # them say nothing. The group fixing the value only *grows* as the descent goes down, so a node need
    # only conjugate by one representative per coset of what its parent already fixes: 16 at the root,
    # then 8, 4, 2. The images are then reused three times over, for the stabiliser, for the search for
    # σ, and for σ(u) itself.
    cosets = _coset_reps(reps, fixed, h)
    imgs = [_sigma(u, a, h, Ψ) for a in cosets]
    stab = _close_group(vcat(collect(fixed),
                             [cosets[i] for i in eachindex(cosets) if imgs[i] == u]), h)
    j = nothing; w = nothing
    for (i, a) in enumerate(cosets)
        a in stab && continue
        _rep(a * a, h) in stab || continue
        j = a; w = imgs[i]; break
    end
    j === nothing && return nothing
    half = QQ(1, 2)
    t = _redq((u + w) * half, Ψ)
    r = _redq((u - w) * half, Ψ)
    s = _redq(r * r, Ψ)
    # both halves lie in the fixed field of ⟨Stab(u), σ⟩, which is what the children may skip over
    below = _close_group(vcat(collect(stab), [j]), h)
    te = _descend(t, h, Ψ, reps, below)
    te === nothing && return nothing
    se = _descend(s, h, Ψ, reps, below)
    se === nothing && return nothing
    # `r` is nonzero (σ moves u), but it can be small: take the sign only from a certified evaluation.
    rv = _eval_at_x(r, h; rtol = 1e-8, maxbits = DESCENT_MAX_BITS)
    rv === nothing && return nothing
    sgn = rv < 0 ? Rational{BigInt}(-1) : Rational{BigInt}(1)
    return RadExpr(te.rat, vcat(te.terms, [(sgn, se)]))
end

"""
    radical_form(v::ExactX; maxlen = 80) -> RadExpr or nothing

The exact value as real nested square roots, or `nothing` when there is none to give: either the level
fails [`has_radical_form`](@ref) — in which case no such expression exists, for any amount of effort — or
the expression exists but is longer than `maxlen` characters, which happens as soon as the descent needs
three levels, or the value's degree is past `degree_limit`, where the descent stops being affordable
(`degree_limit = 0` removes that cap and accepts the cost). [`radical`](@ref) is the same thing with an
unbounded length, a sentence in place of the `nothing`, and the number under it; this is the form to
call when the caller wants to branch on the answer rather than read it.

```julia
radical_form(q6j(Exact(3; form = :x), 1, 1, 1, 1, 1, 1))   # −(3 − √5)/2, the Fibonacci level
```
"""
function radical_form(v::ExactX; maxlen::Int = 80,
                     degree_limit::Int = RADICAL_DESCENT_MAX_DEGREE)
    iszero(v.p) && return RadExpr(0)
    if degree(v.p) <= 0 && (isempty(v.sqclass) || degree(v.r) <= 0)
        c = Rational{BigInt}(coeff(v.p, 0))
        isempty(v.sqclass) && return RadExpr(c)
        r = Rational{BigInt}(coeff(v.r, 0))
        r < 0 && return nothing
        e = simplify(RadExpr(big(0)//big(1), [(c, RadExpr(r))]))
        return maxlen > 0 && length(_rad_str(e)) > maxlen ? nothing : e
    end
    # The *value*, not the level: a symbol can be rational — or quadratic — at a level whose field has an
    # odd prime in its degree, and refusing on the level alone hid exactly those. Measured, 12 of 1632
    # sampled values were rational at k = 5, 7, 9 and were being told no radical form existed.
    has_radical_form(v) || return nothing
    h = v.k + 2
    # Two caps, one minimal polynomial. A degree-2ᵐ value descends to 2ᵐ rational leaves, so with a
    # length budget anything past a handful of them is already far longer than `maxlen` and the descent
    # would only be paying to be thrown away. Without one the cap is the cost itself: the descent is
    # exponential in that degree, and at 64 it exhausts memory. `degree_limit = 0` says the caller has
    # already decided the degree, or accepts whatever it costs.
    cap = maxlen > 0 ? (degree_limit > 0 ? min(degree_limit, RADICAL_MAX_DEGREE) : RADICAL_MAX_DEGREE) :
                       degree_limit
    if cap > 0
        euler_phi(2h) ÷ 2 <= RADICAL_MINPOLY_MAX_DEGREE || return nothing
        _value_degree(_square(v), h) > cap && return nothing
    end
    Ψ = _psiq(h)
    reps = _conj_reps(h)
    e = nothing
    if isempty(v.sqclass)
        e = _descend(v.p, h, Ψ, reps)                    # no root: descend on the value itself
    end
    if e === nothing
        w = _redq(_redq(v.p * v.p, Ψ) * v.r, Ψ)          # v² is in the real field; v = ±√(v²)
        inner = _descend(w, h, Ψ, reps)
        sg = inner === nothing ? nothing : numeric_value(v)
        e = (inner === nothing || sg === nothing) ? nothing :
            RadExpr(Rational{BigInt}(0),
                    [(sg < 0 ? Rational{BigInt}(-1) : Rational{BigInt}(1), inner)])
    end
    e === nothing && return nothing
    e = simplify(e)
    maxlen > 0 && length(_rad_str(e)) > maxlen && return nothing
    return e
end

"""
    NoRadical

What [`radical`](@ref) gives back when there is no nested-radical expression to give, carrying the reason
in `kind`, so that "there is none" and "I did not look" are different answers and not one sentence doing
duty for both:

* `:none` — none exists. `v²` lies in an abelian field, and an abelian field is a tower of quadratic
  extensions exactly when its degree is a power of two, so this is decided, not merely unattempted.
* `:untried` — the field degree is past `degree_limit` and the minimal polynomial that would settle it
  was not computed.
* `:long` — one exists but is longer than the `maxlen` asked for. Only `radical(v; maxlen = n)` with
  `n > 0` can produce this; the default budget is unbounded.
* `:failed` — the descent ran and could not certify a sign within `DESCENT_MAX_BITS`.

Printing it prints the reason. `float` is deliberately not defined: there is no number here.
"""
struct NoRadical
    k::Int
    kind::Symbol
    degree::Int        # of the value, 0 when it was not worth computing
    field::Int         # φ(2h)/2, always known
    limit::Int
    approx::Union{Nothing,Float64}
end

"""
The sentence, and under it the number. A reader who asked for a closed form and cannot have one is owed
the value anyway — that is the whole reason they were looking — and it is the one thing always available.
"""
function _no_radical_str(n::NoRadical)
    body = if n.kind === :none
        (n.degree > 0 ?
            "no radical form: v² has degree " * string(n.degree) * " over ℚ at level " * string(n.k) :
            "no radical form: the degree of ℚ(v²) at level " * string(n.k) *
            " has an odd prime factor") *
        ", and only a power of two is a tower of square roots\n" *
        "  (`v.x_value` gives (P, R), which is then the value's only closed form)"
    elseif n.kind === :long
        "nested square roots, longer than the budget asked for\n" *
        "  (`radical(v)` returns the expression whatever its size)"
    elseif n.kind === :untried
        "radical form not attempted: " *
        (n.degree > n.limit ?
            "the value at level " * string(n.k) * " has degree " * string(n.degree) *
            ", past degree_limit = " * string(n.limit) * "\n" *
            "  (`radical(v; degree_limit = " * string(n.degree) *
            ")` runs it — the descent is exponential in the degree)" :
         n.field > RADICAL_MINPOLY_MAX_DEGREE ?
            "the field at level " * string(n.k) * " has degree " * string(n.field) *
            ", too large to decide the value's own\n" *
            "  (`radical(v; degree_limit = " * string(n.field) * ")` tries anyway)" :
            "the value at level " * string(n.k) * " has degree past degree_limit = " *
            string(n.limit) * "\n" *
            "  (raise `degree_limit` to descend anyway — the cost doubles with every degree)")
    else
        "no radical expression could be built: the descent could not certify a sign within " *
        string(DESCENT_MAX_BITS) * " bits"
    end
    n.approx === nothing && return body
    return body * "\n  ≈ " * string(n.approx)
end

Base.show(io::IO, n::NoRadical) = print(io, "NoRadical(:", n.kind, ", k = ", n.k, ")")
Base.show(io::IO, ::MIME"text/plain", n::NoRadical) = print(io, _no_radical_str(n))

"The number to print under the sentence, or `nothing` when even that cannot be certified."
function _no_radical(v::ExactX, kind::Symbol, degree::Int, limit::Int)
    nv = numeric_value(v)
    return NoRadical(v.k, kind, degree, euler_phi(2 * (v.k + 2)) ÷ 2, limit,
                     nv === nothing ? nothing : Float64(nv))
end

"""
    radical(v::ExactX; maxlen = 0, degree_limit = 32) -> RadExpr or NoRadical

The value written in real nested square roots — `(√5 − 3)/2` rather than the polynomial in `x` that
printing an [`ExactX`](@ref) shows. Most values have no such form: only a level whose field degree
`φ(2h)/2` is a power of two puts every one of its values in a tower of square roots, and at the other
levels a particular value may still land in a 2-power subfield. When there is none the answer is a
[`NoRadical`](@ref) that says which of those it is, rather than a bare `nothing`.

There is no *length* budget by default: asking for the radical is asking for all of it. There is a
**degree** budget, because the descent is exponential in the degree of the value and at degree 64 it
exhausts memory and takes the session with it. Past `degree_limit` the answer names the degree it found
and the call that would run it, and gives the number meanwhile. `maxlen > 0` additionally declines an
expression longer than that many characters.

```julia
radical(q6j(Exact(3), 1, 1, 1, 1, 1, 1))    # (√5 − 3)/2, the Fibonacci level
radical(q6j(Exact(5), 1, 1, 1, 1, 1, 1))    # no radical form: v² has degree 3 over ℚ at level 5, …
radical(q6j(Exact(254), 45, 45, 45, 30, 30, 30))          # degree 64: not attempted, with the number
radical(q6j(Exact(254), 45, 45, 45, 30, 30, 30); degree_limit = 64)   # …and this is the long wait
```
"""
function radical(v::ExactX; maxlen::Int = 0, degree_limit::Int = RADICAL_DESCENT_MAX_DEGREE)
    kind, e, d = _radical_view(v; maxlen = maxlen, degree_limit = degree_limit)
    kind === :zero && return RadExpr(0)
    kind === :ok && return e
    return _no_radical(v, kind, d, degree_limit)
end

# ---------------------------------------------------------------------------------
#  Rendering
# ---------------------------------------------------------------------------------

_num_str(r::Rational{BigInt}) = denominator(r) == 1 ? string(abs(numerator(r))) :
                                string(abs(numerator(r))) * "/" * string(denominator(r))

"`c√m`, or `c√(…)` for a nested radicand, with a unit coefficient left implicit."
function _term_str(c::Rational{BigInt}, inner::RadExpr)
    body = is_rational(inner) ? "√" * _num_str(inner.rat) : "√(" * _rad_str(inner) * ")"
    a = abs(c)
    isone(a) && return body
    n = numerator(a); d = denominator(a)
    lead = isone(n) ? "" : string(n)
    return isone(d) ? lead * body : lead * body * "/" * string(d)
end

"""
A radical expression over a common denominator, so `3/2 − (1/2)√5` prints as `(3 − √5)/2` rather than as a
sum of fractions, and with the positive parts first, so it reads `(√5 − 3)/2` and not `(−3 + √5)/2`.
"""
function _rad_str(e::RadExpr)
    isempty(e.terms) && return (e.rat < 0 ? "−" : "") * _num_str(e.rat)
    den = denominator(e.rat)
    for (c, _) in e.terms
        den = lcm(den, denominator(c))
    end
    items = Tuple{Bool,String}[]                     # (negative?, body)
    r = e.rat * den
    iszero(r) || push!(items, (r < 0, _num_str(r)))
    for (c, inner) in e.terms
        push!(items, (c < 0, _term_str(c * den, inner)))
    end
    sort!(items, by = t -> t[1])                     # positives first; stable, so order is otherwise kept
    body = (items[1][1] ? "−" : "") * items[1][2]
    for (neg, str) in items[2:end]
        body *= (neg ? " − " : " + ") * str
    end
    den == 1 && return body
    return length(items) == 1 ? body * "/" * string(den) : "(" * body * ")/" * string(den)
end

Base.show(io::IO, e::RadExpr) = print(io, _rad_str(e))

"A polynomial in x over a common denominator: `(2x² − 7)/3`."
function _xpoly_str(f)
    degree(f) < 0 && return "0"
    den = BigInt(1)
    for i in 0:degree(f)
        den = lcm(den, denominator(Rational{BigInt}(coeff(f, i))))
    end
    out = IOBuffer()
    first = true
    for i in degree(f):-1:0
        c = Rational{BigInt}(coeff(f, i)) * den
        iszero(c) && continue
        n = numerator(c)
        if first
            n < 0 && print(out, "−")
            first = false
        else
            print(out, n < 0 ? " − " : " + ")
        end
        a = abs(n)
        (isone(a) && i > 0) || print(out, a)
        i == 0 && continue
        print(out, "x")
        i == 1 || print(out, to_superscript(i))
    end
    body = String(take!(out))
    isempty(body) && (body = "0")
    den == 1 && return body
    return occursin(" ", body) ? "(" * body * ")/" * string(den) : body * "/" * string(den)
end

"""
Longest polynomial the display writes out in full. Past it the line is a size, not a formula: at k = 420
the value is a degree-208 polynomial with 240-bit coefficients, and printing it is three screens of
digits that tell a reader nothing `xpolynomial(v)` would not tell them better.
"""
const XPOLY_MAX_DEGREE = 24
const XPOLY_MAX_CHARS = 240

"A polynomial too large to read, described instead: the same courtesy `phi_form` extends."
function _xpoly_summary(f)
    d = degree(f)
    cs = [Rational{BigInt}(coeff(f, i)) for i in 0:d]
    den = BigInt(1)
    for c in cs
        den = lcm(den, denominator(c))
    end
    bits = maximum(c -> ndigits(numerator(c * den); base = 2), cs; init = 0)
    nz = count(!iszero, cs)
    body = "an integer polynomial in x of degree $d, $nz nonzero coefficients, ≤$bits bits"
    return den == 1 ? "(" * body * ")" : "(" * body * ")/" * string(den)
end

"`_xpoly_str` when it fits the budget, `nothing` when it does not."
function _xpoly_fits(f; maxdeg::Int = XPOLY_MAX_DEGREE, maxchars::Int = XPOLY_MAX_CHARS)
    degree(f) > maxdeg && return nothing
    str = _xpoly_str(f)
    return length(str) > maxchars ? nothing : str
end

"`_xpoly_str`, or a description of it when writing it out would be worse than useless."
function _xpoly_show(f; maxdeg::Int = XPOLY_MAX_DEGREE, maxchars::Int = XPOLY_MAX_CHARS)
    str = _xpoly_fits(f; maxdeg = maxdeg, maxchars = maxchars)
    return str === nothing ? _xpoly_summary(f) : str
end

"`degree 208, 105 terms, ≤144 bits`, the shape of a polynomial without the polynomial."
function _xpoly_dims(f)
    d = degree(f)
    d < 0 && return "0"
    cs = [Rational{BigInt}(coeff(f, i)) for i in 0:d]
    den = BigInt(1)
    for c in cs
        den = lcm(den, denominator(c))
    end
    bits = maximum(c -> ndigits(numerator(c * den); base = 2), cs; init = 0)
    out = "degree $d, $(count(!iszero, cs)) terms, ≤$bits bits"
    return den == 1 ? out : "(" * out * ")/" * string(den)
end

"Longest square class the sizes line writes out; past it the count says as much and fits."
const SQCLASS_MAX_SHOWN = 6

"The class under the root, as ψ indices or — past a handful of them — as how many there are."
function _sqclass_str(cls::Vector{Int})
    isempty(cls) && return "1"
    length(cls) > SQCLASS_MAX_SHOWN && return string(length(cls)) * " ψ factors"
    return join(["ψ" * to_subscript(e) for e in cls], "·")
end

"Is either half of `P(x)·√(R(x))` past what the display will write out?"
function _xform_long(v::ExactX)
    iszero(v.p) && return false
    _xpoly_fits(v.p) === nothing && return true
    return !isempty(v.sqclass) && _xpoly_fits(v.r) === nothing
end

"The sizes of `P` and `R`, one line each, for the lines under a schematic form."
function _xform_dims(v::ExactX)
    isempty(v.sqclass) && return ["P: " * _xpoly_dims(v.p) * ",  R = 1"]
    return ["P: " * _xpoly_dims(v.p),
            "R: " * _xpoly_dims(v.r) * ",  " * _sqclass_str(v.sqclass)]
end

"""
The `P(x)·√(R(x))` line: what printing a value shows, always available and never computing anything —
no minimal polynomial, no Lagrange descent, nothing that can fail or take unbounded time. A radicand that
reduced to a bare constant keeps its class (it is not a rational square, or the constructor would have
folded it) and is written `√2` rather than `√(2)`.
"""
function _xform_str(v::ExactX)
    iszero(v.p) && return "0"
    # Past the budget the honest line is the *shape* of the value, `P(x)·√(R(x))`, with the sizes on the
    # line below: a reader takes nothing from three screens of digits, and `v.x_value` gives both halves
    # in full to anyone who wants them.
    _xform_long(v) && return isempty(v.sqclass) ? "P(x)" : "P(x) · √(R(x))"
    ps = _xpoly_show(v.p)
    isempty(v.sqclass) && return ps
    rs = _xpoly_show(v.r)
    root = occursin(" ", rs) || occursin("x", rs) ? "√(" * rs * ")" : "√" * rs
    ps == "1" && return root
    ps == "-1" && return "−" * root
    need = occursin(" ", ps) && !startswith(ps, "(")
    return (need ? "(" * ps * ")" : ps) * " · " * root
end

function Base.show(io::IO, v::ExactX)
    print(io, "ExactX(k = ", v.k, ", ", _xform_str(v), ")")
end

"The value's degree when that is cheap to know, `0` when it is only worth a sentence."
_cheap_degree(u, h::Int, d::Int) = d <= RADICAL_MINPOLY_MAX_DEGREE ? _value_degree(u, h) : 0

"""
    _radical_view(v; maxlen, degree_limit, allowed) -> (kind, expr, degree)

What [`radical`](@ref) can say about this value, decided in one place so that the answer and the sentence
explaining it cannot disagree. `:zero` and `:ok` carry an expression; `:none`, `:untried`, `:long` and
`:failed` become a [`NoRadical`](@ref) and are documented there.
"""
function _radical_view(v::ExactX; maxlen::Int = 80,
                      degree_limit::Int = RADICAL_DESCENT_MAX_DEGREE, allowed::Bool = true)
    iszero(v.p) && return (:zero, nothing, 0)
    h = v.k + 2
    d = euler_phi(2h) ÷ 2
    # A rational value, or a rational multiple of one surd, is its own radical: no field work, no
    # descent, no cap. This must come first, or a degenerate value at a big level is refused for the
    # size of a field it does not live in.
    if degree(v.p) <= 0 && (isempty(v.sqclass) || degree(v.r) <= 0)
        e = radical_form(v; maxlen = maxlen, degree_limit = 0)
        e === nothing && return (maxlen > 0 ? :long : :failed, nothing, 1)
        return (:ok, e, 1)
    end
    allowed || return (:untried, nothing, 0)
    vd = d
    # Three questions, in increasing cost, and none of them asked unless the answer is needed:
    #   1. does a radical exist?     `G_odd ⊆ Stab(v²)` — a couple of conjugations
    #   2. is the degree affordable? the orbit, stopped as soon as it passes the limit
    #   3. what is the degree?       the whole orbit, and only when it is small enough to be cheap
    # The old route asked (3) always, through a minimal polynomial of the field's degree: 164 ms at
    # k = 420 where (1) settles it in 12 ms.
    if !has_radical_form(v.k) || d > degree_limit
        u = _square(v)
        # Existence first, and cheaply: it is the question, and the degree is only the sentence.
        _degree_is_2power(u, h) || return (:none, nothing, _cheap_degree(u, h, d))
        d <= RADICAL_MINPOLY_MAX_DEGREE || return (:untried, nothing, 0)
        vd = _value_degree(u, h)
        vd <= degree_limit || return (:untried, nothing, vd)
    end
    vd <= degree_limit || return (:untried, nothing, vd)
    e = radical_form(v; maxlen = maxlen, degree_limit = 0)   # the degree is already decided
    # with no length budget the only way back is a descent that could not certify a sign, which is not
    # the same statement as "longer than you asked for"
    e === nothing && return (maxlen > 0 ? :long : :failed, nothing, vd)
    return (:ok, e, vd)
end

"""
Printing an `ExactX` does no arithmetic of any kind.

The stored `P(x)·√(R(x))` is exact at every level and free to produce; the radical form is a question
the reader asks (`radical(v)`, which prints the number with it), and the certified `≈` is a Horner pass
this display used to pay for at every value — `IOContext(io, :approximate => true)` asks for it back,
and `Float64(v)` was always the direct way. Past the length budget the line becomes the *shape* of the
value with its sizes underneath, which is what a degree-208 polynomial has to say for itself. Nothing
else is printed: the display is the value, and the rest is what the reader calls for.
"""
function Base.show(io::IO, ::MIME"text/plain", v::ExactX)
    h = v.k + 2
    d = euler_phi(2h) ÷ 2
    lines = ["Exact value at level k = " * string(v.k) * "   (x = 2cos(π/" * string(h) *
             "), degree " * string(d) * " over ℚ)"]
    if iszero(v.p)
        print(io, lines[1], "\n  = 0")
        return
    end
    push!(lines, "  = " * _xform_str(v))
    if _xform_long(v)
        for l in _xform_dims(v)
            push!(lines, "    " * l)
        end
        push!(lines, "    `v.x_value` gives (P, R) in full")
    end
    if get(io, :approximate, false)
        nv = numeric_value(v)
        push!(lines, nv === nothing ?
              "  (the value cancels beyond " * string(EVAL_MAX_BITS) * " bits at x = 2cos(π/" *
              string(h) * "); for a number use q6j(…; k = " * string(v.k) * "))" :
              "  ≈ " * string(Float64(nv)))
    end
    print(io, join(lines, "\n"))
end

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
immediate proof of vanishing. `iszero` uses exact tests where available, but can fall back
to numerical comparison for multiple specialized radical classes.
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
Whether the sum is zero.

An empty term list is zero by construction, and a single term is zero exactly when its coefficient or its
radical vanishes — both free. With several terms the class keys may not be independent, so the norm
decides: a nonzero norm **proves** the value nonzero. A zero norm means *some* sign choice of the roots
vanishes, and which one is settled by evaluating; the candidates differ by the whole size of a term, so
the comparison is not delicate. This is the same contract the cyclotomic layer offers, at half the degree.
"""
function Base.iszero(s::ExactXSum)
    isempty(s.terms) && return true
    is_provably_nonzero(s) && return false
    length(s.terms) == 1 && return true              # its coefficient or its radical vanishes
    v = numeric_value(s)
    v === nothing && return false                    # cannot evaluate: keep the structural answer
    scale = zero(v)
    for S in keys(s.terms)
        t = numeric_value(_term(s, S))
        t === nothing && return false
        scale = max(scale, abs(t))
    end
    return abs(v) <= 1e-20 * max(scale, one(scale))
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

