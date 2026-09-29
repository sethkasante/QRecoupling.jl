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
