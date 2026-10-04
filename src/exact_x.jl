# ---------------------------------------------------------------------------------
#  Exact level values in x = q + q⁻¹, and the display ladder
#
#  A level value is real, so it is stored in the real subfield: v = P(x)·√(R(x)) with x = 2cos(π/h),
#  h = k + 2, and P, R ∈ ℚ[x] reduced modulo Ψ_h (degree φ(2h)/2, half the cyclotomic degree). This is the
#  generic `XValue` of `generic_x.jl` reduced modulo Ψ_h as the sum is formed.
#
#  Display ladder: rational → single surd → nested radical → polynomial in x. A nested-radical form exists
#  iff φ(2h)/2 is a power of 2 (the real cyclotomic field is abelian); it is found by Lagrange's descent
#  u = (u + σu)/2 + √(((u − σu)/2)²), with σ_j: x ↦ C_j(x) the Galois action. Long radical forms fall back
#  to the polynomial.
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

"60-bit primes for the modular division below (more than needed, so unlucky primes cannot exhaust them)."
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

`num · den⁻¹ mod Ψ_h`, from images modulo word-sized primes, or `nothing` if the primes run out. The answer
is far smaller than the inputs, so this beats Euclidean inversion over ℚ[x] (13–23× at large labels). Each
reconstruction is checked exactly, `p·den ≡ num (mod Ψ)`, before it is returned.
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

`sqclass` lists the ψ indices under the root (empty when `v = P(x)`); a class that is a rational square at
this level is folded into `P`. Degrees are below `φ(2h)/2`.

Printing shows the stored `P(x)·√(R(x))`. [`radical`](@ref) gives nested square roots on request, or says
why there are none. Properties:

| | |
|---|---|
| `v.x_value` | the pair `(P, R)`, the value being `P(x)·√(R(x))` |
| `v.rad` | shorthand for `radical(v)` |

Lower-level accessors (not exported): `radical_form`, `xpolynomial`, `radicand`, `has_radical_form`.
"""
struct ExactX
    k::Int
    p::QQPolyRingElem          # P, reduced mod Ψ_h
    sqclass::Vector{Int}       # ψ indices under the root
    r::QQPolyRingElem          # R = ∏_{e ∈ sqclass} ψ_e, reduced mod Ψ_h

    function ExactX(k::Integer, p::QQPolyRingElem, cls::Vector{Int}, r::QQPolyRingElem)
        kk = Int(k)
        S, _ = _qqx()
        # At a level the radicand can reduce to a rational square: fold it into P, so equal values compare
        # equal and no `√(1)` is carried.
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
`v.x_value` is the named pair `(P, R)`, the value being `P(x)·√(R(x))` (`R = 1` when the class is empty).
`v.rad` is [`radical(v)`](@ref radical): nested square roots, or a [`NoRadical`](@ref) saying why not.
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
    # Reduce as the sum is formed. A reduced 0/0 is not a verdict: retry unreduced, where a shared factor
    # cancels exactly, before calling the rule singular.
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
Horner with the running mass `Σ|cᵢ||x|ⁱ`, which bounds the relative error, `(2d+1)·u·mass/|value|`. Needed:
near `x = 2` the polynomial cancels, more so at higher levels.
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
The lower cap for the Lagrange descent, which needs only signs, at every node. A node that cannot be
signed at 4096 bits abandons the descent (no radical form reported).
"""
const DESCENT_MAX_BITS = 4096

"""
Largest degree `radical_form` descends on under a length budget; beyond it the form cannot fit `maxlen`
(degree 16 took 8.5 s at k = 100). `radical_form(v; maxlen = 0)` ignores the cap.
"""
const RADICAL_MAX_DEGREE = 8

"""
Largest field degree at which `radical_levels` refines a level the sufficient condition rejects (an
unbounded sweep to kmax = 400 took 11.7 s).
"""
const REFINE_MAX_DEGREE = 64

"""
Largest value degree [`radical`](@ref) descends on by default. The descent is exponential in it (degree 64
exhausted memory); past the cap `radical` returns a [`NoRadical`](@ref), and `radical(v; degree_limit = …)`
lifts it.
"""
const RADICAL_DESCENT_MAX_DEGREE = 32


"""
Largest field degree at which a value's exact degree (a minimal polynomial over ℚ) is computed; beyond it
the cost grows quickly and the value is past the descent budget anyway.
"""
const RADICAL_MINPOLY_MAX_DEGREE = 64

"""
    _eval_at_x(f, h; rtol) -> BigFloat or nothing

`f(2cos(π/h))` at the precision that certifies `rtol`, or `nothing` past `EVAL_MAX_BITS`. The first pass
measures the cancellation and the second is wide enough for it (a fixed 256 bits gave −64.0 for 1/6 at
k = 420).
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

The number, at whatever precision certifies it, rounded once. Throws when none does; [`numeric_value`](@ref)
returns `nothing` instead.
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

The number an exact value denotes, in `T` (the same call as for the cyclotomic carrier).
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
