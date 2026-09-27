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

`rad` lists the ψ indices under the root — an 𝔽₂ square class over the ψ basis of `generic_x.jl`, empty
when the value needs no root at all, in which case `v = P(x)` outright. Degrees are below `φ(2h)/2`, half
the degree of the cyclotomic field the same value occupies in `Exact(k)`'s default form.

Printing runs the display ladder; [`radical_form`](@ref) asks for the closed form in radicals on its own,
and [`xpolynomial`](@ref) for the coefficients of `P`.
"""
struct ExactX
    k::Int
    p::QQPolyRingElem          # P, reduced mod Ψ_h
    rad::Vector{Int}           # ψ indices under the root
    r::QQPolyRingElem          # R = ∏_{e ∈ rad} ψ_e, reduced mod Ψ_h
end

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
the same cleared ratio recursion as [`generic_value`](@ref), reduced modulo `Ψ_h` at every step so the
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
    for e in v.rad
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
    _eval_at_x(f, h; rtol) -> BigFloat or nothing

`f(2cos(π/h))`, at whatever precision it takes to certify `rtol`, or `nothing` when `EVAL_MAX_BITS` is not
enough. The precision is not guessed twice: the first pass measures the actual cancellation and the second
is made wide enough for it.

This exists because a fixed precision is wrong. At `k = 420` the value of `{1 1 1; 1 1 1}` is a degree-208
polynomial in `x = 2cos(π/422) = 1.99989…`, where 256 bits print **−64.011** for a symbol whose value is
**0.16666**. The expression was exact; only the number under it was not. An exact value that cannot say
what its number is should say nothing, not something.
"""
function _eval_at_x(f, h::Int; rtol::Float64 = 1e-16, maxbits::Int = EVAL_MAX_BITS)
    degree(f) < 0 && return BigFloat(0)
    bits = 256
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
    r = numeric_value(v)
    r === nothing && throw(DomainError(v.k,
        "the exact value at level $(v.k) cancels beyond $(EVAL_MAX_BITS) bits when evaluated at " *
        "x = 2cos(π/$(v.k + 2)); the expression is still exact. For a number use the certified numeric " *
        "path, `q6j(labels...; k = $(v.k))`."))
    return T(r)
end
Base.Float64(v::ExactX) = float(v, Float64)
Base.BigFloat(v::ExactX) = float(v, BigFloat)

"""
    numeric_value(v::ExactX) -> BigFloat or nothing

The number `v` denotes, certified, or `nothing` when no precision below `EVAL_MAX_BITS` certifies it. The
display uses this and simply omits the `≈` line in that case: the package has certified numeric paths of
its own, and a wrong number beside an exact expression is worse than no number at all.
"""
function numeric_value(v::ExactX)
    iszero(v.p) && return BigFloat(0)
    h = v.k + 2
    p = _eval_at_x(v.p, h)
    p === nothing && return nothing
    isempty(v.rad) && return p
    r = _eval_at_x(v.r, h)
    r === nothing && return nothing
    r < 0 && (r = zero(r))          # a real value forces R ≥ 0; guard the last bit of rounding
    return setprecision(() -> p * sqrt(r), BigFloat, max(precision(p), precision(r)))
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
    return count_ones(_value_degree(_square(v), v.k + 2)) == 1
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
while the field degree is at most `refine_max_degree` (default [`REFINE_MAX_DEGREE`](@ref)), and
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
        if iszero(v.p) || (degree(v.p) <= 0 && isempty(v.rad)) || has_radical_form(v)
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
three levels.

```julia
radical_form(q6j(Exact(3; form = :x), 1, 1, 1, 1, 1, 1))   # −(3 − √5)/2, the Fibonacci level
```
"""
function radical_form(v::ExactX; maxlen::Int = 80)
    iszero(v.p) && return RadExpr(0)
    # The *value*, not the level: a symbol can be rational — or quadratic — at a level whose field has an
    # odd prime in its degree, and refusing on the level alone hid exactly those. Measured, 12 of 1632
    # sampled values were rational at k = 5, 7, 9 and were being told no radical form existed.
    has_radical_form(v) || return nothing
    h = v.k + 2
    Ψ = _psiq(h)
    # A degree-2ᵐ value descends to 2ᵐ rational leaves, so anything past a handful of them is already far
    # longer than `maxlen` and the descent would only be paying to be thrown away. `maxlen = 0`, which
    # asks for the expression whatever its size, skips this.
    maxlen > 0 && _value_degree(_square(v), h) > RADICAL_MAX_DEGREE && return nothing
    reps = _conj_reps(h)
    e = nothing
    if isempty(v.rad)
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

"`_xpoly_str`, or a description of it when writing it out would be worse than useless."
function _xpoly_show(f)
    degree(f) > XPOLY_MAX_DEGREE && return _xpoly_summary(f)
    str = _xpoly_str(f)
    return length(str) > XPOLY_MAX_CHARS ? _xpoly_summary(f) : str
end

"The `P(x)·√(R(x))` line, always available."
function _xform_str(v::ExactX)
    iszero(v.p) && return "0"
    ps = _xpoly_show(v.p)
    isempty(v.rad) && return ps
    rs = _xpoly_show(v.r)
    ps == "1" && return "√(" * rs * ")"
    ps == "-1" && return "−√(" * rs * ")"
    need = occursin(" ", ps) && !startswith(ps, "(")
    return (need ? "(" * ps * ")" : ps) * " · √(" * rs * ")"
end

function Base.show(io::IO, v::ExactX)
    print(io, "ExactX(k = ", v.k, ", ", _xform_str(v), ")")
end

function Base.show(io::IO, ::MIME"text/plain", v::ExactX)
    h = v.k + 2
    d = euler_phi(2h) ÷ 2
    println(io, "Exact value at level k = ", v.k, "   (x = 2cos(π/", h, "), degree ", d, " over ℚ)")
    if iszero(v.p)
        println(io, "  = 0")
        return
    end
    rf = radical_form(v)
    xs = _xform_str(v)
    rs = rf === nothing ? nothing : _rad_str(rf)
    rs === nothing || println(io, "  = ", rs)
    rs == xs || println(io, "  = ", xs)
    if rf === nothing
        if !has_radical_form(v)
            println(io, "  no radical form: this value generates a subfield of ℚ(2cos(π/", h,
                        ")) whose degree is not")
            println(io, "  a power of two, so no real nested-square-root expression exists ",
                        "(casus irreducibilis)")
        else
            println(io, "  a radical form exists but is too long to be useful; ",
                        "`radical_form(v; maxlen = 0)` forces it")
        end
    end
    nv = numeric_value(v)
    if nv === nothing
        print(io, "  (the value cancels beyond ", EVAL_MAX_BITS,
                  " bits at x = 2cos(π/", h, "); for a number use q6j(…; k = ", v.k, "))")
    else
        print(io, "  ≈ ", Float64(nv))
    end
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

Base.:-(v::ExactX) = ExactX(v.k, -v.p, v.rad, v.r)
function _one_exactx(k::Integer)
    S, _ = _qqx()
    return ExactX(Int(k), S(1), Int[], S(1))
end
Base.one(::Type{ExactX}, k::Integer) = _one_exactx(k)

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
    inter = intersect(a.rad, b.rad)
    diff = sort!(symdiff(a.rad, b.rad))
    G = _psi_prod(inter, h)
    iszero(G) && return zero(ExactX, a.k)          # a radical that vanishes at this level
    sg = prod(e -> _psi_sign(e, h), inter; init = 1)
    p = _redq(_redq(a.p * b.p, Ψ) * (sg < 0 ? -G : G), Ψ)
    iszero(p) && return zero(ExactX, a.k)
    return ExactX(a.k, p, diff, _psi_prod(diff, h))
end

Base.:*(c::Union{Integer,Rational}, v::ExactX) =
    iszero(c) ? zero(ExactX, v.k) : ExactX(v.k, v.p * QQ(c), v.rad, v.r)
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
    return ExactX(v.k, q, copy(v.rad), v.r)
end
Base.:/(a::ExactX, b::ExactX) = a * inv(b)

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
free proof of vanishing. [`iszero`](@ref) asks for a real one.
"""
struct ExactXSum
    k::Int
    terms::Dict{Vector{Int},QQPolyRingElem}
end

ExactXSum(k::Integer) = ExactXSum(Int(k), Dict{Vector{Int},QQPolyRingElem}())
function ExactXSum(v::ExactX)
    iszero(v.p) && return ExactXSum(v.k)
    return ExactXSum(v.k, Dict(copy(v.rad) => v.p))
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
Base.:+(a::ExactX, b::ExactX) = ExactXSum(a) + ExactXSum(b)
Base.:-(a::ExactX, b::ExactX) = ExactXSum(a) - ExactXSum(b)

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

function Base.show(io::IO, s::ExactXSum)
    isempty(s.terms) && return print(io, "ExactXSum(k = ", s.k, ", 0)")
    parts = [_xform_str(_term(s, S)) for S in sort!(collect(keys(s.terms)))]
    print(io, "ExactXSum(k = ", s.k, ", ", join(parts, " + "), ")")
end

function Base.show(io::IO, ::MIME"text/plain", s::ExactXSum)
    println(io, "Exact sum at level k = ", s.k, "   (", length(s.terms),
                length(s.terms) == 1 ? " square class" : " square classes", ", x = 2cos(π/", s.k + 2, "))")
    if isempty(s.terms)
        print(io, "  = 0")
        return
    end
    for S in sort!(collect(keys(s.terms)))
        println(io, "  + ", _xform_str(_term(s, S)))
    end
    v = numeric_value(s)
    print(io, v === nothing ? "  (the value could not be certified numerically)" : "  ≈ " * string(Float64(v)))
end
