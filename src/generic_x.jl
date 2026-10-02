# ---------------------------------------------------------------------------------
#  Exact values in the reciprocal variable x = q + q⁻¹
#
#  Every q-integer is a polynomial in x with integer coefficients, [n] = U_{n−1}(x/2), so a factorial rule
#  evaluates exactly in ℤ[x] with no root of unity, no level and no number field. Three things make this the
#  right carrier rather than one more backend.
#
#  1. ℤ[x] is a UFD, so a *square class* is canonical. Better: the irreducible factors of the q-integers are
#     explicit. Writing ψ_e for the fold of Φ_e (the minimal polynomial of 2cos(2π/e), monic of degree
#     φ(e)/2),
#
#         [n] = ∏_{e | 2n, e ≥ 3} ψ_e(x) ,
#
#     verified here for n = 1..24. So the multiplicity of each ψ_e in any product of q-factorials is a
#     divisor count, and a radical class is an 𝔽₂ exponent vector over an explicit basis — no polynomial
#     factoring, no number-field square testing, and none of the formal-key trouble of `CyclotomicMonomial`
#     (where Φ₃(q²) = q² at q = ζ₈ collapses a formal radical into the base field).
#
#  2. A level is a quotient: Ψ_h := ψ_{2h} is the minimal polynomial of 2cos(π/h) = x at q = e^{iπ/h}, monic
#     over ℤ of degree exactly φ(2h)/2 — half the cyclotomic field's degree. Reducing an integer polynomial
#     modulo a monic integer polynomial stays integral, so a level value's canonical form is a primitive
#     integer vector and equality is a vector comparison.
#
#  3. An identity proved here is proved for every q, hence for every level at once. The coherence identities
#     are *polynomial* in x once a common radical is divided out, so
#     `prove_identity` is a genuine proof rather than a level-by-level check.
#
#  What this is not: a route for large spin at a fixed level. Degrees grow like ~10j², so the cleared form
#  over ℤ[q]/(q^{2h}−1) stays the computational kernel there, with this module above it for comparison,
#  proof and display.
# ---------------------------------------------------------------------------------

const X_LOCK = ReentrantLock()
const _XRING = Ref{Any}()          # built on first use: FLINT objects are not precompiled into the image

"The polynomial ring ℤ[x] and its generator, x = q + q⁻¹."
function xring()
    isassigned(_XRING) && return _XRING[]::Tuple{ZZPolyRing,ZZPolyRingElem}
    lock(X_LOCK) do
        isassigned(_XRING) || (_XRING[] = polynomial_ring(ZZ, "x"))
    end
    return _XRING[]::Tuple{ZZPolyRing,ZZPolyRingElem}
end

const _QRING = Ref{Any}()

"An auxiliary ring ℤ[q], used only to build cyclotomic polynomials before folding them."
function _qring()
    isassigned(_QRING) && return _QRING[]::Tuple{ZZPolyRing,ZZPolyRingElem}
    lock(X_LOCK) do
        isassigned(_QRING) || (_QRING[] = polynomial_ring(ZZ, "q"))
    end
    return _QRING[]::Tuple{ZZPolyRing,ZZPolyRingElem}
end

const _CHEB = Ref{Any}()           # C_0, C_1, ... with C_n = qⁿ + q⁻ⁿ as polynomials in x
const _PSI = Dict{Int,Any}()       # ψ_e, e ≥ 3
const _QINT_X = Ref{Any}()         # [0], [1], [2], ...

"""
    cheb_x(n) -> ZZPolyRingElem

`C_n(x)` with `C_n = qⁿ + q⁻ⁿ`: `C₀ = 2`, `C₁ = x`, `C_{n+1} = x C_n − C_{n−1}`. This is `2T_n(x/2)`.
The basis combines each pair of reciprocal Laurent monomials into one term with the same coefficient,
giving a compact representation of reciprocal expressions.
"""
function cheb_x(n::Int)
    n >= 0 || throw(DomainError(n, "Chebyshev index must be nonnegative"))
    R, x = xring()
    lock(X_LOCK) do
        isassigned(_CHEB) || (_CHEB[] = Any[R(2), x])
        c = _CHEB[]::Vector{Any}
        while length(c) <= n
            push!(c, x * c[end] - c[end-1])
        end
        return c[n+1]::ZZPolyRingElem
    end
end

"""
    fold_palindromic(f) -> ZZPolyRingElem

A palindromic `f ∈ ℤ[q]` of even degree `2d` satisfies `q^{-d} f(q) = g(x)` for a unique `g ∈ ℤ[x]`; this
returns `g`. Monic in, monic out. Throws if `f` is not palindromic of even degree.
"""
function fold_palindromic(f)
    D = degree(f)
    iseven(D) || throw(ArgumentError("fold_palindromic needs even degree, got $D"))
    d = D ÷ 2
    for i in 0:d
        coeff(f, i) == coeff(f, D - i) ||
            throw(ArgumentError("fold_palindromic needs a palindromic polynomial"))
    end
    R, _ = xring()
    g = R(coeff(f, d))
    for n in 1:d
        g += coeff(f, d + n) * cheb_x(n)
    end
    return g
end

"""
    psi_x(e) -> ZZPolyRingElem

`ψ_e`, the fold of the cyclotomic polynomial `Φ_e`: the minimal polynomial of `2cos(2π/e)`, monic over ℤ of
degree `φ(e)/2`, irreducible. Defined for `e ≥ 3` (`Φ₁` and `Φ₂` are the antisymmetric factors that cancel
out of every balanced ratio and never appear in a q-integer).
"""
function psi_x(e::Int)
    e >= 3 || throw(DomainError(e, "ψ_e is defined for e ≥ 3"))
    lock(X_LOCK) do
        get!(_PSI, e) do
            _, q = _qring()
            fold_palindromic(cyclotomic(e, q))
        end
    end::ZZPolyRingElem
end

"""
    psi_level(h) -> ZZPolyRingElem

`Ψ_h = ψ_{2h}`, the minimal polynomial of `x = q + q⁻¹` at `q = e^{iπ/h}`, i.e. of `2cos(π/h)`. Monic over ℤ
of degree `φ(2h)/2`. Reduction modulo it is the level specialisation, and it keeps integer coefficients
because it is monic.
"""
psi_level(h::Int) = psi_x(2h)

"""
    qint_x(n) -> ZZPolyRingElem

`[n]` as a polynomial in `x = q + q⁻¹`: `[0] = 0`, `[1] = 1`, `[n+1] = x[n] − [n−1]`, i.e. `U_{n−1}(x/2)`.
Negative `n` gives `−[−n]`.
"""
function qint_x(n::Int)
    n < 0 && return -qint_x(-n)
    R, x = xring()
    lock(X_LOCK) do
        isassigned(_QINT_X) || (_QINT_X[] = Any[R(0), R(1)])
        c = _QINT_X[]::Vector{Any}
        while length(c) <= n
            push!(c, x * c[end] - c[end-1])
        end
        return c[n+1]::ZZPolyRingElem
    end
end

"`[n]!` in ℤ[x]."
function qfact_x(n::Int)
    n >= 0 || throw(DomainError(n, "factorial of a negative argument"))
    R, _ = xring()
    p = R(1)
    for m in 2:n
        p *= qint_x(m)
    end
    return p
end


# ---------------------------------------------------------------------------------
#  Radical classes as 𝔽₂ exponent vectors over the ψ basis
# ---------------------------------------------------------------------------------

"""
    psi_exponents(pairs) -> Dict{Int,Int}

Multiplicity of each `ψ_e` in `∏ₙ [n]!^{cₙ}`, from `pairs` a collection of `n => c`. Since
`[n]! = ∏_{m≤n} [m]` and `[m] = ∏_{e | 2m, e ≥ 3} ψ_e`, the multiplicity of `ψ_e` is a divisor count:
no factoring is involved.
"""
function psi_exponents(pairs)
    isempty(pairs) && return Dict{Int,Int}()
    N = maximum(Int(p.first) for p in pairs)
    N <= 0 && return Dict{Int,Int}()
    E = Dict{Int,Int}()
    # ψ_e divides [m] exactly when e divides 2m. Count those multiples in n! directly.
    for e in 3:2N
        step = iseven(e) ? e ÷ 2 : e
        v = sum(Int(c) * (Int(n) ÷ step) for (n,c) in pairs; init=0)
        iszero(v) || (E[e] = v)
    end
    return E
end

"""
Reduce modulo a monic `Ψ` when one is supplied. Every step of the ℤ[x] arithmetic below is a ring
operation, so reducing as it goes gives the same element of ℤ[x]/(Ψ) as reducing at the end, at a degree
that stays below `φ(2h)/2` instead of growing with the labels.
"""
_rx(f, m) = m === nothing ? f : mod(f, m)

# ---- reduced building blocks, cached per level ------------------------------------------------
#
# The operands of the reduced recursion do not have to be the full ℤ[x] polynomials. `[n]` has degree
# n−1, which at j = 20, k = 60 reaches 80 against a modulus of degree 30 — so two thirds of every such
# multiplication is thrown away by the reduction that follows it. Reducing each `[n]` and each `ψ_e`
# once per level and caching them makes every product a degree-<d by degree-<d one.
#
# The cache is keyed by the level, not by the modulus polynomial, so that a caller cannot silently share
# entries between different Ψ; `generic_value` only consults it when it is told the level.

const _RED_CACHE = LRU{Tuple{Symbol,Int,Int},Any}(maxsize = 8192)

function _qint_red(n::Int, h::Int)
    lock(X_LOCK) do
        get!(_RED_CACHE, (:qint, h, n)) do
            mod(qint_x(n), psi_level(h))
        end
    end::ZZPolyRingElem
end

function _psi_red(e::Int, h::Int)
    lock(X_LOCK) do
        get!(_RED_CACHE, (:psi, h, e)) do
            mod(psi_x(e), psi_level(h))
        end
    end::ZZPolyRingElem
end

"`ψ_e` (or `[n]`) as the reduced form when a level is known, and the plain one otherwise."
@inline _psi_at(e::Int, level) = level === nothing ? psi_x(e) : _psi_red(e, level)
@inline _qint_at(n::Int, level) = level === nothing ? qint_x(n) : _qint_red(n, level)

"""
`∏ ψ_e^{v_e}` with the work grouped by exponent: the factors sharing an exponent are multiplied once and
the result raised by binary powering, instead of one reduced multiplication per unit of exponent. At
j = 20, k = 60 the prefactor's exponents sum to 392 and there are 76 distinct factors, so the grouping is
most of the cost of building it.
"""
function _psi_pow_group(pairs, modulus, level)
    R, _ = xring()
    isempty(pairs) && return R(1)
    bym = Dict{Int,Vector{Int}}()
    for (e, m) in pairs
        m > 0 && push!(get!(bym, m, Int[]), e)
    end
    acc = R(1)
    for m in sort!(collect(keys(bym)))
        base = R(1)
        for e in bym[m]
            base = _rx(base * _psi_at(e, level), modulus)
        end
        acc = _rx(acc * _powmod_x(base, m, modulus), modulus)
    end
    return acc
end

"`f^n` by squaring, reducing at every step."
function _powmod_x(f, n::Int, modulus)
    R, _ = xring()
    n <= 0 && return R(1)
    result = R(1); base = f
    while n > 0
        isodd(n) && (result = _rx(result * base, modulus))
        n >>= 1
        n > 0 && (base = _rx(base * base, modulus))
    end
    return result
end

"""
Halve a ψ exponent vector, returning the square part's exponents and the odd-multiplicity indices — the
two halves of `√(∏ ψ_e^{E_e})`, without building any polynomial yet. Keeping it symbolic is the point:
the square part then merges with the first term's exponents before either is expanded.
"""
function _halve_exponents(E::Dict{Int,Int})
    half = Dict{Int,Int}(); rad = Int[]
    for e in sort!(collect(keys(E)))
        v = E[e]
        r = mod(v, 2)
        r == 1 && push!(rad, e)
        f = (v - r) ÷ 2
        iszero(f) || (half[e] = f)
    end
    return half, rad
end

"Merge `b` into `a`, summing exponents and dropping the cancellations."
function _merge_exponents!(a::Dict{Int,Int}, b::Dict{Int,Int})
    for (e, v) in b
        w = get(a, e, 0) + v
        iszero(w) ? delete!(a, e) : (a[e] = w)
    end
    return a
end

"Split `∏ ψ_e^{E_e}` as `(square root of the square part, odd-multiplicity indices)`."
function _split_sqrt(E::Dict{Int,Int}; modulus = nothing, level = nothing)
    rad = Int[]
    nums = Tuple{Int,Int}[]; dens = Tuple{Int,Int}[]
    for e in sort!(collect(keys(E)))
        v = E[e]
        r = mod(v, 2)
        f = (v - r) ÷ 2
        r == 1 && push!(rad, e)
        f > 0 ? push!(nums, (e, f)) : f < 0 && push!(dens, (e, -f))
    end
    return _psi_pow_group(nums, modulus, level), _psi_pow_group(dens, modulus, level), rad
end

"`∏ ψ_e^{E_e}` as a numerator/denominator pair."
function _psi_monomial(E::Dict{Int,Int}; modulus = nothing, level = nothing)
    nums = Tuple{Int,Int}[]; dens = Tuple{Int,Int}[]
    for e in sort!(collect(keys(E)))
        v = E[e]
        v > 0 ? push!(nums, (e, v)) : v < 0 && push!(dens, (e, -v))
    end
    return _psi_pow_group(nums, modulus, level), _psi_pow_group(dens, modulus, level)
end


# ---------------------------------------------------------------------------------
#  Exact values over ℤ[x]
# ---------------------------------------------------------------------------------


"""
    XValue

The exact value `√(∏_{e ∈ rad} ψ_e(x)) · num/den` with `num, den ∈ ℤ[x]`. Build one with `xvalue(rad, num, den)`,
which normalises: `rad` is a sorted list of distinct `e ≥ 3` — an 𝔽₂ exponent vector over the ψ basis, and
therefore a genuine square class in ℚ(x), not a formal key — `gcd(num, den) = 1`, `den` has positive
leading coefficient, and the pair is primitive. Normalised values are equal as numbers exactly when they are
equal as structs, so `==` is a decision procedure.

The branch of the square root is whichever the caller's construction fixes; `prove_identity` compares
expressions built the same way and is therefore branch-independent.
"""
struct XValue
    rad::Vector{Int}
    num::ZZPolyRingElem
    den::ZZPolyRingElem
end

"""
    xvalue(rad, num, den = 1) -> XValue

The `XValue` constructor — not to be confused with [`x_form`](@ref), which expands a symbolic value
into one. This spelling stays; the one-argument `xvalue(v)` is the deprecated name of `x_form`.

Normalising constructor: cancels `gcd(num, den)`, makes `den` primitive with positive leading coefficient,
and sorts the radical indices. A zero numerator gives the canonical zero.
"""
function xvalue(rad::AbstractVector{<:Integer}, num::ZZPolyRingElem,
                den::ZZPolyRingElem = one(parent(num)))
    iszero(den) && throw(DivideError())
    R, _ = xring()
    iszero(num) && return XValue(Int[], R(0), R(1))
    # 𝔽₂ semantics: a repeated index is a square and moves into the numerator (√ψ·√ψ = ψ),
    # it does not cancel. Only the odd-multiplicity indices stay under the root.
    keep = Int[]
    if !isempty(rad)
        cnt = Dict{Int,Int}()
        for e in rad
            e >= 3 || throw(DomainError(e, "radical indices are ψ indices, e ≥ 3"))
            cnt[e] = get(cnt, e, 0) + 1
        end
        for e in sort!(collect(keys(cnt)))
            c = cnt[e]
            isodd(c) && push!(keep, e)
            f = c ÷ 2
            f > 0 && (num *= psi_x(e)^f)
        end
    end
    g = gcd(num, den)
    isone(g) || (num = divexact(num, g); den = divexact(den, g))
    if leading_coefficient(den) < 0
        num = -num; den = -den
    end
    c = gcd(content(num), content(den))
    isone(c) || (num = divexact(num, c); den = divexact(den, c))
    return XValue(keep, num, den)
end

function Base.zero(::Type{XValue})
    R, _ = xring()
    return XValue(Int[], R(0), R(1))
end
function Base.one(::Type{XValue})
    R, _ = xring()
    return XValue(Int[], R(1), R(1))
end
Base.iszero(v::XValue) = iszero(v.num)
Base.:(==)(a::XValue, b::XValue) =
    iszero(a) ? iszero(b) : (!iszero(b) && a.rad == b.rad && a.num == b.num && a.den == b.den)
Base.hash(v::XValue, h::UInt) = hash(v.rad, hash(v.num, hash(v.den, h)))
Base.:-(a::XValue) = XValue(a.rad, -a.num, a.den)

"Is the value free of radicals, i.e. an element of ℚ(x)?"
is_rational_x(v::XValue) = isempty(v.rad)

function Base.:*(a::XValue, b::XValue)
    (iszero(a) || iszero(b)) && return zero(XValue)
    shared = intersect(a.rad, b.rad)                  # √A √B = (∏ shared ψ) √(A △ B)
    R, _ = xring()
    extra = R(1)
    for e in shared
        extra *= psi_x(e)
    end
    return xvalue(symdiff(a.rad, b.rad), a.num * b.num * extra, a.den * b.den)
end

function Base.inv(a::XValue)
    iszero(a) && throw(DivideError())
    R, _ = xring()
    rp = R(1)
    for e in a.rad
        rp *= psi_x(e)                                 # 1/√A = √A / A
    end
    return xvalue(a.rad, a.den, a.num * rp)
end
Base.:/(a::XValue, b::XValue) = a * inv(b)

function Base.:^(a::XValue, n::Integer)
    n < 0 && return inv(a)^(-n)
    r = one(XValue)
    for _ in 1:n
        r = r * a
    end
    return r
end

"Sum of two values with the *same* square class; use [`XSum`](@ref) when they differ."
function Base.:+(a::XValue, b::XValue)
    iszero(a) && return b
    iszero(b) && return a
    a.rad == b.rad || throw(ArgumentError(
        "cannot add values with different square classes $(a.rad) and $(b.rad); use XSum"))
    return xvalue(a.rad, a.num * b.den + b.num * a.den, a.den * b.den)
end
Base.:-(a::XValue, b::XValue) = a + (-b)

Base.:*(c::Integer, a::XValue) = iszero(c) ? zero(XValue) : xvalue(a.rad, c * a.num, a.den)
Base.:*(a::XValue, c::Integer) = c * a


# ---------------------------------------------------------------------------------
#  Sums of values with different square classes
# ---------------------------------------------------------------------------------

"""
    XSum

A formal sum `Σ_S c_S √(∏_{e ∈ S} ψ_e)` with each `c_S ∈ ℚ(x)`, keyed by the square class `S`. Because the
keys are genuine square classes in ℚ(x) and distinct classes give ℚ(x)-linearly independent square roots,
`iszero` is a **decision procedure** here — unlike the cyclotomic analogue, whose formal keys can collapse
after specialisation.
"""
struct XSum
    terms::Dict{Vector{Int},XValue}
end

XSum() = XSum(Dict{Vector{Int},XValue}())
function XSum(v::XValue)
    iszero(v) && return XSum()
    return XSum(Dict(v.rad => v))
end
Base.iszero(s::XSum) = isempty(s.terms)
Base.length(s::XSum) = length(s.terms)
Base.:(==)(a::XSum, b::XSum) = a.terms == b.terms

function Base.:+(a::XSum, b::XSum)
    out = copy(a.terms)
    for (k, v) in b.terms
        if haskey(out, k)
            w = out[k] + v
            iszero(w) ? delete!(out, k) : (out[k] = w)
        else
            out[k] = v
        end
    end
    return XSum(out)
end
Base.:-(a::XSum) = XSum(Dict(k => -v for (k, v) in a.terms))
Base.:-(a::XSum, b::XSum) = a + (-b)
Base.:+(a::XSum, b::XValue) = a + XSum(b)
Base.:+(a::XValue, b::XSum) = XSum(a) + b

function Base.:*(a::XSum, b::XSum)
    out = XSum()
    for (_, u) in a.terms, (_, v) in b.terms
        out = out + XSum(u * v)
    end
    return out
end
Base.:*(a::XSum, b::XValue) = a * XSum(b)
Base.:*(a::XValue, b::XSum) = XSum(a) * b


# ---------------------------------------------------------------------------------
#  Evaluating a factorial rule over ℤ[x]
# ---------------------------------------------------------------------------------

"""
    generic_value(s::FactorialSum) -> XValue

The exact value of a factorial rule as an element of `ℚ(x)` times one square root, `x = q + q⁻¹`.

The prefactor and the first term are products of q-factorials, so they never need polynomial arithmetic:
their ψ multiplicities are divisor counts ([`psi_exponents`](@ref)), and the square root of the prefactor
splits into a ψ-monomial times the odd-multiplicity class. Only the sum itself is a genuine polynomial
computation, done by the cleared ratio recursion

    P ← B·Q ± A·P,   Q ← B·Q,

with `A/B` the adjacent-term ratio as a product of q-integers — the same recursion as the level kernel in
`exact_cleared.jl`, over ℤ[x] instead of ℤ[q]/(q^{2h}−1). There is no level, no root of unity and no
number field anywhere in this path.
"""
function generic_value(s::FactorialSum; modulus = nothing, level = nothing)
    R, _ = xring()
    is_empty_sum(s) && return zero(XValue)

    # --- prefactor: √(∏[n]!^c) or ∏[n]!^c, as ψ exponents and a square class ---
    #
    # The exponents of the prefactor and of the first term are merged *before* either becomes a
    # polynomial. They overlap heavily — measured, the separate expansions need three times as many
    # reduced multiplications as the merged one at every size tried — and a ψ that appears in both
    # simply cancels instead of being built twice and divided out afterwards.
    Epre = psi_exponents(collect(s.pre))
    Enet, rad = s.sqrt_pre ? _halve_exponents(Epre) : (copy(Epre), Int[])

    # --- first term ∏[arg(f, zlo)]!^c ---
    first_pairs = Pair{Int,Int}[]
    for f in s.fac
        Int(f.c) == 0 && continue
        n = _arg(f, s.zlo)
        n >= 0 || throw(DomainError(n, "negative factorial argument in the rule"))
        push!(first_pairs, n => Int(f.c))
    end
    _merge_exponents!(Enet, psi_exponents(first_pairs))
    pn, pd = _psi_monomial(Enet; modulus = modulus, level = level)

    # --- the sum, by the cleared ratio recursion; Q stays a product of q-integers ---
    P = R(1); Q = R(1)
    rsign = s.alternating ? -1 : 1
    for z in (s.zhi - 1):-1:s.zlo
        A = R(1); B = R(1)
        for f in s.fac
            lo, hi, c = _factor_step(f, z)
            (c == 0 || lo > hi) && continue
            blk = R(1)
            for n in lo:hi
                blk = _rx(blk * _qint_at(n, level), modulus)
            end
            if c > 0
                A = _rx(A * _powmod_x(blk, c, modulus), modulus)
            else
                B = _rx(B * _powmod_x(blk, -c, modulus), modulus)
            end
        end
        BQ = _rx(B * Q, modulus)
        P = _rx(BQ + rsign * A * P, modulus)
        Q = BQ
    end

    sgn = Int(s.sign0) * ((s.alternating && isodd(s.zlo)) ? -1 : 1)
    return xvalue(rad, _rx(sgn * _rx(pn * P, modulus), modulus), _rx(pd * Q, modulus))
end

"The generic-q exact value of one 6j symbol, from doubled labels."
generic_sixj(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int) =
    generic_value(sixj_sum(J1, J2, J3, J4, J5, J6))


# ---------------------------------------------------------------------------------
#  Specialisation to a level, and to a number
# ---------------------------------------------------------------------------------

"""
    at_level(v::XValue, k) -> (rad, num, den)

Reduce a generic value modulo `Ψ_{k+2}`, i.e. specialise `x` to `2cos(π/(k+2))`. All three parts stay in
ℤ[x] because `Ψ_h` is monic, and the result is the canonical form of the level value: `num` and `den` are
integer vectors of length at most `φ(2h)/2`, half the cyclotomic field's degree.

A vanishing `den` means the rule is singular at this level; `iszero(num)` with a nonzero `den` is an exact
zero. The square class may *degenerate* under the specialisation (a `rad` that is squarefree in ℤ[x] can
become a square modulo `Ψ_h`), which is why two values with different `rad` must not be declared distinct
at a fixed level without checking that.
"""
function at_level(v::XValue, k::Integer)
    h = Int(k) + 2
    Ψ = psi_level(h)
    R, _ = xring()
    rp = R(1)
    for e in v.rad
        rp *= psi_x(e)
    end
    return (rad = mod(rp, Ψ), num = mod(v.num, Ψ), den = mod(v.den, Ψ))
end

"""
    chebyshev(f) -> Vector{ZZRingElem}

Coefficients of `f ∈ ℤ[x]` in the basis `C_n = qⁿ + q⁻ⁿ`, with the convention that the constant term is the
coefficient of `1` rather than of `C₀ = 2`. Grouping reciprocal powers into one term gives a compact
basis for displaying exact values.
"""
function chebyshev(f::ZZPolyRingElem)
    d = degree(f)
    d < 0 && return ZZRingElem[]
    out = [coeff(f, i) for i in 0:d]
    for n in d:-1:1                      # peel the leading Chebyshev term, highest first
        c = out[n+1]
        iszero(c) && continue
        Cn = cheb_x(n)
        for i in 0:n-1
            out[i+1] -= c * coeff(Cn, i)
        end
    end
    return out
end

"""
Numerical value of `f ∈ ℤ[x]` at the level `k`, i.e. at `x = 2cos(π/(k+2))`.

**This is the one ill-conditioned operation in the module.** The exact algebra is cancellation-free, but a
high-degree integer polynomial evaluated near `x = 2` cancels catastrophically — the conditioning that the
Racah sum has in the *terms* reappears here in the *coefficients*. Use `BigFloat` with enough precision
(degree × coefficient bits is the right scale) whenever the degree is more than a few dozen; `Float64` is
adequate only for small labels, and is provided for spot checks.
"""
function evaluate_at_level(f::ZZPolyRingElem, k::Integer, ::Type{T} = Float64) where {T}
    xv = 2 * cos(T(π) / T(Int(k) + 2))
    acc = zero(T)
    for i in degree(f):-1:0
        acc = acc * xv + T(coeff(f, i))
    end
    return acc
end

"Numerical value of an `XValue` at the level `k` (the positive branch of the square root)."
function evaluate_at_level(v::XValue, k::Integer, ::Type{T} = Float64) where {T}
    iszero(v) && return zero(T)
    r = one(T)
    for e in v.rad
        r *= evaluate_at_level(psi_x(e), k, T)
    end
    return sqrt(r) * evaluate_at_level(v.num, k, T) / evaluate_at_level(v.den, k, T)
end


# ---------------------------------------------------------------------------------
#  Identities, proved once for every level
# ---------------------------------------------------------------------------------

"""
    exceptional_levels(v::XValue; kmax = 4096) -> Vector{Int}

Levels `k ≤ kmax` at which the value's denominator vanishes, so that the generic expression says nothing
there. Because `Ψ_h = ψ_{2h}` is irreducible, `h` is exceptional exactly when `ψ_{2h}` divides the
denominator; the candidates are the `e = 2h` that occur in the denominator's ψ-content, so the search is a
short exact divisibility test rather than a scan.
"""
function exceptional_levels(v::XValue; kmax::Int = 4096)
    out = Int[]
    iszero(v) && return out
    D = degree(v.den)
    D <= 0 && return out
    for h in 2:(kmax + 2)
        deg = euler_phi(2h) ÷ 2
        deg > D && continue
        iszero(mod(v.den, psi_x(2h))) && push!(out, h - 2)
    end
    return out
end

"""
    prove_identity(lhs, rhs; kmax = 4096) -> NamedTuple

Decide `lhs == rhs` for `XValue` or `XSum` arguments, in ℤ[x]. The verdict is `:proved` or `:refuted`, and a
`:proved` verdict is a statement about **every** `q`, hence about every level at once: a polynomial identity
in `x` specialises to `x = 2cos(π/h)` for every `h`.

Returned fields:

  * `verdict`   — `:proved` or `:refuted`
  * `degree`    — the degree of the compared numerator, i.e. the size of the certificate
  * `classes`   — the square classes involved
  * `exceptional` — levels at which a denominator vanishes, where the identity is silent rather than false
  * `witness`   — for `:refuted`, the nonzero difference

Equality of two normalised `XValue`s is a struct comparison, so no tolerance, prime or bound enters. For
`XSum` the difference is taken termwise; the keys are genuine square classes and distinct classes give
linearly independent roots over ℚ(x), so an empty difference is a proof.
"""
function prove_identity(lhs::XValue, rhs::XValue; kmax::Int = 4096)
    d = nothing
    same = lhs == rhs
    if !same
        d = try
            lhs - rhs
        catch e
            e isa ArgumentError ? nothing : rethrow()   # different square classes: not equal
        end
    end
    exc = sort!(unique(vcat(exceptional_levels(lhs; kmax = kmax),
                            exceptional_levels(rhs; kmax = kmax))))
    return (verdict = same ? :proved : :refuted,
            degree = max(degree(lhs.num), degree(rhs.num)),
            classes = sort!(unique(vcat(lhs.rad, rhs.rad))),
            exceptional = exc,
            witness = same ? nothing : d)
end

function prove_identity(lhs::XSum, rhs::XSum; kmax::Int = 4096)
    d = lhs - rhs
    exc = Int[]
    for s in (lhs, rhs), (_, v) in s.terms
        append!(exc, exceptional_levels(v; kmax = kmax))
    end
    deg = 0
    for s in (lhs, rhs), (_, v) in s.terms
        deg = max(deg, degree(v.num))
    end
    return (verdict = iszero(d) ? :proved : :refuted,
            degree = deg,
            classes = sort!(unique(vcat([k for k in keys(lhs.terms)]...,
                                        [k for k in keys(rhs.terms)]...))),
            exceptional = sort!(unique(exc)),
            witness = iszero(d) ? nothing : d)
end
prove_identity(lhs::XSum, rhs::XValue; kw...) = prove_identity(lhs, XSum(rhs); kw...)
prove_identity(lhs::XValue, rhs::XSum; kw...) = prove_identity(XSum(lhs), rhs; kw...)

"""
    verify_by_evaluation(lhs::XValue, rhs::XValue) -> Bool

An independent check of a `:proved` verdict that uses no polynomial algebra: a polynomial identity of degree
`D` is settled by agreement at `D + 1` distinct points, and here the evaluation is exact integer arithmetic.
Only valid when the two sides carry the same square class.
"""
function verify_by_evaluation(lhs::XValue, rhs::XValue)
    lhs.rad == rhs.rad || return false
    f = lhs.num * rhs.den - rhs.num * lhs.den
    D = max(degree(f), 0)
    for t in 0:(D + 1)
        acc = ZZ(0)
        xv = ZZ(t)
        for i in degree(f):-1:0
            acc = acc * xv + coeff(f, i)
        end
        iszero(acc) || return false
    end
    return true
end


# ---------------------------------------------------------------------------------
#  Triangle admissibility without a level — used by the generic machinery above and by the coherence
#  examples, which live in the tutorials and the test suite rather than here.
# ---------------------------------------------------------------------------------

"Triangle admissibility for doubled labels, without a level."
@inline _tri_x(a::Int, b::Int, c::Int) = a + b >= c && b + c >= a && c + a >= b && iseven(a + b + c)
@inline _tet_x(a::Int, b::Int, c::Int, d::Int, e::Int, f::Int) =
    _tri_x(a, b, c) && _tri_x(a, e, f) && _tri_x(d, b, f) && _tri_x(d, e, c)
