# Exact level values, continued from exact_x.jl: radical expressions, the Lagrange descent and `radical`.

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
        "radical_levels needs a function attached to a symbol rule: q6j, q3j_factorial, fsymbol, gsymbol"))
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
