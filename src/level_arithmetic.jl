# Arithmetic of level values in plain Julia, continued from level_field.jl.
#
# Products use the fusion rule. With A = (ζ − ζ⁻¹)a = Σ a_i (ζ^i − ζ^(−i)) and b = Σ_m s_m ζ^m, whose coefficients
# s_m are sums over every second coordinate of b, the product has coordinates r_c = Σ_i a_i (s_|c−i| − s_(c+i)).
# Labels above h/2 fold back by [h − n] = [n], and those above d reduce through the minimal polynomial.
#
# Square roots stay positive, so a product never decides a sign. Sums are kept by square class, and their zero test
# uses outward enclosures: a nonzero algebraic integer of degree D whose conjugates are at most B has size at
# least B^(1−D), so an enclosure entirely below that resolves it.

# ---- the basis ----

"Coordinates on [1], …, [d] of [d+1], …, [h÷2], from Φ_2h(ζ) = ζ^d Ψ(ζ + 1/ζ) and x[n] = [n+1] + [n−1]."
function _basis_rows(F::LevelField)
    (F.h ÷ 2 <= F.d || !isempty(F.basis)) && return F.basis
    lock(LEVEL_FIELD_LOCK) do
        isempty(F.basis) || return F.basis
        d = F.d
        phi = _level_phi(F.h - 2)
        # Ψ = c_d[1] + Σ_j c_(d+j) ([j+1] − [j−1]), monic in [d+1]
        top = zeros(BigInt, d)
        top[1] -= phi[d+1]
        for j in 1:d
            c = phi[d+j+1]
            j <= d - 1 && (top[j+1] -= c)
            j >= 2 && (top[j-1] += c)
        end
        rows = Vector{Vector{BigInt}}()
        push!(rows, top)
        for n in d+2:F.h÷2
            # [n] = x[n−1] − [n−2], with x[m] = [m+1] + [m−1] and [d+1] replaced by its coordinates
            prev = rows[end]; next = zeros(BigInt, d)
            for m in 1:d
                c = prev[m]
                iszero(c) && continue
                m >= 2 && (next[m-1] += c)
                if m < d
                    next[m+1] += c
                else
                    next .+= c .* top
                end
            end
            if n - 2 == d
                next[d] -= 1
            elseif n - 2 < d
                next[n-2] -= 1
            else
                next .-= rows[n-2-d]
            end
            push!(rows, next)
        end
        F.basisbits = maximum(ndigits(c; base = 2) for a in rows for c in a)
        F.basis = [[j => a[j] for j in 1:d if !iszero(a[j])] for a in rows]
        _trim_level_fields!(F)
        return F.basis
    end
end

"The integer vector over its common denominator in lowest terms, with a positive denominator."
function _lowest(num::Vector{BigInt}, den::BigInt)
    g = abs(den)
    for a in num
        isone(g) && break
        iszero(a) || (g = gcd(g, a))
    end
    den < 0 && (g = -g)
    isone(g) && return num, den
    return [a ÷ g for a in num], den ÷ g
end

function _level_make(k::Int, num::Vector{BigInt}, den::BigInt, r::BigInt, rad::Vector{Int})
    n, d = _lowest(num, den)
    all(iszero, n) && return zero(LevelValue, k)
    return LevelValue(k, n, d, r, rad)
end

Base.zero(::Type{LevelValue}, k::Integer) =
    LevelValue(Int(k), [big(0) for _ in 1:_totient(2(Int(k) + 2))÷2], big(1), big(1), Int[])
Base.zero(v::LevelValue) = zero(LevelValue, v.k)
function Base.one(::Type{LevelValue}, k::Integer)
    v = zero(LevelValue, k)
    v.num[1] = big(1)
    return v
end
Base.one(v::LevelValue) = one(LevelValue, v.k)
Base.isone(v::LevelValue) = isone(v.den) && isone(v.r) && isempty(v.rad) && isone(v.num[1]) && all(iszero, @view v.num[2:end])

"The product when every intermediate fits a word pair: most values at small labels."
function _label_mul_small(a::Vector{BigInt}, b::Vector{BigInt}, F::LevelField)
    d = F.d; h = F.h; top = h ÷ 2
    s = Int128.(b)
    for m in d-2:-1:1
        @inbounds s[m] += s[m+2]
    end
    r = zeros(Int128, 2d)
    @inbounds for i in 1:d
        iszero(a[i]) && continue
        ai = Int128(a[i])
        for m in 0:d-1
            t = ai * s[m+1]
            r[i+m] += t
            if m > 0
                c = i - m
                c > 0 && (r[c] += t)
                c < 0 && (r[-c] -= t)
            end
        end
    end
    @inbounds for c in top+1:2d
        c < h && (r[h-c] += r[c])
    end
    if top > d
        rows = _basis_rows(F)
        @inbounds for n in d+1:min(top, 2d)
            c = r[n]
            iszero(c) && continue
            for (j, a2) in rows[n-d]
                r[j] += c * Int128(a2)
            end
        end
    end
    return BigInt[big(r[n]) for n in 1:d]
end

"Exact product of two integer coordinate vectors."
function _label_mul(a::Vector{BigInt}, b::Vector{BigInt}, F::LevelField)
    d = F.d; h = F.h
    count(!iszero, a) > count(!iszero, b) && ((a, b) = (b, a))
    F.h ÷ 2 > d && _basis_rows(F)
    bits(v) = maximum(c -> ndigits(c; base = 2), v; init = 0)
    # 2d terms of size a·s, each s a sum of d coordinates, and a second round through the basis
    bits(a) + bits(b) + 2 * ndigits(2d; base = 2) + F.basisbits + ndigits(max(F.h ÷ 2 - d, 1); base = 2) <= 124 &&
        return _label_mul_small(a, b, F)
    s = copy(b)                                     # s[m+1] is the coefficient of ζ^m in b
    for m in d-2:-1:1
        s[m] += s[m+2]
    end
    r = [big(0) for _ in 1:2d]
    t = big(0)
    @inbounds for i in 1:d
        ai = a[i]
        iszero(ai) && continue
        for m in 0:d-1
            sm = s[m+1]
            iszero(sm) && continue
            MPZ.mul!(t, ai, sm)
            c = i + m                               # ζ^(i+m), and ζ^(i−m) when m is not zero
            MPZ.add!(r[c], r[c], t)
            if m > 0
                c = i - m
                c > 0 && MPZ.add!(r[c], r[c], t)
                c < 0 && MPZ.sub!(r[-c], r[-c], t)
            end
        end
    end
    top = h ÷ 2
    for c in top+1:2d
        iszero(r[c]) || c >= h || (r[h-c] += r[c])
    end
    out = r[1:d]
    if top > d
        rows = _basis_rows(F)
        for n in d+1:min(top, 2d)
            c = r[n]
            iszero(c) && continue
            for (j, a2) in rows[n-d]
                out[j] += c * a2
            end
        end
    end
    return out
end

"The value of integer coordinates at the place `i`, by the recurrence of the q-integers."
function _place_coordinates_value(t::LevelPlaces, a::Vector{UInt64}, i::Int)
    m = t.m; X = t.X[i]; b1 = zero(UInt64); b2 = zero(UInt64)
    @inbounds for n in length(a):-1:1
        b0 = mont_sub(m, mont_add(m, a[n], mont_mul(m, X, b1)), b2)
        b2 = b1; b1 = b0
    end
    return b1
end

"""
Coordinates of the inverse of integer coordinates, over a common denominator, from its images modulo word primes.
A candidate is multiplied back exactly before it is returned.
"""
function _label_inverse(a::Vector{BigInt}, F::LevelField)
    all(iszero, a) && throw(DivideError())
    d = F.d
    if all(iszero, @view a[2:end])
        return [n == 1 ? big(sign(a[1])) : big(0) for n in 1:d], abs(a[1])
    end
    acc = [big(0) for _ in 1:d]; M = big(1)
    c = Vector{UInt64}(undef, d); p = Vector{UInt64}(undef, d); image = Vector{UInt64}(undef, d)
    i = 0; accepted = 0; next = 1
    one_ = [n == 1 ? big(1) : big(0) for n in 1:d]
    while true
        i += 1
        t = _places(F, i, 0)
        q = big(t.m.p)
        for n in 1:d
            image[n] = to_mont(t.m, mod(a[n], q))
        end
        # all the inverses from one, by running products
        run = t.m.one
        for j in 1:d
            c[j] = _place_coordinates_value(t, image, j)
            p[j] = run
            run = mont_mul(t.m, run, c[j])
        end
        # A nonzero field element can vanish modulo a prime dividing its norm.
        # Never include that prime in the CRT modulus.
        iszero(run) && continue
        run = mont_inv(t.m, run)
        for j in d:-1:1
            u = mont_mul(t.m, run, p[j])
            run = mont_mul(t.m, run, c[j])
            c[j] = u
        end
        _place_coordinates!(c, p, t)
        Minv = invmod(mod(M, q), q)
        for n in 1:d
            acc[n] += M * mod((big(p[n]) - mod(acc[n], q)) * Minv, q)
        end
        M *= q
        accepted += 1
        accepted == next || continue
        next *= 2
        rec = _rational_coordinates(acc, M)
        rec === nothing && continue
        w, D = rec
        _label_mul(a, w, F) == D .* one_ && return w, D
    end
end

"Coordinates of 1/[n], kept for the level."
function _label_inverse(n::Int, F::LevelField)
    lock(LEVEL_FIELD_LOCK) do
        haskey(F.inverses, n) && return F.inverses[n]
        value = _label_inverse(_label_coordinates(n, F), F)
        F.inverses[n] = value
        _trim_level_fields!(F)
        return value
    end
end

"Coordinates of the single label [n], n ≤ h/2, on the basis."
function _label_coordinates(n::Int, F::LevelField)
    v = [big(0) for _ in 1:F.d]
    if n <= F.d
        v[n] = big(1)
    else
        for (j, c) in _basis_rows(F)[n-F.d]
            v[j] = c
        end
    end
    return v
end

"Multiply integer coordinates over a denominator by ∏[n]^e."
function _times_labels(num::Vector{BigInt}, den::BigInt, front, F::LevelField)
    for (n, e) in front
        if e > 0
            u = _label_coordinates(n, F)
            for _ in 1:e
                num = _label_mul(num, u, F)
            end
        else
            u, D = _label_inverse(n, F)
            for _ in 1:-e
                num = _label_mul(num, u, F); den *= D
            end
        end
    end
    return num, den
end

# ---- square classes ----

"The Kronecker symbol (c/a), for a > 0."
function _kronecker(c::Int, a::Int)
    out = 1
    while iseven(a)
        iseven(c) && return 0
        mod(c, 8) in (3, 5) && (out = -out)
        a ÷= 2
    end
    c = mod(c, a)
    while c != 0
        while iseven(c)
            c ÷= 2
            mod(a, 8) in (3, 5) && (out = -out)
        end
        a, c = c, a
        (mod(a, 4) == 3 && mod(c, 4) == 3) && (out = -out)
        c = mod(c, a)
    end
    return a == 1 ? out : 0
end

"""
The real quadratic subfields: `(D, G)` with G² = D for squarefree D > 1 whose discriminant c divides 2h. G is the
Gauss sum Σ_a (c/a) ζ_c^a, which is positive, written with ζ^m + ζ^(−m) = [m+1] − [m−1].
"""
function _quadratic_roots(F::LevelField)
    h = F.h; out = Tuple{Int,Vector{BigInt}}[]
    primes = _prime_factors(2h)
    for mask in 1:(1<<length(primes))-1
        D = prod(primes[i] for i in eachindex(primes) if (mask >> (i - 1)) & 1 == 1)
        c = mod(D, 4) == 1 ? D : 4D
        rem(2h, c) == 0 || continue
        sum = [big(0) for _ in 1:F.d]
        add!(n, t) = begin
            rem(n, h) == 0 && return
            m, s = _fold_label(n, h)
            sum .+= (s * t) .* _label_coordinates(m, F)
        end
        for a in 1:(c-1)÷2
            chi = _kronecker(c, a)
            iszero(chi) && continue
            m = 2h * a ÷ c
            add!(m + 1, chi)
            m >= 2 && add!(m - 1, -chi)
            m == 0 && add!(1, chi)
        end
        # the Gauss sum squares to c, and D differs from c by the square 1 or 4
        scale = c == D ? 1 : 2
        G, den = _lowest(sum, big(scale))
        isone(den) || continue
        push!(out, (D, G))
    end
    return out
end

"""
The relations of square classes in echelon form modulo 2: first over the labels from the largest down, then over the
primes dividing 2h. The second part holds the relations with even exponents, which say that √c lies in the field,
together with the quadratic subfields, whose roots are not products of q-integers.
"""
function _class_rows(F::LevelField)
    F.classes === nothing || return F.classes
    lock(LEVEL_FIELD_LOCK) do
        F.classes === nothing || return F.classes
        h = F.h; top = h ÷ 2
        gens = _relation_generators(h)
        small = _prime_factors(h)
        # A row is a product of distinct generators, so its entries and its constant stay small whatever the
        # elimination does. Parity decides everything in the first part, and needs nothing but these masks.
        parity = [falses(top) for _ in gens]
        used = [falses(length(gens)) for _ in gens]
        for (i, (_, g, _)) in enumerate(gens)
            used[i][i] = true
            for (m, c) in g
                parity[i][m] = isodd(c)
            end
        end
        function product(mask)
            e = zeros(Int, top); c = 1 // big(1)
            for i in findall(mask)
                for (m, a) in gens[i][2]
                    e[m] += a
                end
                c *= small[gens[i][3]]
            end
            return e, c
        end
        out = ClassRow[]
        for piv in top:-1:2
            i = findfirst(row -> row[piv], parity)
            i === nothing && continue
            pr = popat!(parity, i); ur = popat!(used, i)
            for j in eachindex(parity)
                parity[j][piv] || continue
                parity[j] .⊻= pr; used[j] .⊻= ur
            end
            e, c = product(ur)
            push!(out, ClassRow(piv, e, c, BigInt[]))
        end
        # What is left has even exponents throughout, so it says that √c lies in the field. So do the quadratic subfields.
        rows = Tuple{Vector{Int},Rational{BigInt},Vector{BigInt}}[(product(mask)..., BigInt[]) for mask in used]
        for (D, G) in _quadratic_roots(F)
            push!(rows, (zeros(Int, top), D // big(1), G))
        end
        times(a, b) = isempty(a) ? b : isempty(b) ? a : _label_mul(a, b, F)
        for l in _prime_factors(2h)
            has(row) = isodd(_rational_valuation(row[2], l))
            i = findfirst(has, rows)
            i === nothing && continue
            e, c, G = popat!(rows, i)
            for j in eachindex(rows)
                has(rows[j]) && (rows[j] = (rows[j][1] .+ e, rows[j][2] * c, times(rows[j][3], G)))
            end
            push!(out, ClassRow(-l, e, c, G))
        end
        F.classes = out
        _trim_level_fields!(F)
        return out
    end
end

function _rational_valuation(c::Rational{BigInt}, l::Int)
    v = 0; n = numerator(c); d = denominator(c)
    while n % l == 0
        n ÷= l; v += 1
    end
    while d % l == 0
        d ÷= l; v -= 1
    end
    return v
end

"""
The radicand ∏[n]^E · R rewritten by the relations of the level as (∏[n]^front · outside · ∏G)² · r ∏_{n ∈ rad}[n].
Returns `(front, outside, r, rad, roots)`, with `front` an exponent vector over the labels and `roots` the G.
"""
function _canonical_class(E::Vector{Int}, R::Rational{BigInt}, h::Int)
    top = h ÷ 2
    length(E) < top && (E = vcat(E, zeros(Int, top - length(E))))
    F = _level_field(h)
    roots = Vector{BigInt}[]
    for row in _class_rows(F)
        use = row.pivot > 0 ? isodd(E[row.pivot]) : isodd(_rational_valuation(R, -row.pivot))
        use || continue
        for n in 1:top
            E[n] += row.e[n]
        end
        R /= row.c
        isempty(row.root) || push!(roots, row.root)
    end
    front = [fld(a, 2) for a in E]
    rad = [n for n in 2:top if isodd(E[n])]
    outside, under = _split_square(R)
    # √(a/b) = √(ab)/b, with the root over an integer
    return front, outside / denominator(under), numerator(under) * denominator(under), rad, roots
end

# ---- values ----

_level_field(v::LevelValue) = _level_field(v.k + 2)
_same_level(a::LevelValue, b::LevelValue) = a.k == b.k ||
    throw(ArgumentError("level mismatch: $(a.k) and $(b.k) are values at different roots of unity"))

Base.:-(v::LevelValue) = LevelValue(v.k, -v.num, v.den, v.r, v.rad)
function Base.:*(c::Union{Integer,Rational}, v::LevelValue)
    iszero(c) && return zero(v)
    q = big(numerator(c)) // big(denominator(c))
    return _level_make(v.k, v.num .* numerator(q), v.den * denominator(q), v.r, v.rad)
end
Base.:*(v::LevelValue, c::Union{Integer,Rational}) = c * v
function Base.:/(v::LevelValue, c::Union{Integer,Rational})
    iszero(c) && throw(DivideError())
    return (1 // big(c)) * v
end

function Base.:*(a::LevelValue, b::LevelValue)
    _same_level(a, b)
    (iszero(a) || iszero(b)) && return zero(a)
    F = _level_field(a)
    num = _label_mul(a.num, b.num, F); den = a.den * b.den
    (isempty(a.rad) && isempty(b.rad) && isone(a.r) && isone(b.r)) && return _level_make(a.k, num, den, big(1), Int[])
    E = zeros(Int, F.h ÷ 2)
    for n in a.rad; E[n] += 1; end
    for n in b.rad; E[n] += 1; end
    front, outside, r, rad, roots = _canonical_class(E, (a.r * b.r) // big(1), F.h)
    num, den = _times_labels(num, den, _nonzero_pairs(front), F)
    for G in roots
        num = _label_mul(num, G, F)
    end
    return _level_make(a.k, num .* numerator(outside), den * denominator(outside), r, rad)
end

"1/(P√S) = (1/(P·S))·√S: the same class, and no sign to decide."
function Base.inv(v::LevelValue)
    iszero(v) && throw(DivideError())
    F = _level_field(v)
    num, den = _times_labels(v.num, big(1), [n => 1 for n in v.rad], F)
    w, D = _label_inverse(num, F)
    return _level_make(v.k, w .* v.den, D * den * v.r, v.r, v.rad)
end
Base.:/(a::LevelValue, b::LevelValue) = a * inv(b)
Base.:/(c::Union{Integer,Rational}, v::LevelValue) = c * inv(v)

function Base.:^(v::LevelValue, n::Integer)
    n < 0 && return inv(v)^(-n)
    out = one(v); base = v
    while n > 0
        isodd(n) && (out = out * base)
        n >>= 1
        n > 0 && (base = base * base)
    end
    return out
end

# ---- sums, by square class ----

"""
    LevelValueSum

A sum of values at one level with different square classes. Roots of different classes are independent over the
field unless their product is a square, which the relations of a level almost always show. Nothing rests on that:
`iszero` proves its answer either way.
"""
struct LevelValueSum
    k::Int
    terms::Dict{Tuple{BigInt,Vector{Int}},LevelValue}
end

_class_key(v::LevelValue) = (v.r, v.rad)
LevelValueSum(k::Integer) = LevelValueSum(Int(k), Dict{Tuple{BigInt,Vector{Int}},LevelValue}())
LevelValueSum(v::LevelValue) = iszero(v) ? LevelValueSum(v.k) : LevelValueSum(v.k, Dict(_class_key(v) => v))

"Two values of one class, added."
function _add_same(a::LevelValue, b::LevelValue)
    a.den == b.den && return _level_make(a.k, a.num .+ b.num, a.den, a.r, a.rad)
    g = gcd(a.den, b.den)
    ca = b.den ÷ g; cb = a.den ÷ g
    return _level_make(a.k, a.num .* ca .+ b.num .* cb, a.den * ca, a.r, a.rad)
end

_same_class(a::LevelValue, b::LevelValue) = a.r == b.r && a.rad == b.rad

function Base.:+(a::LevelValueSum, b::LevelValueSum)
    a.k == b.k || throw(ArgumentError("level mismatch: $(a.k) and $(b.k) are values at different roots of unity"))
    out = copy(a.terms)
    for (S, v) in b.terms
        if haskey(out, S)
            w = _add_same(out[S], v)
            iszero(w) ? delete!(out, S) : (out[S] = w)
        else
            out[S] = v
        end
    end
    return LevelValueSum(a.k, out)
end
Base.:-(s::LevelValueSum) = LevelValueSum(s.k, Dict(S => -v for (S, v) in s.terms))
Base.:-(a::LevelValueSum, b::LevelValueSum) = a + (-b)

function Base.:*(a::LevelValueSum, b::LevelValueSum)
    out = LevelValueSum(a.k)
    for u in values(a.terms), v in values(b.terms)
        out = out + LevelValueSum(u * v)
    end
    return out
end
Base.:*(c::Union{Integer,Rational}, s::LevelValueSum) =
    iszero(c) ? LevelValueSum(s.k) : LevelValueSum(s.k, Dict(S => c * v for (S, v) in s.terms))
Base.:*(s::LevelValueSum, c::Union{Integer,Rational}) = c * s

_collapse(s::LevelValueSum) = isempty(s.terms) ? zero(LevelValue, s.k) : length(s.terms) == 1 ? first(values(s.terms)) : s

function Base.:+(a::LevelValue, b::LevelValue)
    _same_level(a, b)
    _same_class(a, b) && return _add_same(a, b)
    return _collapse(LevelValueSum(a) + LevelValueSum(b))
end
Base.:-(a::LevelValue, b::LevelValue) = a + (-b)
Base.:+(v::LevelValue, c::Union{Integer,Rational}) = v + c * one(v)
Base.:+(c::Union{Integer,Rational}, v::LevelValue) = v + c
Base.:-(v::LevelValue, c::Union{Integer,Rational}) = v + (-c)
Base.:-(c::Union{Integer,Rational}, v::LevelValue) = (-v) + c
for op in (:+, :-, :*)
    @eval Base.$op(a::LevelValue, b::LevelValueSum) = $op(LevelValueSum(a), b)
    @eval Base.$op(a::LevelValueSum, b::LevelValue) = $op(a, LevelValueSum(b))
end

# ---- numbers, and the proof of a zero ----

"An enclosure of a level value, reusing the positive q-integers of this embedding."
function _value_interval(v::LevelValue, qn)
    total = _enclose(0)
    for (n, a) in enumerate(v.num)
        iszero(a) || (total = _iadd(total, _imul(_enclose(a), qn[n])))
    end
    root = _enclose(v.r)
    for n in v.rad
        root = _imul(root, qn[n])
    end
    return _idiv(_imul(total, _isqrt(root)), _enclose(v.den))
end

"The value in BigFloat and a certified absolute-error bound."
function _value_ball(v::LevelValue, h::Int)
    qn = [_qint_interval(n, 1, h) for n in 1:max(length(v.num), maximum(v.rad; init = 1))]
    return _interval_ball(_value_interval(v, qn))
end

function _sum_ball(s::LevelValueSum)
    h = s.k + 2
    N = maximum((max(length(v.num), maximum(v.rad; init = 1)) for v in values(s.terms)); init = 1)
    qn = [_qint_interval(n, 1, h) for n in 1:N]
    total = _enclose(0)
    for v in values(s.terms)
        total = _iadd(total, _value_interval(v, qn))
    end
    return _interval_ball(total)
end

"""
An integer L such that every nonzero sum satisfies |s| > 2^(-L).
D*s is integral, has degree at most d*2^m, and all conjugates have size at most B.
The norm gives |s| >= 1/(D*B^(degree-1)). Use |[n]| <= n and integer bounds throughout.
"""
function _sum_separation_bits(s::LevelValueSum)
    d = _totient(2(s.k + 2)) ÷ 2
    D = reduce(lcm, (v.den for v in values(s.terms)); init = big(1))
    B = big(0)
    for v in values(s.terms)
        rad = v.r * prod(big(n) for n in v.rad; init = big(1))
        root = isqrt(rad)
        root * root < rad && (root += 1)
        B += (D ÷ v.den) * sum(abs(a) * n for (n, a) in enumerate(v.num)) * root
    end
    degree = big(d) << length(s.terms)
    return Int(ndigits(D; base = 2) + (degree - 1) * ndigits(max(B, big(1)); base = 2) + 1)
end

"""
Certify distinct square classes by quadratic characters at split primes. Square-equivalent nonzero
radicands have the same character at every place. Distinct signatures therefore prove distinct classes,
whose square roots are linearly independent over the level field. Failure to separate proves nothing.
"""
function _distinct_root_classes(s::LevelValueSum; maxprimes::Int = 2)
    vals = collect(values(s.terms))
    length(vals) <= 1 && return true
    F = _level_field(s.k + 2)
    signatures = [big(0) for _ in vals]
    for i in 1:maxprimes
        t = _places(F, i, 0); m = t.m
        # A character is defined only for nonzero residues. Skip primes dividing scalar radicands.
        any(v -> iszero(mod(v.r, m.p)), vals) && continue
        for a in 1:F.d
            for (j, v) in enumerate(vals)
                u = to_mont(m, mod(v.r, m.p))
                for n in v.rad
                    u = mont_mul(m, u, _place_label(t, n, a, F.h))
                end
                character = mont_pow(m, u, (m.p - 1) ÷ 2)
                signatures[j] = (signatures[j] << 1) + (character == m.one ? 0 : 1)
            end
            length(Set(signatures)) == length(vals) && return true
        end
    end
    return false
end

"Decide zero by an enclosing interval and an integer algebraic separation bound."
function Base.iszero(s::LevelValueSum)
    isempty(s.terms) && return true
    length(s.terms) == 1 && return false
    _distinct_root_classes(s) && return false
    return _sum_iszero_enclosed(s)
end

function _sum_iszero_enclosed(s::LevelValueSum)
    limit = 0; prec = 256
    while true
        total, err = setprecision(() -> _sum_ball(s), BigFloat, prec)
        abs(total) > err && return false
        limit == 0 && (limit = _sum_separation_bits(s))
        # Precision alone is not a certificate: the entire enclosure must be below the bound.
        upper = nextfloat(abs(total) + err)
        upper < ldexp(BigFloat(1), -limit) && return true
        prec = max(2prec, limit + 64)
    end
end

function Base.:(==)(a::LevelValue, b::LevelValue)
    a.k == b.k || return false
    (iszero(a) || iszero(b)) && return iszero(a) && iszero(b)
    _same_class(a, b) && return a.den == b.den && a.num == b.num    # lowest terms, and one basis
    return iszero(LevelValueSum(a) - LevelValueSum(b))
end
Base.:(==)(a::LevelValueSum, b::LevelValueSum) = a.k == b.k && iszero(a - b)
Base.:(==)(a::LevelValue, b::LevelValueSum) = LevelValueSum(a) == b
Base.:(==)(a::LevelValueSum, b::LevelValue) = a == LevelValueSum(b)
Base.:(==)(v::Union{LevelValue,LevelValueSum}, c::Union{Integer,Rational}) = v == c * one(LevelValue, v.k)
Base.:(==)(c::Union{Integer,Rational}, v::Union{LevelValue,LevelValueSum}) = v == c
Base.:(==)(v::Union{LevelValue,LevelValueSum}, x::AbstractFloat) = isfinite(x) && v == Rational{BigInt}(x)
Base.:(==)(x::AbstractFloat, v::Union{LevelValue,LevelValueSum}) = v == x

# Equal values can have different keys and roots, so all that equal values are sure to share is their level.
Base.hash(v::Union{LevelValue,LevelValueSum}, h::UInt) = hash((:LevelValue, v.k), h)
Base.isequal(a::Union{LevelValue,LevelValueSum}, b::Union{LevelValue,LevelValueSum}) = a.k == b.k && a == b
Base.isequal(::Union{LevelValue,LevelValueSum}, ::Number) = false
Base.isequal(::Number, ::Union{LevelValue,LevelValueSum}) = false

"The number a value denotes, correct to `bits` even when its coordinates cancel strongly."
function numeric_value(x::Union{LevelValue,LevelValueSum}; bits::Integer = precision(BigFloat))
    s = x isa LevelValue ? LevelValueSum(x) : x
    iszero(s) && return setprecision(() -> BigFloat(0), BigFloat, Int(bits))
    prec = Int(bits) + 64 + maximum((ndigits(a; base = 2) for v in values(s.terms) for a in v.num); init = 0)
    while true
        total, err = setprecision(() -> _sum_ball(s), BigFloat, prec)
        if !iszero(total) && abs(total) >= err * BigFloat(2)^(Int(bits) + 2)
            return setprecision(() -> BigFloat(total), BigFloat, Int(bits))
        end
        prec *= 2
    end
end
Base.float(x::Union{LevelValue,LevelValueSum}) = Float64(numeric_value(x; bits = 128))
Base.Float64(x::Union{LevelValue,LevelValueSum}) = float(x)

# ---- display ----

function _level_value_str(v::LevelValue; most::Int = QFORM_SHOWN_TERMS)
    iszero(v) && return "0"
    terms = [(a // v.den, n == 1 ? Pair{Int,Int}[] : [n => 1]) for (n, a) in enumerate(v.num) if !iszero(a)]
    text = _q_sum_str(reverse!(terms); most = most)
    (isone(v.r) && isempty(v.rad)) && return _unwrap(text)
    return text * " √(" * _q_fraction(v.r // big(1), [n => 1 for n in v.rad]) * ")"
end

Base.show(io::IO, v::LevelValue) = print(io, "LevelValue(k = ", v.k, ", ", _level_value_str(v), ")")
Base.show(io::IO, s::LevelValueSum) =
    print(io, "LevelValueSum(k = ", s.k, ", ", join((_level_value_str(v) for v in values(s.terms)), " + "), ")")
