# ---------------------------------------------------------------------------------
#  The closed form a user should see: Φ factors in q
#
#  A DCR is already the exact symbolic value for generic q — it is a radical times a sum of monomials
#  σ q^P ∏_d Φ_d(q²)^{e_d}. What it *shows* is the recipe; what a reader wants is the value. Summing those
#  monomials over ℤ[q] and factoring gives the closed form, and two facts make it compact:
#
#   * every Φ_d(q²) splits into irreducible cyclotomics over ℤ[q] — Φ_d(q)Φ_{2d}(q) for odd d, Φ_{2d}(q)
#     for even d — so the exponent bookkeeping is over an explicit basis of irreducibles;
#   * the **denominator is always a product of Φ's**, because it divides a product of q-integers and every
#     q-integer is a product of cyclotomics. It therefore needs no factoring at all.
#
#  Only the numerator can carry something else, and measurement says it carries at most one irreducible
#  non-cyclotomic factor (`dev/results/user_facing_exact.md` §4). So the rendering is
#
#      √(…) · (sign) q^a Φ… (remainder) / (Φ…)
#
#  which is level-independent: one expression valid at every k, and strictly more informative than a
#  cyclotomic field element. `dev/results/canonical_form.md` has the square-class theory this uses.
# ---------------------------------------------------------------------------------

const _PHIQ = Dict{Int,Any}()          # Φ_e(q), e ≥ 1, as elements of ℤ[q]

"Φ_e(q) for e ≥ 1, in the auxiliary ring ℤ[q]."
function phi_q(e::Int)
    e >= 1 || throw(DomainError(e, "cyclotomic index must be positive"))
    lock(X_LOCK) do
        get!(_PHIQ, e) do
            _, q = _qring()
            cyclotomic(e, q)
        end
    end
end

"""
Φ_d(q²) as a multiset of irreducible cyclotomic indices over ℤ[q]: `Φ_d(q)Φ_{2d}(q)` for odd `d`, and
`Φ_{2d}(q)` for even `d`. (`Φ_1(q²) = q²−1 = Φ_1Φ_2`, `Φ_2(q²) = q²+1 = Φ_4`.)
"""
_phi_sq_indices(d::Int) = isodd(d) ? (d, 2d) : (2d,)

"Exponents of the irreducible Φ_e(q) in a `CyclotomicMonomial`, with its sign and q-power."
function _mono_phi_q(m::CyclotomicMonomial)
    E = Dict{Int,Int}()
    for (d, e) in m.phi_exps
        for i in _phi_sq_indices(Int(d))
            E[i] = get(E, i, 0) + Int(e)
        end
    end
    for (i, v) in collect(E)
        iszero(v) && delete!(E, i)
    end
    return Int(m.sign), Int(m.q_pow), E
end

"""
    PhiForm

The rendered closed form of an exact symbolic value:
`sign · q^qpow · √(radical) · (Φ^num_phis · remainder) / (Φ^den_phis)`, with `content` an integer factor.
`radical` is a square class over the irreducible cyclotomics, together with its own sign and q-power, which
stay inside the root. `truncated` records that the numerator was too large to factor, in which case
`summary` carries its size instead.
"""
struct PhiForm
    basis::Symbol
    sign::Int
    qpow::Int
    rad_sign::Int
    rad_qpow::Int
    radical::Vector{Int}
    content_num::BigInt
    content_den::BigInt
    num_phis::Vector{Int}
    den_phis::Vector{Int}
    remainder::Vector{Tuple{Any,Int}}
    truncated::Bool
    summary::String
end

"""
    _dcr_ratio(dcr) -> (N, Dexp, qmin) or nothing

The deferred sum, carried out: the numerator `N ∈ ℤ[q]` over the common Φ denominator `Dexp`, with `qmin`
the q-shift that made every term a polynomial. `nothing` for a structural zero.

This is the part of the closed form that costs only the sum — no factoring — and both views below are
built on it. Separating it is what makes the cheap view cheap: at j = 8 the whole of `phi_form` is 33.7 ms
and essentially all of that is `factor(N)`.
"""
function _dcr_ratio(dcr::DCR)
    Rq, q = _qring()
    (dcr.base.sign == 0 || dcr.radical.sign == 0) && return nothing

    # --- the monomials of the sum: root·base, then fused with each ratio ---
    monos = CyclotomicMonomial[]
    buf = CycloBuffer(max(dcr.max_d, 1))
    mul!(buf, dcr.root, dcr.base)
    cur = snapshot(buf)
    push!(monos, cur)
    for r in dcr.ratios
        reset!(buf); mul!(buf, cur, r); cur = snapshot(buf)
        push!(monos, cur)
    end

    parts = [_mono_phi_q(m) for m in monos]
    # common Φ denominator, and a common q shift so every term is a polynomial
    Dexp = Dict{Int,Int}()
    for (_, _, E) in parts, (i, v) in E
        v < 0 && (Dexp[i] = max(get(Dexp, i, 0), -v))
    end
    qmin = minimum(p[2] for p in parts)

    N = Rq(0)
    for (sg, P, E) in parts
        t = Rq(sg) * q^(P - qmin)
        for (i, v) in Dexp
            w = v + get(E, i, 0)
            w > 0 && (t *= phi_q(i)^w)
        end
        for (i, v) in E
            v > 0 && !haskey(Dexp, i) && (t *= phi_q(i)^v)
        end
        N += t
    end
    return N, Dexp, qmin
end

"""
    phi_form(dcr::DCR; maxdeg = 400, basis = :q) -> PhiForm

The exact closed form of a DCR in q, as Φ factors times at most one irreducible remainder. Works entirely
from the DCR, so it is available for any symbolic value the package builds, and it is independent of the
level.

**This is the expensive expansion**, and [`xvalue`](@ref) is the cheap one. Both carry out the same
deferred sum; `phi_form` then *factors* the numerator over ℤ[q], and that factorisation is essentially the
whole cost — measured on `{6 6 6; 6 6 6}`, 10.1 ms in total of which the sum is 0.07 ms. Against
`xvalue`'s 0.26 ms on the same symbol that is **38×**, and the two reasons are the two differences: `q`
has twice the degree of `x` (the numerator here is degree 228 where the x-form is 82), and `xvalue` never
factors anything.

`maxdeg` caps the degree at which the numerator is factored; above it the numerator is reported by size
rather than expanded, because a degree-800 polynomial is not a closed form anyone reads. **A truncated
call is not a cheaper closed form, it is a different answer**, and it is why `phi_form` can look faster
than `xvalue` on large symbols: at `{9 9 9; 9 9 9}` the numerator reaches degree 504, the default
`maxdeg = 400` skips the factorisation, and the call returns in 0.40 ms against `xvalue`'s 3.0 ms. Ask for
the factors it declined (`maxdeg = 3000`) and the same call takes 52.8 ms. `truncated` records which
happened.

`basis = :x` renders the irreducible remainder in `x = q + q⁻¹`, which halves *its* degree. This is not
what `xvalue` returns: the Φ factors stay in `q` either way and only the remainder moves, whereas `xvalue`
puts the whole value in `x` over a ψ radical. Everything else is unchanged and the display says what `x`
is; the default `:q` needs no such explanation, so it is the default.
"""
function phi_form(dcr::DCR; maxdeg::Int = 400, basis::Symbol = :q)
    basis in (:q, :x) || throw(ArgumentError("basis must be :q or :x, got :$basis"))
    Rq, q = _qring()
    parts0 = _dcr_ratio(dcr)
    if parts0 === nothing
        return PhiForm(basis, 0, 0, 1, 0, Int[], big(0), big(1), Int[], Int[], Tuple{Any,Int}[], false, "")
    end
    N, Dexp, qmin = parts0

    # --- net Φ exponents: numerator from N's own factors, denominator from Dexp ---
    net = Dict{Int,Int}()
    for (i, v) in Dexp
        net[i] = get(net, i, 0) - v
    end

    # --- the radical: split its square class off, the square part joins the rational factor ---
    rsg, rP, rE = _mono_phi_q(dcr.radical)
    radical = Int[]
    for i in sort!(collect(keys(rE)))
        v = rE[i]
        r = mod(v, 2); f = (v - r) ÷ 2
        r == 1 && push!(radical, i)
        f == 0 || (net[i] = get(net, i, 0) + f)
    end
    rad_qpow = isempty(radical) && rsg == 1 ? 0 : rP
    rad_sign = isempty(radical) && rP == 0 ? 1 : rsg
    if isempty(radical) && rP == 0 && rsg == 1
        rad_sign = 1
    end

    # --- content, a bare q power, and the sign of N ---
    sgn = 1
    cN = big(1)
    if iszero(N)
        return PhiForm(basis, 0, 0, 1, 0, Int[], big(0), big(1), Int[], Int[], Tuple{Any,Int}[], false, "")
    end
    c = content(N)
    N = divexact(N, c)
    cN = BigInt(c)
    if cN < 0
        cN = -cN; sgn = -sgn
    end
    if leading_coefficient(N) < 0
        N = -N; sgn = -sgn
    end
    sh = 0
    while sh <= degree(N) && iszero(coeff(N, sh))
        sh += 1
    end
    sh > 0 && (N = shift_right(N, sh))
    qtot = qmin + sh

    # --- factor the numerator and absorb any cyclotomic factors ---
    rest = Tuple{Any,Int}[]
    trunc = false
    summ = ""
    if degree(N) == 0
        # a bare unit: nothing to factor
    elseif degree(N) <= maxdeg
        fa = factor(N)
        emax = 4 * max(dcr.max_d, 1) + 4
        for (f, m) in fa
            hit = nothing
            for e in 1:emax
                if phi_q(e) == f
                    hit = e; break
                end
            end
            hit === nothing ? push!(rest, (f, m)) : (net[hit] = get(net, hit, 0) + m)
        end
    else
        trunc = true
        cs = [coeff(N, i) for i in 0:degree(N)]
        bits = maximum(c2 -> ndigits(BigInt(c2), base = 2), cs; init = 0)
        summ = "an integer polynomial in q of degree $(degree(N)), " *
               "$(count(!iszero, cs)) nonzero coefficients, ≤$bits bits"
    end

    numφ = Int[]; denφ = Int[]
    for i in sort!(collect(keys(net)))
        v = net[i]
        v > 0 ? append!(numφ, fill(i, v)) : v < 0 && append!(denφ, fill(i, -v))
    end
    if basis === :x && !trunc && !isempty(rest)
        # fold_palindromic returns g with g(x) = q^{-D/2} f(q), so replacing f by g owes a q^{D/2}
        # per factor. Forgetting that is a silent q-power error, which is what the round-trip test caught.
        folded = Tuple{Any,Int}[]
        owed = 0
        ok = true
        for (f, m) in rest
            g = try
                fold_palindromic(f)
            catch
                nothing
            end
            if g === nothing
                ok = false
                break
            end
            owed += m * (degree(f) ÷ 2)
            push!(folded, (g, m))
        end
        if ok
            rest = folded
            qtot += owed
        else
            basis = :q
        end
    end
    return PhiForm(basis, sgn * Int(dcr.base.sign == 0 ? 0 : 1), qtot, rad_sign, rad_qpow, radical,
                   cN, big(1), numφ, denφ, rest, trunc, summ)
end

# ---- rendering ----

_qpow_str(a::Int) = a == 0 ? "" : a == 1 ? "q" : "q" * to_superscript(a)

"""
A polynomial rendered the way a reader writes one: superscripts, implicit multiplication, and the sign
carried by the operator rather than by the coefficient.
"""
function _poly_str(p, var::String)
    d = degree(p)
    d < 0 && return "0"
    out = IOBuffer()
    first = true
    for i in d:-1:0
        c = coeff(p, i)
        iszero(c) && continue
        neg = c < 0
        a = abs(c)
        if first
            neg && print(out, "−")
            first = false
        else
            print(out, neg ? " − " : " + ")
        end
        (isone(a) && i > 0) || print(out, a)
        i == 0 && continue
        print(out, var)
        i == 1 || print(out, to_superscript(i))
    end
    return String(take!(out))
end

function _phi_prod_str(ps::Vector{Int})
    isempty(ps) && return ""
    out = String[]
    i = 1
    while i <= length(ps)
        j = i
        while j < length(ps) && ps[j+1] == ps[i]; j += 1; end
        m = j - i + 1
        push!(out, "Φ" * to_subscript(ps[i]) * (m == 1 ? "" : to_superscript(m)))
        i = j + 1
    end
    return join(out, " ")
end

function Base.show(io::IO, f::PhiForm)
    if f.sign == 0
        print(io, "0")
        return
    end
    pieces = String[]
    f.sign < 0 && push!(pieces, "−")
    if !isempty(f.radical)
        inner = strip(join(filter(!isempty,
            [f.rad_sign < 0 ? "−" : "", _qpow_str(f.rad_qpow), _phi_prod_str(f.radical)]), " "))
        push!(pieces, "√(" * inner * ")")
    end
    num = String[]
    f.content_num == 1 || push!(num, string(f.content_num))
    s = _qpow_str(f.qpow); isempty(s) || push!(num, s)
    s = _phi_prod_str(f.num_phis); isempty(s) || push!(num, s)
    if f.truncated
        push!(num, "(" * f.summary * ")")
    else
        v = f.basis === :x ? "x" : "q"
        for (p, m) in f.remainder
            push!(num, "(" * _poly_str(p, v) * ")" * (m == 1 ? "" : to_superscript(m)))
        end
    end
    isempty(num) && push!(num, "1")
    push!(pieces, join(num, " "))
    body = join(filter(!isempty, pieces), " ")
    dstr = _phi_prod_str(f.den_phis)
    print(io, isempty(dstr) ? body : body * " / (" * dstr * ")")
end

"Does the value split completely into cyclotomic factors, with no irreducible remainder?"
splits_completely(f::PhiForm) = !f.truncated && isempty(f.remainder)

# Displaying a DCR never expands, sums or factors it.
Base.show(io::IO, ::MIME"text/plain", dcr::DCR) = show(io, dcr)

"""
    evaluate_phi_form(f::PhiForm, q) -> number

Evaluate a rendered closed form at a numeric `q`. This exists so that what is *displayed* can be checked
against what is *computed*: the two go through different code, and the test suite compares them. The square
root uses the principal branch, which for real `q > 0` is the package's convention; at complex `q` see
`dev/results/user_facing_exact.md` §2 on why a per-factor root is the branch-consistent choice.
"""
function evaluate_phi_form(f::PhiForm, qv::Number)
    f.sign == 0 && return zero(qv) * 0
    T = typeof(one(qv) / one(qv))
    acc = T(f.sign)
    acc *= T(f.content_num) / T(f.content_den)
    acc *= qv^f.qpow
    ev(p, v) = begin
        s = zero(T)
        for i in degree(p):-1:0
            s = s * v + T(Rational(coeff(p, i)))
        end
        s
    end
    for e in f.num_phis
        acc *= ev(phi_q(e), qv)
    end
    for e in f.den_phis
        acc /= ev(phi_q(e), qv)
    end
    v = f.basis === :x ? qv + inv(qv) : qv
    for (p, m) in f.remainder
        acc *= ev(p, v)^m
    end
    if !isempty(f.radical) || f.rad_qpow != 0 || f.rad_sign < 0
        r = T(f.rad_sign) * qv^f.rad_qpow
        for e in f.radical
            r *= ev(phi_q(e), qv)
        end
        acc *= sqrt(r)
    end
    return acc
end
