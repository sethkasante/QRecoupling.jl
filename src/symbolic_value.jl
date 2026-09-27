# ---------------------------------------------------------------------------------
#  What `Symbolic()` hands back
#
#  The DCR is a good deferred structure and a poor thing to show someone. It builds in about a
#  microsecond *because it does not carry out the sum* — it keeps the prefactor, the first term and the
#  update ratios as cyclotomic monomials — and that is exactly the property worth displaying, not
#  spending. The old display spent it twice over: it ran `phi_form`, which sums the terms over ℤ[q] and
#  then factors the result, and printed the internal tree underneath. At j = 8 that was **33.7 ms**
#  against 1.1 µs to build the thing.
#
#  The default view uses x = q + q⁻¹. A bounded small rule is expanded once; other rules are printed
#  as finite products/sums of F(n) = product(U_{r-1}(x/2), r=1:n). The expensive computed views remain
#  explicit: `xvalue(v)` expands and caches the rational x-form, while `phi_form(v)` factors in q.
#  The helpers for the deferred cyclotomic view also remain available internally.
# ---------------------------------------------------------------------------------

"""
    SymbolicValue

A recoupling symbol as an exact expression in `q`, with no level and no evaluation: what `Symbolic()`
returns. It carries the factorial `rule` and its `dcr`, and is built in microseconds because the sum is
deferred in both.

Displayed in `x = q + q⁻¹`: small expressions are simplified automatically; larger expressions remain
finite factorial sums in x. Printing never forces a large polynomial expansion.

Three further views, in increasing cost and all of them optional:

| | | |
|---|---|---|
| `symbolic_terms(v)` | the summands as `CyclotomicMonomial`s | exponent arithmetic |
| `xvalue(v)` | the value as `P(x)/Q(x)·√(∏ψ)`, `x = q + q⁻¹` | carries out the sum |
| `phi_form(v)` | the same value factored over the cyclotomics in `q` | and factors it |

`Exact(k)` specialises the same expression to a level.
"""
struct SymbolicValue{R<:FactorialSum}
    rule::R
    dcr::DCR
    _xcache::Base.RefValue{Union{Nothing,XValue}}
    _xlock::ReentrantLock
end

SymbolicValue(s::FactorialSum, d::DCR) =
    SymbolicValue(s, d, Ref{Union{Nothing,XValue}}(nothing), ReentrantLock())
SymbolicValue(s::FactorialSum) = SymbolicValue(s, _factorial_dcr(s))

"""
    xvalue(v::SymbolicValue) -> XValue

The exact value for generic `q` as `√(∏_{e ∈ rad} ψ_e(x)) · num(x)/den(x)`, `x = q + q⁻¹`.
Expansion is explicit and may be expensive. The result is cached on the symbolic value; callers receive
an owned copy, so mutating its polynomial coefficients cannot corrupt later requests or display.
Square roots describe the formal algebraic expression; numerical evaluation uses the rule's branch convention.
"""
function _cached_xvalue(v::SymbolicValue; expand::Bool=true)
    lock(v._xlock) do
        if v._xcache[] === nothing && expand
            v._xcache[] = generic_value(v.rule)
        end
        return v._xcache[]
    end
end
xvalue(v::SymbolicValue) = deepcopy(_cached_xvalue(v))
xvalue(s::FactorialSum) = generic_value(_validate_rule(s))

phi_form(v::SymbolicValue; kwargs...) = phi_form(v.dcr; kwargs...)
splits_completely(v::SymbolicValue) = splits_completely(phi_form(v))

# The DCR is what the projections consume; forward rather than make callers reach inside.
project_exact(v::SymbolicValue, k::Integer; kwargs...) = project_exact(v.dcr, k; kwargs...)
project_discrete(v::SymbolicValue, k::Integer, ::Type{T} = Float64) where {T} =
    project_discrete(v.dcr, k, T)
project_analytic(v::SymbolicValue, q) = project_analytic(v.dcr, q)
Base.iszero(v::SymbolicValue) = is_empty_sum(v.rule) || v.dcr.base.sign == 0

"""
Largest generic polynomial the display writes out. Tighter than the level display's limit on purpose:
without a level the degrees grow with the labels rather than with `φ(2h)/2`, so `{3 3 3; 3 3 3}` already
has a degree-22 denominator that fills a line and says nothing a size would not.
"""
const GENERIC_MAX_DEGREE = 16
const GENERIC_MAX_CHARS = 120

"""
    symbolic_terms(v::SymbolicValue) -> Vector{CyclotomicMonomial}

The summands of the deferred sum: `t₀ = base`, `tᵢ₊₁ = tᵢ·Rᵢ` over the DCR's update ratios. Monomial
multiplication is exponent addition, so this is integer work and nothing is expanded.
"""
function symbolic_terms(v::SymbolicValue)
    d = v.dcr
    d.base.sign == 0 && return CyclotomicMonomial[]
    out = CyclotomicMonomial[d.base]
    buf = CycloBuffer(max(d.max_d, 1))
    cur = d.base
    for r in d.ratios
        reset!(buf); mul!(buf, cur, r); cur = snapshot(buf)
        push!(out, cur)
    end
    return out
end

"Exponent of `Φ_d` in a monomial, zero when absent."
function _exp_of(m::CyclotomicMonomial, d::Int)
    for (i, e) in m.phi_exps
        i == d && return e
    end
    return 0
end

"""
The largest monomial dividing every term, and the terms with it removed — the coefficientwise minimum
of the q-power and of each Φ exponent. This is what turns a list of monomials into `c·(t₀ + t₁ + ⋯)`,
and it is the only "simplification" the display does, because it is the only one that costs nothing.
"""
function _factor_common(ts::Vector{CyclotomicMonomial})
    isempty(ts) && return ONE_MONOMIAL, ts
    qc = minimum(t.q_pow for t in ts)
    ds = Int[]
    for t in ts, (d, _) in t.phi_exps
        d in ds || push!(ds, d)
    end
    sort!(ds)
    ce = Pair{Int,Int}[]
    for d in ds
        m = minimum(_exp_of(t, d) for t in ts)
        iszero(m) || push!(ce, d => m)
    end
    common = CyclotomicMonomial(1, qc, ce, isempty(ce) ? 0 : maximum(first, ce))
    rest = map(ts) do t
        ex = Pair{Int,Int}[]
        for d in ds
            e = _exp_of(t, d) - _exp_of(common, d)
            iszero(e) || push!(ex, d => e)
        end
        CyclotomicMonomial(t.sign, t.q_pow - qc, ex, isempty(ex) ? 0 : maximum(first, ex))
    end
    return common, rest
end

"Product of two monomials, by adding exponents — no buffer, no allocation beyond the result."
function _mono_mul(a::CyclotomicMonomial, b::CyclotomicMonomial)
    (a.sign == 0 || b.sign == 0) && return ZERO_MONOMIAL
    ds = Int[]
    for (d, _) in a.phi_exps; push!(ds, d); end
    for (d, _) in b.phi_exps; d in ds || push!(ds, d); end
    sort!(ds)
    ex = Pair{Int,Int}[]
    for d in ds
        e = _exp_of(a, d) + _exp_of(b, d)
        iszero(e) || push!(ex, d => e)
    end
    return CyclotomicMonomial(a.sign * b.sign, a.q_pow + b.q_pow, ex,
                              isempty(ex) ? 0 : maximum(first, ex))
end

"""
A monomial written the way one writes it: positive powers inline, negative ones gathered into a single
denominator. `sign = false` drops the leading minus, which the sum prints as an operator instead.
"""
function _mono_str(m::CyclotomicMonomial; sign::Bool = true)
    m.sign == 0 && return "0"
    num = String[]
    den = String[]
    if m.q_pow > 0
        push!(num, m.q_pow == 1 ? "q" : "q" * to_superscript(m.q_pow))
    elseif m.q_pow < 0
        push!(den, m.q_pow == -1 ? "q" : "q" * to_superscript(-m.q_pow))
    end
    for (d, e) in m.phi_exps
        t = "Φ" * to_subscript(d)
        e > 0 ? push!(num, e == 1 ? t : t * to_superscript(e)) :
                push!(den, e == -1 ? t : t * to_superscript(-e))
    end
    body = isempty(num) ? "1" : join(num, " ")
    isempty(den) || (body *= "/" * (length(den) == 1 ? den[1] : "(" * join(den, " ") * ")"))
    return (sign && m.sign < 0 ? "−" : "") * body
end

"The bracketed sum `(t₀ ± t₁ ± ⋯)`, truncated in the middle when there are too many to read."
function _sum_str(ts::Vector{CyclotomicMonomial}, maxterms::Int)
    isempty(ts) && return "0"
    n = length(ts)
    idx = n <= maxterms ? collect(1:n) :
          vcat(collect(1:maxterms-2), [0], [n])          # 0 marks the ellipsis
    out = IOBuffer()
    for (i, k) in enumerate(idx)
        if k == 0
            print(out, " + ⋯ ")
            continue
        end
        t = ts[k]
        if i == 1
            t.sign < 0 && print(out, "−")
        else
            print(out, idx[i-1] == 0 ? "+ " : (t.sign < 0 ? " − " : " + "))
        end
        print(out, _mono_str(t; sign = false))
    end
    return String(take!(out))
end

"Numerator, denominator and radical of the x-form, rendered — or described when too large to read."
function _xvalue_str(xv::XValue)
    iszero(xv.num) && return "0"
    n = _xpoly_show(_toqq(xv.num); maxdeg = GENERIC_MAX_DEGREE, maxchars = GENERIC_MAX_CHARS)
    d = _xpoly_show(_toqq(xv.den); maxdeg = GENERIC_MAX_DEGREE, maxchars = GENERIC_MAX_CHARS)
    body = d == "1" ? n : (occursin(" ", n) && !startswith(n, "(") ? "(" * n * ")" : n) * " / " *
                          (occursin(" ", d) && !startswith(d, "(") ? "(" * d * ")" : d)
    isempty(xv.rad) && return body
    rad = join(["ψ" * to_subscript(e) for e in xv.rad], "")
    return "√(" * rad * ") · " * body
end

"""
How many summands the display writes before it elides the middle, and how many characters the sum may
take. Both are needed: seven terms is few, but seven terms of twenty cyclotomic factors each is not a
line anyone reads. Past the character budget the display drops to three terms, and past that it says
what the sum is instead of writing it. `symbolic_terms(v)` returns all of them either way.
"""
const SYMBOLIC_MAX_TERMS = 6
const SYMBOLIC_MAX_CHARS = 150

"A sum too wide to write, described: how many monomials and over which cyclotomic indices."
function _sum_summary(ts::Vector{CyclotomicMonomial})
    ds = Int[]
    for t in ts, (d, _) in t.phi_exps
        d in ds || push!(ds, d)
    end
    isempty(ds) && return "(a sum of $(length(ts)) powers of q)"
    return "(a sum of $(length(ts)) monomials over Φ" * to_subscript(minimum(ds)) *
           "…Φ" * to_subscript(maximum(ds)) * ")"
end

"The whole value as it is stored: `√(radical) · lead · (t₀ ± t₁ ± ⋯)`, with nothing carried out."
function _deferred_str(v::SymbolicValue; maxterms::Int = SYMBOLIC_MAX_TERMS)
    ts = symbolic_terms(v)
    isempty(ts) && return "0", 0
    common, rest = _factor_common(ts)
    lead = _mono_mul(v.dcr.root, common)
    n = length(rest)
    body = if n == 1
        _mono_str(only(rest))
    else
        b = "(" * _sum_str(rest, maxterms) * ")"
        if length(b) > SYMBOLIC_MAX_CHARS
            b = "(" * _sum_str(rest, 3) * ")"
        end
        length(b) > SYMBOLIC_MAX_CHARS ? _sum_summary(rest) : b
    end
    parts = String[]
    rad = v.dcr.radical
    (rad.sign == 1 && isempty(rad.phi_exps) && rad.q_pow == 0) ||
        push!(parts, "√(" * _mono_str(rad) * ")")
    ls = _mono_str(lead)
    ls == "1" || push!(parts, ls)
    ls == "−1" && (parts[end] = "−")
    push!(parts, body)
    return join(parts, " "), n
end

# Bound the work before expansion, not just the length of the eventual output.
function _small_x_display(s::FactorialSum)
    is_empty_sum(s) && return true
    s.zhi - s.zlo < 8 && max_argument(s) <= 12 || return false
    cost = sum(abs(Float64(c))*n*(n-1)/2 for (n,c) in s.pre; init=0.0)
    for f in s.fac
        n = max(_arg(f,s.zlo),_arg(f,s.zhi))
        cost += abs(Float64(f.c))*n*(n-1)/2
    end
    return cost <= 256
end

function _xfactor_product(items)
    num=String[]; den=String[]
    powers=Dict{String,Int}()
    for (arg,c) in items
        powers[arg]=get(powers,arg,0)+c
    end
    for arg in sort!(collect(keys(powers)))
        c=powers[arg]
        iszero(c) && continue
        arg in ("0","1") && continue
        term="F(" * arg * ")"
        abs(c)==1 || (term *= "^" * string(abs(c)))
        push!(c>0 ? num : den,term)
    end
    n=isempty(num) ? "1" : join(num," · ")
    return isempty(den) ? n : n * " / (" * join(den," · ") * ")"
end

function _affine_xarg(f::AffineFactorial)
    a,b=Int(f.a),Int(f.b)
    a==0 && return string(b)
    str=a==1 ? "z" : a == -1 ? "-z" : string(a)*"z"
    return b==0 ? str : str * (b>0 ? "+" : "") * string(b)
end

function _deferred_x_str(s::FactorialSum)
    is_empty_sum(s) && return "0"
    pre=_xfactor_product((string(n),Int(c)) for (n,c) in s.pre)
    term=_xfactor_product((_affine_xarg(f),Int(f.c)) for f in s.fac)
    length(pre)+length(term)>240 && return "$(s.zhi-s.zlo+1)-term factorial sum in x (see v.rule)"
    pre=pre=="1" ? "" : (s.sqrt_pre ? "√("*pre*")" : "("*pre*")") * " · "
    sign=s.sign0<0 ? "−" : ""
    alt=s.alternating ? "(-1)^z · " : ""
    return sign*pre*"Σ[z="*string(s.zlo)*":"*string(s.zhi)*"] "*alt*"("*term*")"
end

function Base.show(io::IO, v::SymbolicValue)
    # Compact display remains strictly deferred, including inside arrays.
    print(io,"SymbolicValue(",_deferred_x_str(v.rule),")")
end

function Base.show(io::IO, ::MIME"text/plain", v::SymbolicValue)
    println(io,"Exact value for generic q   (x = q + q⁻¹)")
    xv=_cached_xvalue(v;expand=_small_x_display(v.rule))
    if xv !== nothing
        println(io,"  = ",_xvalue_str(xv))
        isempty(xv.rad) || println(io,"  ψ_e(x) is the minimal polynomial of 2cos(2π/e)")
    else
        println(io,"  = ",_deferred_x_str(v.rule))
        println(io,"  F(n) = ∏[r=1:n] U_{r-1}(x/2), with F(0) = 1")
        println(io,"  ",v.rule.zhi-v.rule.zlo+1," terms; the sum has not been expanded")
    end
    print(io,"  `xvalue(v)` expands the x-form; `phi_form(v)` requests cyclotomic factors in q")
end

Base.show(io::IO,v::XValue) = print(io,"XValue(",_xvalue_str(v),")")
function Base.show(io::IO,::MIME"text/plain",v::XValue)
    println(io,"Exact expression in x = q + q⁻¹")
    print(io,"  = ",_xvalue_str(v))
    isempty(v.rad) || print(io,"\n  ψ_e(x) is the minimal polynomial of 2cos(2π/e)")
end

"Evaluating a symbolic value evaluates its rule, which is the accurate route; the DCR is the fallback."
qeval(v::SymbolicValue; kwargs...) = qeval(v.rule; kwargs...)
