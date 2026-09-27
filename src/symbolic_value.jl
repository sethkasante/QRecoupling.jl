# ---------------------------------------------------------------------------------
#  What `Symbolic()` hands back
#
#  Construction retains only the factorial rule; compatibility DCRs are built on demand.
#
#  The display is the rule itself, always: a finite factorial sum in x, written over
#  F(n) = product(U_{r-1}(x/2), r=1:n). It never sums, never factors and never depends on the size of
#  the labels, so what you see for j = 1 and for j = 30 is the same kind of object.
#
#  Expanding is a separate, explicit request, and there are two of them because there are two bases:
#  `xvalue(v)` carries out the sum in x = q + q⁻¹, and `phi_form(v)` carries it out in q and factors the
#  result over the cyclotomics. The second is the expensive one — see its docstring.
# ---------------------------------------------------------------------------------

"""
    SymbolicValue

A recoupling symbol as an exact expression in `q`, with no level and no evaluation: what `Symbolic()`
returns. It retains the factorial `rule` without evaluating the sum. The compatibility property
`v.dcr` constructs and caches a DCR only when explicitly accessed.

**Displayed as the rule**, a finite factorial sum in `x = q + q⁻¹`, whatever the labels are. Printing
carries out no sum, so it costs the same for `j = 1` as for `j = 30` and a small symbol is not silently
treated differently from a large one.

Expanding is a separate request, and which of the two you want depends on the basis:

| | | cost at j = 6 |
|---|---|---|
| `symbolic_terms(v)` | the summands as `CyclotomicMonomial`s, unsummed | 1 µs, exponent arithmetic |
| `xvalue(v)` | the sum carried out in `x`: `√(∏ψ)·P(x)/Q(x)` | 0.3 ms |
| `phi_form(v)` | the sum carried out in `q` and **factored** over the cyclotomics | 10 ms |

`phi_form` is thirty to fifty times the cost of `xvalue` on the same symbol, and essentially all of the
difference is one call to `factor` over ℤ[q]; see its docstring. `xvalue` is cached on the value,
`phi_form` is not.

`Exact(k)` specialises the same expression to a level.
"""
struct SymbolicValue{R<:FactorialSum}
    rule::R
    _dcache::Base.RefValue{Union{Nothing,DCR}}
    _xcache::Base.RefValue{Union{Nothing,XValue}}
    _xlock::ReentrantLock
end

SymbolicValue(s::FactorialSum, d::DCR) =
    SymbolicValue(s, Ref{Union{Nothing,DCR}}(d), Ref{Union{Nothing,XValue}}(nothing), ReentrantLock())
SymbolicValue(s::FactorialSum) =
    SymbolicValue(s, Ref{Union{Nothing,DCR}}(nothing), Ref{Union{Nothing,XValue}}(nothing), ReentrantLock())

function Base.getproperty(v::SymbolicValue, name::Symbol)
    name === :dcr || return getfield(v,name)
    lock(getfield(v,:_xlock)) do
        cache=getfield(v,:_dcache)
        cache[] === nothing && (cache[]=_factorial_dcr(getfield(v,:rule)))
        return cache[]::DCR
    end
end
Base.propertynames(v::SymbolicValue, private::Bool=false) =
    private ? (:rule,:dcr,:_dcache,:_xcache,:_xlock) : (:rule,:dcr)

"""
    xvalue(v::SymbolicValue) -> XValue

The exact value for generic `q` as `√(∏_{e ∈ rad} ψ_e(x)) · num(x)/den(x)`, `x = q + q⁻¹`.

This is the **cheaper of the two expansions** — `phi_form(v)` carries out the same sum in `q` and then
factors it, which costs 30–50× as much (0.26 ms against 10.1 ms on `{6 6 6; 6 6 6}`). Use `xvalue` unless
the cyclotomic factors are what you are after.

Expansion is explicit and may still be expensive at large labels; nothing about displaying `v` triggers
it. The result is cached on the symbolic value, so asking twice is free; callers receive an owned copy, so
mutating its polynomial coefficients cannot corrupt later requests or display. Square roots describe the
formal algebraic expression; numerical evaluation uses the rule's branch convention.
"""
function _cached_xvalue(v::SymbolicValue)
    lock(v._xlock) do
        v._xcache[] === nothing && (v._xcache[] = generic_value(v.rule))
        return v._xcache[]::XValue
    end
end
xvalue(v::SymbolicValue) = deepcopy(_cached_xvalue(v))
xvalue(s::FactorialSum) = generic_value(_validate_rule(s))

phi_form(v::SymbolicValue; kwargs...) = phi_form(v.dcr; kwargs...)
splits_completely(v::SymbolicValue) = splits_completely(phi_form(v))

# Numerical projections use the rule directly; the legacy exact projector consumes a DCR.
project_exact(v::SymbolicValue, k::Integer; kwargs...) = project_exact(v.dcr, k; kwargs...)
project_discrete(v::SymbolicValue, k::Integer, ::Type{T} = Float64) where {T} =
    qeval(v.rule;k=k,T=T)
project_analytic(v::SymbolicValue, q) = qeval(v.rule;q=q)
Base.iszero(v::SymbolicValue) = is_empty_sum(v.rule)

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

"""
The rule as it stands, in `x`: a prefactor over the ψ basis times `Σ_z (−1)^z ∏ F(a z + b)^c`.

Nothing is summed, so this costs the same for two terms and for two hundred, and it is what `Symbolic()`
shows for *every* symbol. It used to expand small rules automatically, which made `{1 1 1; 1 1 1}` and
`{8 8 8; 8 8 8}` print different kinds of object and hid where the cost of expansion begins; `xvalue(v)`
asks for that explicitly now.

A one-term sum is written without the `Σ`, since `Σ[z=0:0]` in front of a single product is noise —
`qdim` and `theta_value` are the rules that hit it.
"""
function _deferred_x_str(s::FactorialSum)
    is_empty_sum(s) && return "0"
    pre = if s.sqrt_pre
        half,rad=_halve_exponents(psi_exponents(s.pre))
        num=String[]; den=String[]
        for e in sort!(collect(keys(half)))
            c=half[e]; t="ψ"*to_subscript(e)*(abs(c)==1 ? "" : to_superscript(abs(c)))
            push!(c>0 ? num : den,t)
        end
        r=isempty(rad) ? "" : "√("*join(("ψ"*to_subscript(e) for e in rad),"")*")"
        isempty(r) || pushfirst!(num,r)
        p=isempty(num) ? "1" : join(num,"")
        isempty(den) ? p : p*" / ("*join(den,"")*")"
    else
        _xfactor_product((string(n),Int(c)) for (n,c) in s.pre)
    end
    single = s.zlo == s.zhi
    term = single ? _xfactor_product((string(_arg(f,s.zlo)),Int(f.c)) for f in s.fac) :
                    _xfactor_product((_affine_xarg(f),Int(f.c)) for f in s.fac)
    length(pre)+length(term)>240 && return "$(s.zhi-s.zlo+1)-term factorial sum in x (see `v.rule`)"
    neg = s.sign0<0
    single && s.alternating && isodd(s.zlo) && (neg = !neg)
    sign = neg ? "" : ""
    body = if single
        term
    else
        "Σ[z="*string(s.zlo)*":"*string(s.zhi)*"] "*(s.alternating ? "(-1)^z · " : "")*"("*term*")"
    end
    # nothing left but the prefactor: write it bare rather than wrapping it for a product it is not in
    body=="1" && return sign*(pre=="1" ? "1" : pre)
    pre=="1" && return sign*body
    return sign*(occursin(" / ",pre) ? "("*pre*")" : pre)*" · "*body
end

function Base.show(io::IO, v::SymbolicValue)
    # Compact display remains strictly deferred, including inside arrays.
    print(io,"SymbolicValue(",_deferred_x_str(v.rule),")")
end

function Base.show(io::IO, ::MIME"text/plain", v::SymbolicValue)
    println(io,"Exact value for generic q   (x = q + q⁻¹)")
    print(io,"  = ",_deferred_x_str(v.rule))
    if is_empty_sum(v.rule)
        return                                   # an inadmissible symbol needs no legend
    end
    println(io)
    println(io,"  F(n) = ∏[r=1:n] U_{r-1}(x/2), with F(0) = 1")
    v.rule.sqrt_pre && println(io,"  ψ_e(x) is the minimal polynomial of 2cos(2π/e)")
    n = v.rule.zhi-v.rule.zlo+1
    println(io,"  ",n,n == 1 ? " term" : " terms","; the sum has not been expanded")
    print(io,"  `xvalue(v)` carries out the sum in x; `phi_form(v)` carries it out in q and factors it")
end

Base.show(io::IO,v::XValue) = print(io,"XValue(",_xvalue_str(v),")")
function Base.show(io::IO,::MIME"text/plain",v::XValue)
    println(io,"Exact expression in x = q + q⁻¹")
    print(io,"  = ",_xvalue_str(v))
    isempty(v.rad) || print(io,"\n  ψ_e(x) is the minimal polynomial of 2cos(2π/e)")
end

"Evaluating a symbolic value evaluates its rule, which is the accurate route; the DCR is the fallback."
qeval(v::SymbolicValue; kwargs...) = qeval(v.rule; kwargs...)
