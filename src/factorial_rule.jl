
# ---------------------------------------------------------------------------------
#  Factorial sums: the rule of a terminating q-hypergeometric sum
#
#  value = σ₀ · (√) Π_n [n]!^{c_n} · Σ_{z = zlo}^{zhi} (±1)^z Π_i [a_i z + b_i]!^{c_i}
#
#  A DCR expands such a sum into one cyclotomic monomial per term ratio. The rule itself is a few
#  integers (16 prefactor and 8 term factorials for a 6j symbol, at any spin), and it answers the
#  questions asked before evaluating anything:
#
#  * Valuations. The Φ_h exponent of [n]! is ⌊n/h⌋ (Legendre), so a term's valuation is
#    V(z) = Σ_i c_i ⌊(a_i z + b_i)/h⌋, piecewise constant between points where an argument crosses a
#    multiple of h. Poles, vanishing terms and contributing z-ranges follow without expanding.
#    Convention (as for DCR projections): with the prefactor under a square root, the doubled valuation
#    of a term is v_pre + 2V(z); negative is a pole, positive means the term vanishes.
#  * Cancellation zeros. The contributing sum vanishes in ℚ(ζ) iff it vanishes at every Galois
#    conjugate; testing two conjugates modulo one prime errs with probability about p⁻².
#  * Numbers. A ratio loop over q-integer tables on the contributing ranges only.
# ---------------------------------------------------------------------------------


"The affine factorial [a z + b]!^c."
struct AffineFactorial
    a::Int32
    b::Int32
    c::Int32
end

@inline _arg(f::AffineFactorial, z::Int) = Int(f.a) * z + Int(f.b)

"""
    FactorialSum

`sign0 · (√)Π pre · Σ_{z=zlo}^{zhi} (−1)^{z·alternating} Π fac(z)`, where `pre` holds pairs n => c meaning
[n]!^c and sits under a square root when `sqrt_pre` is set.

`pre` and `fac` are tuples for the fixed-shape symbols (6j, 3j, F, G), so a rule is a plain bits value:
building one allocates nothing and the kernels are compiled for its exact shape. Any iterable collection
works for other shapes.
"""
struct FactorialSum{P,F}
    pre::P
    sqrt_pre::Bool
    fac::F
    alternating::Bool
    sign0::Int8
    zlo::Int
    zhi::Int
end

const EMPTY_FACTORIAL_SUM = FactorialSum((), false, (), false, Int8(0), 0, -1)

is_empty_sum(s::FactorialSum) = s.sign0 == 0 || s.zlo > s.zhi

"q-integer factors in [a(z+1)+b]!^c / [az+b]!^c: indices lo:hi, exponent p."
@inline function _factor_step(f::AffineFactorial, z::Int)
    n = _arg(f, z)
    a, c = Int(f.a), Int(f.c)
    return a > 0 ? (n + 1, n + a, c) : (n + a + 1, n, -c)
end

"Number of q-integer multiplications in an adjacent-term ratio."
@inline _ratio_cost(s::FactorialSum) = sum(f -> abs(Int(f.a)) * abs(Int(f.c)), s.fac; init=0)

"""
    FactorialSum(range; factors=(), prefactor=(), sqrt_prefactor=false,
                 alternating=false, sign=1)

Construct a finite factorial-product series
`sign * (sqrt or identity)(prod([n]!^c for n=>c in prefactor)) *
sum((-1)^(z*alternating) * prod([a*z+b]!^c for (a,b,c) in factors))`.

Factors may be `AffineFactorial(a,b,c)` objects or integer triples. Slopes and powers
may be positive, negative or zero. Factorial arguments must be nonnegative throughout
the range. Repeated factors are combined and zero powers removed. Use `qeval` to evaluate.
This represents finite factorial-product sums, not arbitrary functions or infinite series.
The positional constructor is the unchecked internal representation.
"""
function FactorialSum(r::UnitRange{<:Integer}; factors=(), prefactor=(),
                      sqrt_prefactor::Bool=false, alternating::Bool=false, sign::Integer=1)
    sign in (-1, 0, 1) || throw(ArgumentError("sign must be -1, 0 or 1"))
    lo, hi = Int(first(r)), Int(last(r))
    pre = Dict{Int,Int}()
    fac = Dict{Tuple{Int,Int},Int}()
    for (n, c) in prefactor
        n isa Integer && c isa Integer || throw(ArgumentError("prefactor entries must be integer pairs"))
        nn, cc = Int(n), Int(c)
        pre[nn] = Base.checked_add(get(pre, nn, 0), cc)
    end
    for f in factors
        a, b, c = f isa AffineFactorial ? (f.a, f.b, f.c) : f
        all(x -> x isa Integer, (a,b,c)) || throw(ArgumentError("factor entries must be integers"))
        key = (Int(a), Int(b))
        fac[key] = Base.checked_add(get(fac, key, 0), Int(c))
    end
    p = Pair{Int32,Int32}[Int32(n)=>Int32(c) for (n,c) in sort!(collect(pre)) if c != 0 && n != 0 && n != 1]
    f = AffineFactorial[AffineFactorial(a,b,c) for ((a,b),c) in sort!(collect(fac))
                       if c != 0 && !(a == 0 && b in (0,1))]
    s = FactorialSum(Tuple(p), sqrt_prefactor, Tuple(f), alternating, Int8(sign), lo, hi)
    _validate_rule(s)
    return s
end

function _validate_rule(s::FactorialSum)
    s.sign0 in (-1,0,1) || throw(ArgumentError("rule sign must be -1, 0 or 1"))
    is_empty_sum(s) && return s
    for (n,c) in s.pre
        c == 0 && continue
        n >= 0 || throw(DomainError(n, "factorial arguments must be nonnegative"))
    end
    for f in s.fac
        f.c == 0 && continue
        for z in (s.zlo,s.zhi)
            n = Base.checked_add(Base.checked_mul(Int(f.a),z),Int(f.b))
            n >= 0 || throw(DomainError(n, "factorial arguments must be nonnegative throughout the range"))
        end
    end
    return s
end


# ---- recoupling symbols as factorial rules (doubled labels) ----

@inline function _sort4(a, b, c, d)
    a, b = minmax(a, b); c, d = minmax(c, d)
    a, c = minmax(a, c); b, d = minmax(b, d)
    b, c = minmax(b, c)
    return a, b, c, d
end

@inline function _sort3(a, b, c)
    a, b = minmax(a, b); b, c = minmax(b, c); a, b = minmax(a, b)
    return a, b, c
end

"Sorted triangle sums α and quadrilateral sums β of a 6j symbol (doubled labels); they fix its value."
@inline function racah_sums(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    α = _sort4((J1 + J2 + J3) ÷ 2, (J1 + J5 + J6) ÷ 2, (J2 + J4 + J6) ÷ 2, (J3 + J4 + J5) ÷ 2)
    β = _sort3((J1 + J2 + J4 + J5) ÷ 2, (J1 + J3 + J4 + J6) ÷ 2, (J2 + J3 + J5 + J6) ÷ 2)
    return α, β
end

"""
The 6j symbol as a factorial rule. The rule depends on the labels only through the sorted sums α and β, so
it is the same for all 24 tetrahedral relabellings and needs no canonical form. Inadmissible labels give a
rule with `sign0 = 0` of the same type (an empty sum).
"""
function sixj_sum(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    ok = _δtet(J1, J2, J3, J4, J5, J6)
    α, β = racah_sums(J1, J2, J3, J4, J5, J6)
    pre = (ntuple(i -> Int32(β[(i - 1) ÷ 4 + 1] - α[(i - 1) % 4 + 1]) => Int32(1), Val(12))...,  # [β_j − α_i]!
           ntuple(i -> Int32(α[i] + 1) => Int32(-1), Val(4))...)                            # 1/[α_i + 1]!
    fac = (AffineFactorial(1, 1, 1),                                                         # [z + 1]!
           ntuple(i -> AffineFactorial(1, -α[i], -1), Val(4))...,                            # 1/[z − α_i]!
           ntuple(i -> AffineFactorial(-1, β[i], -1), Val(3))...)                            # 1/[β_j − z]!
    return FactorialSum(pre, true, fac, true, Int8(ok ? 1 : 0), α[4], ok ? β[1] : α[4] - 1)
end

"The 3j symbol as a factorial rule, including the Wigner phase (−1)^{j1−j2−m3}."
function threej_sum(J1::Int, J2::Int, J3::Int, M1::Int, M2::Int, M3::Int = -M1 - M2)
    ok = _δ(J1, J2, J3) && _mproj_ok(J1, J2, J3, M1, M2, M3)
    pre = (map(n -> Int32(n) => Int32(1),
               ((J1 + J2 - J3) ÷ 2, (J1 - J2 + J3) ÷ 2, (-J1 + J2 + J3) ÷ 2, (J1 + M1) ÷ 2, (J1 - M1) ÷ 2,
                (J2 + M2) ÷ 2, (J2 - M2) ÷ 2, (J3 + M3) ÷ 2, (J3 - M3) ÷ 2))...,
           Int32((J1 + J2 + J3) ÷ 2 + 1) => Int32(-1))
    α1 = (J3 - J2 + M1) ÷ 2; α2 = (J3 - J1 - M2) ÷ 2
    β1 = (J1 + J2 - J3) ÷ 2; β2 = (J1 - M1) ÷ 2; β3 = (J2 + M2) ÷ 2
    fac = (AffineFactorial(1, 0, -1), AffineFactorial(1, α1, -1), AffineFactorial(1, α2, -1),
           AffineFactorial(-1, β1, -1), AffineFactorial(-1, β2, -1), AffineFactorial(-1, β3, -1))
    sign0 = !ok ? Int8(0) : iseven((J1 - J2 - M3) ÷ 2) ? Int8(1) : Int8(-1)
    zlo = max(0, -α1, -α2)
    return FactorialSum(pre, true, fac, true, sign0, zlo, ok ? min(β1, β2, β3) : zlo - 1)
end

@inline _dim_pairs(J::Int) = (Int32(J + 1) => Int32(1), Int32(J) => Int32(-1))
@inline _dim_pairs(J::Int, Js::Int...) = (_dim_pairs(J)..., _dim_pairs(Js...)...)

"Multiply the square-rooted prefactor by the quantum dimensions [J+1] = [J+1]!/[J]! (doubled labels)."
with_dimensions(s::FactorialSum, Js::Int...) =
    FactorialSum((s.pre..., _dim_pairs(Js...)...), s.sqrt_pre, s.fac, s.alternating, s.sign0, s.zlo, s.zhi)

negated(s::FactorialSum) = FactorialSum(s.pre, s.sqrt_pre, s.fac, s.alternating, -s.sign0, s.zlo, s.zhi)

"F-symbol (−1)^{j1+j2+j4+j5} √([2j3+1][2j6+1]) {6j} as a factorial rule."
function fsymbol_sum(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    s = with_dimensions(sixj_sum(J1, J2, J3, J4, J5, J6), J3, J6)
    return iseven((J1 + J2 + J4 + J5) ÷ 2) ? s : negated(s)
end

"G-symbol √(Π [2j_i+1]) {6j} as a factorial rule; dimensions in sorted order, so the rule (and every bit
of the value) is the same for all relabellings."
function gsymbol_sum(J1::Int, J2::Int, J3::Int, J4::Int, J5::Int, J6::Int)
    a, b, c, d = _sort4(J1, J2, J3, J4)
    e, f = minmax(J5, J6)
    return with_dimensions(sixj_sum(J1, J2, J3, J4, J5, J6), _merge_sorted((a, b, c, d), (e, f))...)
end

"Merge two sorted tuples."
@inline _merge_sorted(x::Tuple{}, y::Tuple) = y
@inline _merge_sorted(x::Tuple, y::Tuple{}) = x
@inline _merge_sorted(x::Tuple, y::Tuple) =
    first(x) <= first(y) ? (first(x), _merge_sorted(Base.tail(x), y)...) : (first(y), _merge_sorted(x, Base.tail(y))...)


# ---- valuations at q = e^{iπ/h} ----
#
#  Two facts about the Φ_d multiplicity E_d of a term shape the whole classification.
#
#  * **Bounded multiplicities.** For the 6j the seven denominator arguments n_r satisfy Σ_r n_r = z (because
#    Σ_j β_j = Σ_i α_i), so writing n_r = d ℓ_r + ρ_r gives E_d(z) = ⌊(1 + Σ_r ρ_r)/d⌋, hence 0 ≤ E_d(z) ≤ 6
#    (carry the counting). In particular a *term* never has negative multiplicity: it either contributes 
#    (E = 0) or vanishes, to order at most 6.
#  * **No poles for admissible labels.** Level admissibility gives α_i ≤ k, so ⌊(α_i+1)/h⌋ = 0 with h = k+2
#    and the prefactor multiplicity Σ_{i,j} ⌊(β_j − α_i)/h⌋ ≥ 0. With the bound above, no contributing term
#    can have negative valuation: an admissible symbol at a level is always finite (possibly zero). The
#    `:pole` status below is therefore reachable only for labels outside the level's admissible set, where the
#    symbol has no representation in the first place. Checked in the test suite.
#  * **Segment count.** E_d changes only where one of the seven arguments wraps modulo d, so at most 7 times
#    per period, and the summation range of a 6j is short: min_j β_j − max_i α_i ≤ Σβ/3 − Σα/4 = Σα/12 ≤ k/3
#    < h. So a level 6j has at most 8 contributing segments, which is why the fixed-size `SegList{8}` path
#    below is complete rather than heuristic (measured: it is always 1 segment for admissible labels).

"⌊n/h⌋, the Φ_h exponent of [n]!; negative arguments never reach a contributing term."
@inline fld(n::Int, h::Int) = Base.fld(n, h)

"Φ_h exponent of Π [n]!^c."
prefactor_valuation(s::FactorialSum, h::Int) = sum((Int(c) * fld(Int(n), h) for (n, c) in s.pre); init = 0)

"Φ_h exponent of the term at z."
@inline function term_valuation(s::FactorialSum, z::Int, h::Int)
    v = 0
    @inbounds for f in s.fac
        v += Int(f.c) * fld(_arg(f, z), h)
    end
    return v
end

"Points z in (zlo, zhi] where some ⌊(a z + b)/h⌋ changes value."
function valuation_breakpoints(s::FactorialSum, h::Int)
    pts = Int[]
    lo, hi = s.zlo, s.zhi
    for f in s.fac
        f.c == 0 && continue
        a, b = Int(f.a), Int(f.b)
        if a == 1 || a == -1
            r = a == 1 ? mod(-b, h) : mod(b + 1, h)
            z = lo + 1 + mod(r - (lo + 1), h)
            while z <= hi
                push!(pts, z)
                z += h
            end
        elseif a != 0
            z = lo
            while z < hi
                z = _next_break(f, z, h)
                z <= hi && push!(pts, z)
            end
        end
    end
    sort!(pts)
    unique!(pts)
    return pts
end

function _merge_ranges(rs::Vector{UnitRange{Int}})
    out = UnitRange{Int}[]
    for r in rs
        if !isempty(out) && first(r) == last(out[end]) + 1
            out[end] = first(out[end]):last(r)
        else
            push!(out, r)
        end
    end
    return out
end

"""
    classify_at_level(s, k) -> (status, segments)

Status at q = e^{iπ/(k+2)}: `:empty`, `:pole` (some term has negative valuation), `:zero` (every term
vanishes) or `:finite`; `segments` are the z-ranges of contributing terms. Cancellation zeros are not
detected here; see `is_cancellation_zero`.
"""
function classify_at_level(s::FactorialSum, k::Int)
    is_empty_sum(s) && return (:empty, UnitRange{Int}[])
    h = k + 2
    vp = (s.sqrt_pre ? 1 : 2) * prefactor_valuation(s, h)
    starts = vcat(s.zlo, valuation_breakpoints(s, h))
    segs = UnitRange{Int}[]
    for i in eachindex(starts)
        z0 = starts[i]
        z1 = i < length(starts) ? starts[i+1] - 1 : s.zhi
        v2 = vp + 2 * term_valuation(s, z0, h)
        v2 < 0 && return (:pole, UnitRange{Int}[])
        v2 == 0 && push!(segs, z0:z1)
    end
    isempty(segs) && return (:zero, segs)
    return (:finite, _merge_ranges(segs))
end

"""
Do all factorial arguments over the whole range lie in 0..h-1? Then no factor carries Φ_h: there is one
contributing segment, no pole and no vanishing term, and every factorial is in the level tables. Arguments
are affine in z, so the two ends of the range decide it. This is the common case, settled without
allocating; `classify_at_level` handles the rest.
"""
@inline function _valuation_free(s::FactorialSum, h::Int)
    for (n, _) in s.pre
        0 <= n < h || return false
    end
    for f in s.fac
        a0 = _arg(f, s.zlo); a1 = _arg(f, s.zhi)
        (0 <= min(a0, a1) && max(a0, a1) < h) || return false
    end
    return true
end

"Up to N contributing z-ranges in a bits value: what `classify_at_level` returns, without allocating."
struct SegList{N}
    r::NTuple{N,UnitRange{Int}}
    n::Int
end

SegList{N}() where {N} = SegList{N}(ntuple(_ -> 0:-1, Val(N)), 0)
Base.iterate(l::SegList, i::Int = 1) = i > l.n ? nothing : (@inbounds(l.r[i]), i + 1)
Base.length(l::SegList) = l.n
Base.eltype(::Type{<:SegList}) = UnitRange{Int}
Base.isempty(l::SegList) = l.n == 0
Base.collect(l::SegList) = UnitRange{Int}[x for x in l]

"Next point after z where ⌊(a z + b)/h⌋ changes, for any nonzero integer slope."
@inline function _next_break(f::AffineFactorial, z::Int, h::Int)
    a = Int(f.a)
    (a == 0 || f.c == 0) && return typemax(Int)
    r = mod(_arg(f,z),h)
    return z + (a > 0 ? cld(h-r,a) : fld(r,-a)+1)
end

"""
    _classify_small(s, k) -> (status, SegList) or nothing

`classify_at_level` without allocation, by sweeping z from piece to piece: at each piece the valuation is
evaluated once and the next breakpoint is the nearest one over the factors. Supports arbitrary integer
slopes; returns `nothing` after 8 disjoint contributing pieces (then `classify_at_level` is used).
"""
function _classify_small(s::FactorialSum, k::Int)
    h = k + 2
    vp = (s.sqrt_pre ? 1 : 2) * prefactor_valuation(s, h)
    segs = SegList{8}()
    z = s.zlo
    while z <= s.zhi
        nxt = s.zhi + 1
        for f in s.fac
            nxt = min(nxt, _next_break(f, z, h))
        end
        v2 = vp + 2 * term_valuation(s, z, h)
        v2 < 0 && return (:pole, SegList{8}())
        if v2 == 0
            if segs.n > 0 && last(segs.r[segs.n]) == z - 1          # merge with the previous piece
                segs = SegList{8}(Base.setindex(segs.r, first(segs.r[segs.n]):(nxt - 1), segs.n), segs.n)
            else
                segs.n == 8 && return nothing
                segs = SegList{8}(Base.setindex(segs.r, z:(nxt - 1), segs.n + 1), segs.n + 1)
            end
        end
        z = nxt
    end
    return (segs.n == 0 ? :zero : :finite, segs)
end

"Do all factorials needed on the contributing segments lie in the level tables (arguments 0..k+1)?"
function _within_level_tables(s::FactorialSum, segs, k::Int)
    all(p -> 0 <= p.first <= k + 1, s.pre) || return false
    for seg in segs, f in s.fac
        (0 <= _arg(f, first(seg)) <= k + 1 && 0 <= _arg(f, last(seg)) <= k + 1) || return false
    end
    return true
end

"Closed tetrahedron in the existing normalization: the 6j rule with its full triangle prefactor."
function tetrahedron_sum(Js::Vararg{Int,6})
    s = sixj_sum(Js...)
    return FactorialSum(s.pre,false,s.fac,s.alternating,s.sign0,s.zlo,s.zhi)
end
