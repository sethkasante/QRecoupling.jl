# ----------------------------------------------------------------------------------
#     --- Project DCR to discrete level k ----
#  Evaluate at SU(2)_k roots of unity, q = exp(iπ/(k+2)).
#  Magnitudes are summed in log space (Log-Sum-Exp); signs and phases are exact integers
#  (see cyclotomic_values.jl), and vanishing at Φ_{k+2} is tracked as a valuation.
# ----------------------------------------------------------------------------------


const ROU_TABLE_CACHE = Dict{Tuple{DataType, Int, Int}, Any}()
const ROU_TABLE_LOCK = ReentrantLock()
const ROU_TABLE_MAXSIZE = 7_000

"""
    get_rou_table(k::Int, max_d::Int, ::Type{T})
Cached `RootOfUnityTable` for q = exp(iπ/(k+2)) covering d = 1..max_d. Thread-safe.
"""
function get_rou_table(k::Int, max_d::Int, ::Type{T}) where {T}
    key = (T, k, T === BigFloat ? precision(BigFloat) : 0)
    tab = @lock ROU_TABLE_LOCK begin
        cached = get(ROU_TABLE_CACHE, key, nothing)
        if cached === nothing || length(cached.logmag) < max_d
            length(ROU_TABLE_CACHE) >= ROU_TABLE_MAXSIZE && empty!(ROU_TABLE_CACHE)
            old = cached === nothing ? 0 : length(cached.logmag)
            cached = build_rou_table(max(max_d, 2old), k, T)
            ROU_TABLE_CACHE[key] = cached
        end
        cached
    end
    return tab::RootOfUnityTable{T}
end


#  --- Projection to discrete level k ----

"""
    project_discrete(m::CyclotomicMonomial, k::Int, ::Type{T}=Float64)
Value of a monomial at q = exp(iπ/(k+2)). Returns a real `T` when the value is real and a
`Complex{T}` otherwise. Throws a `DomainError` at a pole (negative exponent of Φ_{k+2}).
"""
function project_discrete(m::CyclotomicMonomial, k::Int, ::Type{T}=Float64) where {T}
    m.sign == 0 && return zero(T)
    tab = get_rou_table(k, m.max_d, T)
    lm, s, w, v = mono_at_root(m, tab)
    v > 0 && return zero(T)
    v < 0 && throw(DomainError(k, "Topological pole at level k=$k."))
    return apply_phase(s * exp(lm), w, tab.h)
end


"""
Log-Sum-Exp over the terms of a DCR relative to its first term. Terms whose Φ_h valuation
(`v2`, doubled) is positive vanish; a negative valuation is a pole. With a real `S`, returns
`ok = false` as soon as a contributing term has a phase other than 0 or π.
Returns (max_log, scaled_sum, ok).
"""
function _dcr_lse(::Type{S}, ratios::Vector{CyclotomicMonomial}, tab::RootOfUnityTable{T},
                  l0::T, s0::Int, v2::Int) where {S, T}
    h = tab.h
    max_l = typemin(T)
    acc = zero(S)
    cur_l, cur_s, cur_w, cur_v2 = l0, s0, 0, v2
    n = length(ratios)
    i = 0
    @inbounds while true
        if cur_v2 == 0
            if S <: Real
                cur_w % (2h) == 0 || return max_l, acc, false
                u = S(cur_w == 0 ? cur_s : -cur_s)
            else
                u = cur_s * cispi(real(S)(cur_w) / (2h))
            end
            if cur_l > max_l
                acc = acc * exp(max_l - cur_l) + u
                max_l = cur_l
            else
                acc += u * exp(cur_l - max_l)
            end
        elseif cur_v2 < 0
            throw(DomainError(h - 2, "Topological pole at level k=$(h - 2)."))
        end
        i == n && break
        i += 1
        rl, rs, rw, rv = mono_at_root(ratios[i], tab)
        cur_l += rl
        cur_s *= rs
        cur_w += rw
        cur_w >= 4h && (cur_w -= 4h)
        cur_v2 += 2rv
    end
    return max_l, acc, true
end


"""
    project_discrete(res::DCR, k::Int, ::Type{T}=Float64)
Value of a DCR at q = exp(iπ/(k+2)) by single-pass Log-Sum-Exp summation. Signs and phases are
exact; terms vanishing at Φ_{k+2} are skipped by valuation. The radical is square-rooted on the
balanced branch q^{P/2} √(ΠΨ_d). Returns a real `T` when the value is real.
"""
function project_discrete(res::DCR, k::Int, ::Type{T}=Float64) where {T}
    (res.base.sign == 0 || res.root.sign == 0 || res.radical.sign == 0) && return zero(T)
    tab = get_rou_table(k, res.max_d, T)
    h = tab.h
    isodd(_phi_exponent(res.radical, 1)) &&
        throw(ArgumentError("a radical containing Φ₁ to an odd power has no real balanced square root"))

    lr, sr, wr, vr = mono_at_root(res.root, tab)
    lq, _, wrad, vq = mono_at_root(res.radical, tab)
    lb, sb, wb, vb = mono_at_root(res.base, tab)

    # √radical = q^{P/2} √(ΠΨ); ΠΨ has phase wrad - 2P, and a negative ΠΨ contributes a factor i
    P = balanced_phase(res.radical, tab)
    wq = P + (mod(wrad - 2P, 4h) == 2h ? h : 0)

    l0 = lr + lq / 2 + lb
    s0 = sr * sb
    v2 = vq + 2(vr + vb)
    max_l, acc, ok = _dcr_lse(T, res.ratios, tab, l0, s0, v2)
    if !ok
        max_l, acc, _ = _dcr_lse(Complex{T}, res.ratios, tab, l0, s0, v2)
    end
    return apply_phase(exp(max_l) * acc, wr + wq + wb, h)
end
