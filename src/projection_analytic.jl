
# ---------------------------------------------------------------------------
#                   --- Analytic Continuation ----
# Project Deferred Cyclotomic Representations (DCR) analytically for q ∈ ℂ
# Convention: [n]_q =q^{n}-q^{-n} / (q-q^{-1}) = q^{1-n} * \prod \Phi_d(q^2)
# Optimized for complex phase tracking and branch-cut stability.
# Thread-Safe
# ---------------------------------------------------------------------------



# ---- On-The-Fly: cyclotomic table builder --- --

"""
    build_analytic_table(max_d, q_sq [, q]) -> Vector

Numerical values of `Φ_d(q²)` up to `max_d`.

`(q²)^n − 1` is accumulated by the recurrence `wₙ₊₁ = q²·wₙ + w₁` rather than formed as a difference of a
power and 1. That difference loses everything as `q → 1`: `qdim(1/2; q = 1 − 1e−7)` came back with a
relative error of **8e−11**, against 1e−16 now, while every other numeric path in the package held 1e−15
there. The recurrence only ever adds quantities of the same size, so nothing cancels after the seed.

The seed `w₁ = q² − 1` is the one place cancellation can still bite, and passing `q` itself avoids it:
`(q − 1)(q + 1)` is accurate to an ulp for every `q`, because near 1 the subtraction `q − 1` is exact.
Without `q` the plain difference is used, which is what a caller that only has `q²` can do.
"""
function build_analytic_table(max_d::Int, q_sq::T, q = nothing) where T
    max_d == 0 && return Vector{T}(undef, 0)
    table = Vector{T}(undef, max_d)

    w1 = q === nothing ? q_sq - one(T) : (T(q) - one(T)) * (T(q) + one(T))
    val = w1
    @inbounds for n in 1:max_d
        table[n] = val                       # (q²)^n − 1 = ∏_{d|n} Φ_d(q²)
        val = val * q_sq + w1
    end
    # Sieve over multiples instead of testing every possible divisor. At step d,
    # table[d] is finalized; each entry is still divided in ascending divisor order.
    # This preserves the arithmetic of the old loop in O(N log N) visits.
    @inbounds for d in 1:(max_d ÷ 2)
        td = table[d]
        iszero(td) && throw(DomainError(q_sq,
            "q² is a primitive $(d)-th root of unity; Φ_$(2d) is singular here."))
        for n in (2d):d:max_d
            table[n] /= td
        end
    end
    return table
end



# -- Internal computations ---


"""
    _eval_mono_analytic(m::CyclotomicMonomial, q::T, table::Vector{T}) where T
Evaluates the monomial directly using exact multiplication.
M =  q^{q_pow} * ∏ Φ_d(q²)^e
"""
@inline function _eval_mono_analytic(m::CyclotomicMonomial, q::T, table::Vector{T}) where T
    m.sign == 0 && return zero(T)
    
    val = one(T)
    pe = m.phi_exps
    
    
    @inbounds for i in 1:length(pe)
        p = pe[i]
        d, e = p.first, p.second
        
        # unroll exponent to bypass complex `^` operator
        td   = table[d]
        
        if e == 1
            val *= td
        elseif e == -1
            iszero(td) && throw(DomainError(d, "Division by zero: Φ_$d(q²) = 0."))
            val /= td
        elseif e == 2
            val *= (td * td)
        elseif e == -2
            iszero(td) && throw(DomainError(d, "Division by zero: Φ_$d(q²) = 0."))
            val /= (td * td)
        elseif e > 0
            val *= td^e
        else
            iszero(td) && throw(DomainError(d, "Division by zero: Φ_$d(q²) = 0."))
            val /= td^(-e)
        end
    end
    
    return m.sign * (q^m.q_pow) * val
end


# -----  Core Evaluator (Single-Pass) --- 

"""
    _radical_sqrt(rad, q, table)

Square root of a square-free radical monomial on the balanced branch. Writing
`rad = σ q^P Π Ψ_d^{e_d}` with `Ψ_d(q) = q^{-φ(d)} Φ_d(q²)` and `P = balanced_phase(rad)`, this returns

    √σ · q^{P/2} · Π (√Ψ_d)^{e_d} ,

i.e. **one principal root per Ψ_d**, never a single root of the assembled product. That distinction is not
cosmetic: √ is not multiplicative across its branch cut, so a root of the product differs from the product
of roots by a sign that depends on how the individual phases add. Two sides of a coherence identity then
assemble different products and their radicals stop cancelling. Taking one root per factor makes the branch
depend only on the *set* of factors, which both sides share. Measured: Biedenharn--Elliott 40/40 on and off
the unit circle, against 9–36/40 for a root of the product.
Each `Ψ_d` is real on `|q| = 1`, and is snapped to the real axis
when its imaginary part is at rounding level so that the branch is reproducible.

For positive real `q` the ordinary square root of the assembled value is kept.
Negative-real DCR projections enter through the complex branch, as in the factorial-rule evaluator.
"""
function _radical_sqrt(rad::CyclotomicMonomial, q::T, table) where T
    T <: Real && return sqrt(_eval_mono_analytic(rad, q, table))
    P = balanced_phase(rad)
    acc = sqrt(T(rad.sign))
    tol = sqrt(eps(real(T)))
    for (d, e) in rad.phi_exps
        e == 0 && continue
        ψ = table[d] * q^(-_totient(d))
        if abs(imag(ψ)) <= tol * abs(ψ)
            ψ = T(real(ψ))
        end
        iszero(ψ) && e < 0 && throw(DomainError(d, "Division by zero: Φ_$d(q²) = 0."))
        r = sqrt(ψ)
        acc *= e > 0 ? r^e : inv(r)^(-e)
    end
    phase = _negative_real_axis(q) ? _avalue(_aqhalfpow(q,P)) : exp((P // 2) * log(q))
    return acc * phase
end

function _eval_analytic_dcr(res::DCR, q::T, q_sq::T) where T
    # for cyclotomic polynomials
    table = build_analytic_table(res.max_d, q_sq, q)

    # project prefactors
    val_root = _eval_mono_analytic(res.root, q, table)
    val_base = _eval_mono_analytic(res.base, q, table)

    # root * √(radical) * base; the radical is rooted factor by factor (see `_radical_sqrt`)
    pref_val = val_root * _radical_sqrt(res.radical, q, table) * val_base

    # This evaluator multiplies cyclotomic values directly: there is no split-exponent scaling, no
    # compensation and no error bound, so its intermediates overflow or underflow well inside the range
    # of ordinary inputs. In numerical checks, at j = 20 it returned NaN or a
    # spurious zero for five of eight sampled q, and at j = 30 for six of eight, while the factorial-rule
    # kernel is exact against a 512-bit reference throughout. Rather than hand back a silently wrong
    # number, say so and name the path that works. A *structurally* zero DCR is a different thing and
    # still returns zero.
    structural_zero = res.base.sign == 0 || res.radical.sign == 0 || res.root.sign == 0
    structural_zero && return zero(T)
    if iszero(pref_val) || !isfinite(abs(pref_val))
        throw(DomainError(q, "the DCR evaluator underflowed or overflowed at q = $q: it multiplies " *
            "cyclotomic values without scaling. Use the factorial-rule kernel for a certified value — " *
            "`q6j(labels...; q = $q)` for a symbol, or `qeval(rule; q = $q)` for a general rule. The DCR " *
            "is the symbolic view; see `phi_form` for its closed form."))
    end

    sum_val = pref_val
    curr_term = pref_val

    ratios = res.ratios
    @inbounds for i in 1:length(ratios)
        r_val = _eval_mono_analytic(ratios[i], q, table)
        curr_term *= r_val
        iszero(curr_term) && break
        sum_val += curr_term
    end

    isfinite(abs(sum_val)) || throw(DomainError(q,
        "the DCR evaluator overflowed while summing at q = $q; use the factorial-rule kernel " *
        "(`q6j(labels...; q = $q)`) for a certified value."))
    return sum_val
end


# ---  Projection of cyclotomic monomials  -----

"""
    project_analytic(m::CyclotomicMonomial, q::Number)
Evaluates a single quantum monomial
"""
function project_analytic(m::CyclotomicMonomial, q::Number)
    T = typeof(float(q)) 
    
    m.sign == 0 && return zero(T)
    # Only q = 1 itself is refused. The old guard rejected a whole `√eps` neighbourhood because the table
    # formed `(q²)^n − 1` as a difference and lost its digits there; the recurrence in
    # `build_analytic_table` does not, so `q = 1 + 1e−9` is now an ordinary point and is answered as
    # accurately as any other — and consistently with `q6j`, which always accepted it.
    q == one(T) && throw(ArgumentError("q = 1: use classical mode instead (set q=1)."))
    
    q_T = T(q)
    q_sq = q_T * q_T
    
    table = build_analytic_table(m.max_d, q_sq, q_T)
    return _eval_mono_analytic(m, q_T, table)
end

# ---  Projection of DCR -----

"""
    project_analytic(dcr::DCR, q::Number)
Fast, thread-safe for evaluating TQFT symbols at any generic q.
"""
function project_analytic(dcr::DCR, q::Number)
    _negative_real_axis(q) && (q = complex(real(q)))
    T = typeof(q * 1.0)
    
    dcr.base.sign == 0 && return zero(T)
    abs(q - 1.0) < 1e-12 && throw(ArgumentError("For q=1, use mode=:classical instead."))
    
    q_T = T(q)
    q_sq = q_T * q_T

    return _eval_analytic_dcr(dcr, q_T, q_sq)
end
