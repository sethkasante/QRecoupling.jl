# ---------------------------------------------------------------------
# Analytic Evaluator: Generic Complex Plane Evaluation
# Optimized for: Complex phase tracking and branch-cut stability.
# ---------------------------------------------------------------------

using LRUCache

# --- Basis Caches ---
# Key: (Type, q_sq_value)
const ANALYTIC_MAG_CACHE = LRU{Tuple{DataType, Any}, Vector}(maxsize=100)
const ANALYTIC_PHS_CACHE = LRU{Tuple{DataType, Any}, Vector}(maxsize=100)

"""
    get_complex_basis(D_max::Int, q_sq::Complex{T})
Unified accessor for complex cyclotomic basis log|Φ_d(q²)| and arg(Φ_d(q²)).
"""
function get_complex_basis(D_max::Int, q_sq::Complex{T}) where T
    key = (T, q_sq)
    if !haskey(ANALYTIC_MAG_CACHE, key) || length(ANALYTIC_MAG_CACHE[key]) < D_max
        mags, phs = build_complex_basis(D_max, q_sq)
        ANALYTIC_MAG_CACHE[key] = mags
        ANALYTIC_PHS_CACHE[key] = phs
    end
    return ANALYTIC_MAG_CACHE[key]::Vector{T}, ANALYTIC_PHS_CACHE[key]::Vector{T}
end

"""
    build_complex_basis(D_max::Int, q_sq::Complex{T})
Möbius sieve in the complex log-domain. 
"""
function build_complex_basis(D_max::Int, q_sq::Complex{T}) where T
    V_mag = Vector{T}(undef, D_max)
    V_phs = Vector{T}(undef, D_max)
    
    # Pre-calculate log(q²) and arg(q²) to avoid cumulative multiplication error
    ln_q2 = log(q_sq)
    one_q = one(q_sq)

    @inbounds for n in 1:D_max
        # xⁿ - 1 = exp(n * ln(x)) - 1
        term = exp(n * ln_q2) - one_q
        if iszero(term)
            V_mag[n] = -T(Inf)
            V_phs[n] = zero(T)
        else
            V_mag[n] = log(abs(term))
            V_phs[n] = angle(term)
        end
    end
    
    # Multiplicative Sieve: log Φ_n(x) = log(xⁿ-1) - Σ log Φ_d(x)
    @inbounds for d in 1:D_max
        m_d, p_d = V_mag[d], V_phs[d]
        isinf(m_d) && continue
        for m in (2d):d:D_max
            if !isinf(V_mag[m])
                V_mag[m] -= m_d
                V_phs[m] -= p_d
            end
        end
    end
    return V_mag, V_phs
end




# ==============================================================================
# 1. MONOMIAL ANALYTIC PROJECTOR
# ==============================================================================

function project_analytic(m::CyclotomicMonomial, q::Complex{T}) where T
    m.sign == 0 && return zero(Complex{T})
    
    # Φ_d(q²) basis
    table_m, table_p = get_complex_basis(m.max_d, q * q)

    # Magnitude and Phase from q^P
    # q^{P/2} logic:
    ln_q = log(q)
    lm = m.q_pow * real(ln_q) * 0.5
    lp = m.q_pow * imag(ln_q) * 0.5
    
    # Parity sign
    (m.sign == -1) && (lp += T(π))

    @inbounds for (d, e) in m.phi_exps
        lm += e * table_m[d]
        lp += e * table_p[d]
    end
    return exp(lm) * cis(lp)
end

# ==============================================================================
# 2. DCR ANALYTIC PROJECTOR
# ==============================================================================

function project_analytic(dcr::DCR, q::Complex{T}; ws=nothing) where T
    table_m, table_p = get_complex_basis(dcr.max_d, q * q)
    ln_q = log(q)
    l_abs_q, arg_q = real(ln_q), imag(ln_q)

    # Helper for inlined monomial evaluation
    @inline function _eval_m(mono)
        lm, lp = zero(T), zero(T)
        pe = mono.phi_exps
        @inbounds for j in eachindex(pe)
            d, e = pe[j]
            lm += e * table_m[d]
            lp += e * table_p[d]
        end
        lm += mono.q_pow * l_abs_q * 0.5
        lp += mono.q_pow * arg_q * 0.5
        (mono.sign == -1) && (lp += T(π))
        return lm, lp
    end

    # Prefactors
    lm_root, lp_root = _eval_m(dcr.root)
    lm_rad,  lp_rad  = _eval_m(dcr.radical)
    lm_base, lp_base = _eval_m(dcr.base)

    # Initial magnitude and phase
    # sqrt(radical) logic applied here
    curr_lm = lm_root + (0.5 * lm_rad) + lm_base
    curr_lp = lp_root + (0.5 * lp_rad) + lp_base
    
    isinf(curr_lm) && return zero(Complex{T})

    # Workspace
    nt = length(dcr.ratios) + 1
    _ws = (ws === nothing && T === Float64) ? _WS_C64 : ws
    if _ws === nothing || length(_ws.log_mags) < nt
        _ws = AnalyticBuffer{T}(max(nt, 512))
    end
    
    l_mags, l_phs = _ws.log_mags, _ws.phases
    max_lm = l_mags[1] = curr_lm
    l_phs[1] = curr_lp

    # HOT LOOP
    for i in 1:length(dcr.ratios)
        m = dcr.ratios[i]
        
        # Ratio update
        rlm, rlp = zero(T), zero(T)
        re = m.phi_exps
        @inbounds for j in eachindex(re)
            d, e = re[j]
            rlm += e * table_m[d]
            rlp += e * table_p[d]
        end
        rlm += m.q_pow * l_abs_q * 0.5
        rlp += m.q_pow * arg_q * 0.5
        (m.sign == -1) && (rlp += T(π))

        curr_lm += rlm
        curr_lp += rlp
        
        l_mags[i+1], l_phs[i+1] = curr_lm, curr_lp
        (curr_lm > max_lm) && (max_lm = curr_lm)
    end

    # Complex Log-Sum-Exp Summation
    acc = zero(Complex{T})
    @inbounds for i in 1:nt
        # Shift magnitude to max_lm to prevent overflow
        acc += exp(l_mags[i] - max_lm) * cis(l_phs[i])
    end

    return exp(max_lm) * acc
end




"""
    is_singular(q::Complex{T}, max_d::Int; threshold=1e-12) where T
Checks if q is dangerously close to a root of unity Φ_d(q²) = 0.
Returns (is_unsafe, d_culprit)
"""
function is_singular(q::Complex{T}, max_d::Int; threshold=1e-10) where T
    q_sq = q * q
    ln_q2 = log(q_sq)
    for d in 1:max_d
        # Root check: q^{2d} ≈ 1
        dist = abs(exp(d * ln_q2) - 1)
        if dist < threshold
            return true, d
        end
    end
    return false, 0
end
