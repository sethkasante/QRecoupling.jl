
# ---------------------------------------------------------------------
# Generic Skeletal Assembler for Deferred Cyclotomic Representations
# ---------------------------------------------------------------------



# --- q-Hypergeometric Primitives ---

"""
    add_qint!(buf::CycloBuffer, n::Int, p::Int=1)
Algebraically multiplies the buffer by ([n]_q)^p.
"""
function add_qint!(buf::CycloBuffer, n::Int, p::Int=1)
    n <= 1 && return
    buf.q_pow += p * (1 - n) 
    ensure_capacity!(buf, n)
    @inbounds for d in 2:n
        (n % d == 0) && (buf.exps[d] += p)
    end
    n > buf.max_d && (buf.max_d = n)
end

"""
    add_qfact!(buf::CycloBuffer, n::Int, p::Int=1)
Algebraically multiplies the buffer by ([n]_q!)^p.
"""
function add_qfact!(buf::CycloBuffer, n::Int, p::Int=1)
    n <= 1 && return
    buf.q_pow += p * (n * (1 - n) ÷ 2)
    ensure_capacity!(buf, n)
    @inbounds for d in 2:n
        buf.exps[d] += p * (n ÷ d)
    end
    n > buf.max_d && (buf.max_d = n)
end


"""
    add_monomial!(buf::CycloBuffer, m::CyclotomicMonomial, p::Int=1)
Algebraically multiplies the buffer by the given monomial raised to power `p`.
This acts as the bridge between explicit monomial objects and the DCR compiler.
"""
function add_monomial!(buf::CycloBuffer, m::CyclotomicMonomial, p::Int=1)
    p == 0 && return buf
    m.sign == 0 && (buf.sign = 0; return)
    
    # Flip the sign if the monomial is negative and the applied power is odd
    if m.sign == -1 && isodd(p)
        buf.sign = -buf.sign
    end
    
    buf.q_pow += m.q_pow * p
    ensure_capacity!(buf, m.max_d)
    
    @inbounds for (d, e) in m.phi_exps
        update_exps!(buf, d, e * p)
    end
end



"""
    qint(n::Int, p::Int=1)

Returns the algebraic representation of `([n]_q)^p` as a `CyclotomicMonomial`.
"""
function qint(n::Int, p::Int=1)
    n < 0 && return ZERO_MONOMIAL
    n <= 1 && return ONE_MONOMIAL
    buf = CycloBuffer(n)
    add_qint!(buf, n, p)
    return snapshot(buf)
end

"""
    qfact(n::Int, p::Int=1)

Returns the algebraic representation of `([n]_q!)^p` as a `CyclotomicMonomial`.
"""
function qfact(n::Int, p::Int=1)
    n < 0 && return ZERO_MONOMIAL
    n <= 1 && return ONE_MONOMIAL
    buf = CycloBuffer(n)
    add_qfact!(buf, n, p)
    return snapshot(buf)
end

"""
    qbinomial(n::Int, k::Int)
Returns the cyclotomic monomial for the quantum binomial coefficient `[n]_q! / ([k]_q! [n-k]_q!)`.
"""
function qbinomial(n::Int, k::Int)
    (k < 0 || k > n) && return ZERO_MONOMIAL
    (k == 0 || k == n) && return ONE_MONOMIAL
    buf = CycloBuffer(n)
    add_qfact!(buf, n, 1)
    add_qfact!(buf, k, -1)
    add_qfact!(buf, n - k, -1)
    return snapshot(buf)
end


"""
    CyclotomicMonomial(sign::Int, q_pow::Int, dense_exps::Vector{Int})

Build sparse cyclotomic monomials directly from dense exponent arrays.
"""
function CyclotomicMonomial(sign::Int, q_pow::Int, dense_exps::Vector{Int})
    max_d = length(dense_exps)
    max_d == 0 && return ONE_MONOMIAL
    
    phi_exps = Pair{Int,Int}[]
    sizehint!(phi_exps, count(!iszero, dense_exps)) 
    for (d, exp) in enumerate(dense_exps)
        exp != 0 && push!(phi_exps, d => exp)
    end
    
    return CyclotomicMonomial(sign, q_pow, phi_exps, max_d)
end


# --- utilities for DCR construction ---- 

"""
    fuse_radical(res::DCR, m::CyclotomicMonomial)
Algebraically fuses a new prefactor `m` into the `radical` (square root) 
part of a DCR. Extracts resulting perfect squares into the `root`.
"""
function fuse_radical(res::DCR, m::CyclotomicMonomial)
    res.base.sign == 0 && return res
    
    m_d = max(res.max_d, m.max_d)
    buf = CycloBuffer(m_d)

    mul!(buf, res.radical, m)
    
    # extract squares: radical * extra = (bonus_root)^2 * final_rad
    root_bonus, final_rad = snapshot_square_root(buf)
    
    # update prefactors
    final_root = res.root * root_bonus
    
    return DCR(final_root, final_rad, res.base, res.ratios, res.z_range, m_d)
end

"""
    fuse_root(res::DCR, m::CyclotomicMonomial)
Multiplies the 'root' (rational) part of the DCR by `m`.
Useful for quantum dimensions or integer phases.
"""
function fuse_root(res::DCR, m::CyclotomicMonomial)
    new_root = res.root * m
    new_max_d = max(res.max_d, m.max_d)
    return DCR(new_root, res.radical, res.base, res.ratios, new_max_d)
end


# --- Universal DCR & series assembler  ---

"""
    build_series(summand::Function, z_range::UnitRange; prefactor::CyclotomicMonomial=ONE_MONOMIAL, extract_radical::Bool=false)

Constructs a `DCR` from a user-provided summand function.
`summand(z)` must return a `CyclotomicMonomial` representing the z-th term in the series.
"""
function build_series(summand::Function, z_range::UnitRange{Int};
                      prefactor::CyclotomicMonomial = ONE_MONOMIAL,
                      extract_radical::Bool = false)
    z_min = first(z_range)
    z_max = last(z_range)

    (z_min > z_max || prefactor.sign == 0) && return ZERO_DCR

    # We allocate a single buffer to handle both the prefactor and ratio divisions
    buf = CycloBuffer(prefactor.max_d)

    # evaluate Prefactor
    if extract_radical
        reset!(buf)
        add_monomial!(buf, prefactor, 1)
        m_root, m_rad = snapshot_square_root(buf)
    else
        m_root = prefactor
        m_rad = ONE_MONOMIAL
    end

    g_max_d = max(m_root.max_d, m_rad.max_d)

    # Base Term
    m_base = summand(z_min)
    m_base.sign == 0 && return ZERO_DCR  # Trivial zero series
    g_max_d = max(g_max_d, m_base.max_d)

    # ratios
    ratios = Vector{CyclotomicMonomial}(undef, z_max - z_min)
    curr_term = m_base

    for (i, z) in enumerate(z_min:(z_max - 1))
        next_term = summand(z + 1)
        # stop and truncate sum if we hit a structural zero
        if next_term.sign == 0 
            resize!(ratios, i - 1)
            z_range = z_min:(z_min + i - 1)
            break 
        end
        
        # calculate next_term / curr_term 
        reset!(buf)
        add_monomial!(buf, next_term, 1)
        add_monomial!(buf, curr_term, -1)
        
        ratio = snapshot(buf)
        ratios[i] = ratio
        g_max_d = max(g_max_d, ratio.max_d)

        curr_term = next_term
    end

    return DCR(m_root, m_rad, m_base, ratios, z_range, g_max_d)
end


"""
    build_dcr!(pref_func, base_func, ratio_func, z_min, z_max; kwargs...)

The master DCR compiler. Generates the full combinatorial skeleton of an arbitrary 
sequence by orchestrating closures over a highly optimized `CycloBuffer`.

**Keyword Arguments:**
- `extract_radical::Bool` (default `false`): If true, searches the prefactor for perfect squares and separates them into the `root` and `radical` fields.
- `alternating_sign::Bool` (default `false`): If true, automatically injects a (-1)^z term into the sequence (standard for recoupling invariants).
"""
function build_dcr!(buf::CycloBuffer, 
                            pref_func::Function, 
                            base_func::Function, 
                            ratio_func::Function, 
                            z_min::Int, 
                            z_max::Int;
                            extract_radical::Bool=false,
                            alternating_sign::Bool=false)
    
    (z_min > z_max || buf.sign == 0) && return ZERO_DCR
    
    # prefactor evaluation
    reset!(buf)
    pref_func(buf)
    
    if extract_radical
        m_root, m_rad = snapshot_square_root(buf)
    else
        m_root = snapshot(buf)
        m_rad  = ONE_MONOMIAL
    end
    g_max_d = max(m_root.max_d, m_rad.max_d)
    
    # base term at z_min
    # inject initial alternating sign if requested
    base_sign = alternating_sign ? (iseven(z_min) ? 1 : -1) : 1
    reset!(buf, base_sign)
    
    base_func(buf, z_min)
    m_base = snapshot(buf)
    g_max_d = max(g_max_d, m_base.max_d)
    
    # sequence of update ratios
    ratios = Vector{CyclotomicMonomial}(undef, z_max - z_min)
    ratio_sign = alternating_sign ? -1 : 1
    
    for (i, z) in enumerate(z_min:(z_max - 1))
        reset!(buf, ratio_sign)
        ratio_func(buf, z)
        ratios[i] = snapshot(buf)
        g_max_d = max(g_max_d, ratios[i].max_d)
    end
    
    return DCR(m_root, m_rad, m_base, ratios, z_min:z_max, g_max_d)
end


"""
    qseries(summand::Function, z_range::UnitRange{Int}; 
            prefactor::CyclotomicMonomial = ONE_MONOMIAL, 
            extract_radical::Bool = false)

A lightweight, user-friendly API for compiling a finite q-hypergeometric series 
into a Deferred Cyclotomic Representation (DCR). 

# Example
```julia
my_series = qseries(3:10) do z
    return (-1)^z * (qfact(z) / qfact(z - 3))
end
```
"""
function qseries(summand::Function, z_range::UnitRange{Int}; kwargs...)
    return build_series(summand, z_range; kwargs...)
end


"""
    qeval(m::CyclotomicMonomial; k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)

Universal evaluation API for a single CyclotomicMonomial. 
- Classical Limit: Pass `q = 1`.
- Root of Unity (TQFT): Pass `k` (integer level).
- Complex Analytic: Pass `q` (complex or real parameter).
- Precision is controlled by `exact` (Float64 vs Rational/Cyclotomic).
"""
function qeval(m::CyclotomicMonomial; 
               k=nothing, 
               q=nothing, 
               exact::Bool=false, 
               T::Type=Float64)
    
    # ---  classical limit (q -> 1) ---
    if !isnothing(q) && (q == 1 || q == 1.0)
        return exact ? project_classical_exact(m) : project_classical(m, T)
    end

    # --- complex/analytic projection ---
    if !isnothing(q)
        return project_analytic(m, q)
    end

    # --- root of unity projection ---
    if !isnothing(k)
        return exact ? project_exact(m, k) : project_discrete(m, k, T)
    end

    throw(ArgumentError("Projection target missing. Specify `k` (integer level) or `q` (parameter)."))
end


"""
    qeval(dcr::DCR; k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)

Universal evaluation API for a DCR. 
- Classical Limit: Pass `q = 1`.
- Root of Unity (TQFT): Pass `k` (integer level).
- Complex Analytic: Pass `q` (complex or real parameter).
- Precision is controlled by `exact` (Float64 vs Rational/Cyclotomic).
"""
function qeval(dcr::DCR; 
                     k=nothing, 
                     q=nothing, 
                     exact::Bool=false, 
                     T::Type=Float64)
    
    # ---  classical limit (q -> 1) ---
    if !isnothing(q) && (q == 1 || q == 1.0)
        return exact ? project_classical_exact(dcr) : project_classical(dcr, T)
    end

    # --- complex/analytic projection ---
    if !isnothing(q)
        return project_analytic(dcr, q)
    end

    # --- root of unity projection ---
    if !isnothing(k)
        return exact ? project_exact(dcr, k) : project_discrete(dcr, k, T)
    end

    throw(ArgumentError("Projection target missing. Specify `k` (integer level) or `q` (parameter)."))
end