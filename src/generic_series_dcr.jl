
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




# --- Universal DCR assembler  ---
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
    
    # 1. Prefactor Evaluation
    reset!(buf)
    pref_func(buf)
    
    if extract_radical
        m_root, m_rad = snapshot_square_root(buf)
    else
        m_root = snapshot(buf)
        m_rad  = ONE_MONOMIAL
    end
    g_max_d = max(m_root.max_d, m_rad.max_d)
    
    # 2. Base Term at z_min
    # Inject initial alternating sign if requested
    base_sign = alternating_sign ? (iseven(z_min) ? 1 : -1) : 1
    reset!(buf, base_sign)
    
    base_func(buf, z_min)
    m_base = snapshot(buf)
    g_max_d = max(g_max_d, m_base.max_d)
    
    # 3. Sequence of Dynamic Update Ratios
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
    build_dcr(pref_func, base_func, ratio_func, z_min, z_max; kwargs...)

The master DCR compiler. Generates the full combinatorial skeleton of an arbitrary 
sequence. 
"""
function build_dcr(pref_func::Function, 
                           base_func::Function, 
                           ratio_func::Function, 
                           z_min::Int, 
                           z_max::Int;
                           extract_radical::Bool=false,
                           alternating_sign::Bool=false)
    
    # Estimate initial capacity based on the summation bounds
    initial_capacity = max(20, z_max + 10)
    buf = CycloBuffer(initial_capacity)

    return _build_dcr!(buf, pref_func, base_func, ratio_func, z_min, z_max;
                              extract_radical=extract_radical, 
                              alternating_sign=alternating_sign)
end



"""
    project_dcr(dcr::DCR; k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)

Universal evaluation API for a DCR. 
- Classical Limit: Pass `q = 1`.
- Root of Unity (TQFT): Pass `k` (integer level).
- Complex Analytic: Pass `q` (complex or real parameter).
- Precision is controlled by `exact` (Float64 vs Rational/Cyclotomic).
"""
function project_dcr(dcr::DCR; 
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