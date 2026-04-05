
#assembler_dcr.jl
# ---------------------------------------------------------------------
# High-level builders for DCR Symbols 
# ---------------------------------------------------------------------


#zero DCR 
const ZERO_DCR = DCR(ZERO_MONOMIAL, ONE_MONOMIAL, ZERO_MONOMIAL, CyclotomicMonomial[], 0:-1, 0)


# --- utilities for DCR construction ---- 


"""
    fuse_radical(res::DCR, m::CyclotomicMonomial)
Algebraically fuses a new prefactor `m` into the `radical` (square root) 
part of a DCR. 

**Notes:**
1. Multiplies the current radical by the new monomial.
2. Extracts perfect squares using `snapshot_square_root`.
3. Moves extracted squares into the `root` (rational) prefactor.
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
Use this for quantum dimensions or integer phases.
"""
function fuse_root(res::DCR, m::CyclotomicMonomial)
    new_root = mul(res.root, m)
    new_max_d = max(res.max_d, m.max_d)
    return DCR(new_root, res.radical, res.base, res.ratios, new_max_d)
end




# --- DCR assembler  ---

"""
    build_dcr!(buf, pre_func, base_func, ratio_func, z_min, z_max)
Generic assembler for symbols defined by a q-hypergeometric series:
    Symbol = Prefactor * Σ_{z} [Base(z) * Π Ratio(i)]

**Arguments:**
- `pre_func`: Populates `buf` with terms for the square-root prefactor.
- `base_func`: Populates `buf` with the starting term of the sum at `z_min`.
- `ratio_func`: Populates `buf` with the update ratio R(z) = T(z+1)/T(z).
"""
function build_dcr!(buf, pre_func, base_func, ratio_func, z_min::Int, z_max::Int)

    # Topological Zero or Empty Range
    (z_min > z_max || buf.sign == 0) && return ZERO_DCR

    # prefactors 
    reset!(buf)
    pre_func(buf)
    m_root, m_rad = snapshot_square_root(buf)
    g_max_d = buf.max_d

    
    # Base Term at z_min
    reset!(buf, iseven(z_min) ? 1 : -1)
    base_func(buf, z_min)
    m_base = snapshot(buf)
    g_max_d = max(g_max_d, buf.max_d)

    # Ratio Sequence
    ratios = Vector{CyclotomicMonomial}(undef, z_max - z_min)
    for (i, z) in enumerate(z_min:z_max-1)
        reset!(buf, -1) # Ratios in these series usually carry a -1 flip
        ratio_func(buf, z)
        ratios[i] = snapshot(buf)
        g_max_d = max(g_max_d, buf.max_d)
    end

    return DCR(m_root, m_rad, m_base, ratios, z_min:z_max, g_max_d)
end


