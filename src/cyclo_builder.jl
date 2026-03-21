
# -------------------------------------------------------------------------
# `CycloMonomial` construction in terms of cyclotomic polynomials (Φ_d) 
#  symbolic construction of the Racah-Wigner 3j and 6j symbols 
# -------------------------------------------------------------------------

const ZERO_MONOMIAL = CycloMonomial(0, 0, Pair{Int,Int}[])


# Ensure allocation-Free Operations -----

# resize with zeros if necessary to ensure capacity
@inline function _ensure_capacity!(buf::SymbolicBuffer, n::Int)
    old_len = length(buf.exps)
    if n > old_len
        resize!(buf.exps, n)
        fill!(view(buf.exps, (old_len+1):n), 0)
    end
    return nothing
end

# multiply term ny qintegers -> add to exponent vector
function add_qint!(buf::SymbolicBuffer, n::Int, power::Int=1)
    n <= 1 && return nothing
    # z^(1-n) : contribution from [n]
    buf.z_pow += power * (1 - n)
    _ensure_capacity!(buf, n)
    @inbounds for d in 2:n
        if n % d == 0
            buf.exps[d] += power
        end
    end
    return nothing
end

function add_qfact!(buf::SymbolicBuffer, n::Int, power::Int=1)
    n <= 1 && return nothing
    # z^(-n(n-1)/2) : contribution from [n]!
    buf.z_pow += power * ((n*(1 - n)) ÷ 2)
    _ensure_capacity!(buf, n)
    @inbounds for d in 2:n
        buf.exps[d] += power * div(n, d)
    end
    return nothing
end

# helper to avoid closure allocations inside the main builders
@inline function get_max_d(buf::SymbolicBuffer)
    idx = findlast(!iszero, buf.exps)
    return isnothing(idx) ? 0 : idx
end


#---- Structural Combinatorics ------- 

# square of triangle coefficients 
function qdelta2_symb!(buf::SymbolicBuffer, j1::Spin, j2::Spin, j3::Spin)
    add_qfact!(buf, Int(j1+j2-j3), 1)
    add_qfact!(buf, Int(j1-j2+j3), 1)
    add_qfact!(buf, Int(-j1+j2+j3), 1)
    add_qfact!(buf, Int(j1+j2+j3+1), -1)
    return nothing
end

# products of 4 triangle coefficients forming a tetrahedron
function qtricoeff2_symb!(buf::SymbolicBuffer, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    # triangle 1 (j1,j2,j3)
    add_qfact!(buf, Int(j1+j2-j3), 1) 
    add_qfact!(buf, Int(j1-j2+j3), 1)
    add_qfact!(buf, Int(-j1+j2+j3), 1)
    add_qfact!(buf, Int(j1+j2+j3+1), -1)
    # triangle 2 (j1, j5, j6)
    add_qfact!(buf, Int(j1+j5-j6), 1)
    add_qfact!(buf, Int(j1-j5+j6), 1)
    add_qfact!(buf, Int(-j1+j5+j6), 1)
    add_qfact!(buf, Int(j1+j5+j6+1), -1)
    # triangle 3 (j2, j4, j6)
    add_qfact!(buf, Int(j2+j4-j6), 1)
    add_qfact!(buf, Int(j2-j4+j6), 1)
    add_qfact!(buf, Int(-j2+j4+j6), 1)
    add_qfact!(buf, Int(j2+j4+j6+1), -1)
    # triangle 4 (j3, j4, j5)
    add_qfact!(buf, Int(j3+j4-j5), 1)
    add_qfact!(buf, Int(j3-j4+j5), 1)
    add_qfact!(buf, Int(-j3+j4+j5), 1)
    add_qfact!(buf, Int(j3+j4+j5+1), -1)
    return nothing
end




# --- Main construction of quantum 6j symbol (algebraic) --- 


"""
    q6j_cyclo(j1, j2, j3, j4, j5, j6)

Constructs the exact algebraic representation of the quantum 6j-symbol.
Formulates the Racah-Wigner sum as a highly optimized hypergeometric ratio sequence `(1 + R_1 + R_1 R_2 + ...)`
to guarantee O(1) mathematical divisions during downstream evaluation.
"""
function q6j_cyclo(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    
    α = (Int(j1+j2+j3), Int(j1+j5+j6), Int(j2+j4+j6), Int(j3+j4+j5))
    β = (Int(j1+j2+j4+j5), Int(j1+j3+j4+j6), Int(j2+j3+j5+j6))
    
    z_min, z_max = max(α...), min(β...)
    cap = max(z_max + 2, β...) 
    buf = SymbolicBuffer(cap) # initialize
    max_d_tracker = 0
    
    # Prefactor extraction -> split into root and remainder
    qtricoeff2_symb!(buf, j1, j2, j3, j4, j5, j6)
    max_d_tracker = max(max_d_tracker, get_max_d(buf))
    root, radical = snapshot_square_root(buf)
    
    # topological zeros (triangle inequality violation)
    if z_min > z_max
        return CycloResult(root, radical, ZERO_MONOMIAL, CycloMonomial[], 0:-1, max_d_tracker)
    end
    
    # ---- construct the hypergeometric representation of the sum ----
    
    # -- get base term M(z_min) in sum --
    # reset buffer 
    buf.sign = iseven(z_min) ? 1 : -1
    buf.z_pow = 0
    fill!(buf.exps, 0)
    
    add_qfact!(buf, z_min+1, 1) #[z+1]!
    for a in α
        add_qfact!(buf, z_min-a, -1)
    end

    for b in β
        add_qfact!(buf, b-z_min, -1)
    end
    
    max_d_tracker = max(max_d_tracker, get_max_d(buf))
    m_min = snapshot(buf)
    
    # hypergeometric ratios: R_z
    ratios = CycloMonomial[]
    sizehint!(ratios, max(0, z_max - z_min))
    
    for z in z_min : z_max-1

        buf.sign = -1
        buf.z_pow = 0
        fill!(buf.exps, 0)
        add_qint!(buf, z + 2, 1) 
        for b in β
            add_qint!(buf, b - z, 1)
        end
        for a in α
            add_qint!(buf, z + 1 - a, -1)
        end
        
        max_d_tracker = max(max_d_tracker, get_max_d(buf))
        push!(ratios, snapshot(buf))
    end
    
    return CycloResult(root, radical, m_min, ratios, z_min:z_max, max_d_tracker)
end


# ----- Quantum 3j Symbol (CycloResult) ----- 

"""
    q3j_cyclo(j1, j2, j3, m1, m2, m3)

Constructs the exact algebraic representation of the quantum 3j-symbol.
"""
function q3j_cyclo(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    cap_pref = Int(j1 + j2 + j3 + abs(m1) + abs(m2) + 2)
    buf = SymbolicBuffer(cap_pref) #start 
    max_d_tracker = 0

    # Prefactor^2 
    qdelta2_symb!(buf, j1, j2, j3) 
    add_qfact!(buf, Int(j1+m1), 1)
    add_qfact!(buf, Int(j1-m1), 1)
    add_qfact!(buf, Int(j2+m2), 1)
    add_qfact!(buf, Int(j2-m2), 1)
    add_qfact!(buf, Int(j3-m1-m2), 1)
    add_qfact!(buf, Int(j3+m1+m2), 1)
    
    max_d_tracker = max(max_d_tracker, get_max_d(buf))
    # split to root and remainder 
    root, radical = snapshot_square_root(buf)

    # Bounds for sum 
    α1 = Int(j3 - j2 + m1)
    α2 = Int(j3 - j1 - m2)

    β1 = Int(j1 + j2 - j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    z_min = max(-α1, -α2, 0)
    z_max = min(β1, β2, β3)

    if z_min > z_max
        empty_m = CycloMonomial(0, 0, Pair{Int,Int}[])
        return CycloResult(root, radical, empty_m, CycloMonomial[], 0:-1, max_d_tracker)
    end

    # --- initial summand M(z_min) ---
    cap_m = max(z_max + max(α1, α2), β1, β2, β3) + 1
    _ensure_capacity!(buf, cap_m) 
    
    buf.sign = isodd(Int(z_min + α1 - α2)) ? -1 : 1
    buf.z_pow = 0
    fill!(buf.exps, 0)
    
    add_qfact!(buf, z_min, -1)
    add_qfact!(buf, α1 + z_min, -1)
    add_qfact!(buf, α2 + z_min, -1)
    add_qfact!(buf, β1 - z_min, -1)
    add_qfact!(buf, β2 - z_min, -1)
    add_qfact!(buf, β3 - z_min, -1)
    
    max_d_tracker = max(max_d_tracker, get_max_d(buf))
    m_min = snapshot(buf)

    # --- 4. Recursive Ratios R_z ---
    ratios = CycloMonomial[]
    sizehint!(ratios, z_max - z_min)
    
    for z in z_min+1 : z_max
        buf.sign = -1; buf.z_pow = 0; fill!(buf.exps, 0)
        
        add_qint!(buf, β1 - z + 1, 1)
        add_qint!(buf, β2 - z + 1, 1)
        add_qint!(buf, β3 - z + 1, 1)
        add_qint!(buf, z, -1)
        add_qint!(buf, α1 + z, -1)
        add_qint!(buf, α2 + z, -1)
        
        max_d_tracker = max(max_d_tracker, get_max_d(buf))
        push!(ratios, snapshot(buf))
    end

    return CycloResult(root, radical, m_min, ratios, z_min:z_max, max_d_tracker)
end
