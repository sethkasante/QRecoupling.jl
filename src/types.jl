# types.jl

# -----------------------------------------------------
# -- Core Symbolic Structures ---
# The Cyclotomic representation of 3j and 6j symbols
# -----------------------------------------------------


"""
    CycloMonomial

Represents a quantum prime factorization of basic objects (like q-integers and q-factorials) 
in terms of products of  `q` and irreducible cyclotomic polynomials Φ_d(q). 
Uses a sparse representation, storing ONLY non-zero exponents as a strictly ordered `Vector` 
of Pairs to guarantee O(1) cache-friendly loops during evaluation.
"""
struct CycloMonomial
    sign::Int
    z_pow::Int 
    exps::Vector{Pair{Int, Int}} #(d => exponent)
end


"""
    SymbolicBuffer

A dense, mutable accumulator used during the algorithmic construction phase.
"""
mutable struct SymbolicBuffer
    sign::Int
    z_pow::Int
    exps::Vector{Int} # dense for fast direct indexing
end

SymbolicBuffer(capacity::Int) = SymbolicBuffer(1, 0, zeros(Int, capacity))


"""
    snapshot(buf::SymbolicBuffer)

Freezes `SymbolicBuffer` into sparse `CycloMonomial`.
"""
function snapshot(buf::SymbolicBuffer)
    nnz = 0
    @inbounds for e in buf.exps
        if e != 0; nnz += 1; end
    end
    
    sparse_exps = Vector{Pair{Int, Int}}(undef, nnz)
    
    idx = 1
    @inbounds for d in 1:length(buf.exps)
        e = buf.exps[d]
        if e != 0
            sparse_exps[idx] = d => e
            idx += 1
        end
    end
    
    return CycloMonomial(buf.sign, buf.z_pow, sparse_exps)
end



"""
    snapshot_square_root(buf::SymbolicBuffer)

Intercepts the dense buffer and factorizes algebraic perfect squares in a single pass.
It isolates the perfect square root from the strictly square-free radical.
"""
function snapshot_square_root(buf::SymbolicBuffer)
    nnz_sq = 0
    nnz_rad = 0
    
    @inbounds for e in buf.exps
        if e != 0
            q, r = divrem(e, 2)
            if q != 0; nnz_sq += 1; end
            if r != 0; nnz_rad += 1; end
        end
    end
    
    sq_exps = Vector{Pair{Int, Int}}(undef, nnz_sq)
    rad_exps = Vector{Pair{Int, Int}}(undef, nnz_rad)
    
    idx_sq = 1
    idx_rad = 1
    @inbounds for d in 1:length(buf.exps)
        e = buf.exps[d]
        if e != 0
            q, r = divrem(e, 2)
            if q != 0
                sq_exps[idx_sq] = d => q
                idx_sq += 1
            end
            if r != 0
                rad_exps[idx_rad] = d => r
                idx_rad += 1
            end
        end
    end
    
    q_z, r_z = divrem(buf.z_pow, 2)
    
    m_root = CycloMonomial(1, q_z, sq_exps)
    m_rad  = CycloMonomial(buf.sign, r_z, rad_exps)
    
    return m_root, m_rad # returns root and radical
end


# ------- CycloMonomial Arithmetic Overloads ---- 

# ---- Base functions for products and divisions ---- 

function Base.:*(a::CycloMonomial, b::CycloMonomial)
    a.sign == 0 && return a
    b.sign == 0 && return b
    
    max_d = 0
    !isempty(a.exps) && (max_d = max(max_d, a.exps[end].first)) # Array is sorted!
    !isempty(b.exps) && (max_d = max(max_d, b.exps[end].first))
    
    buf = SymbolicBuffer(max_d)
    buf.sign = a.sign * b.sign
    buf.z_pow = a.z_pow + b.z_pow
    
    @inbounds for (d, e) in a.exps; buf.exps[d] += e; end
    @inbounds for (d, e) in b.exps; buf.exps[d] += e; end
    
    return snapshot(buf)
end

function Base.:/(a::CycloMonomial, b::CycloMonomial)
    a.sign == 0 && return a
    b.sign == 0 && throw(DivideError())
    
    max_d = 0
    !isempty(a.exps) && (max_d = max(max_d, a.exps[end].first))
    !isempty(b.exps) && (max_d = max(max_d, b.exps[end].first))
    
    buf = SymbolicBuffer(max_d)
    buf.sign = a.sign * b.sign
    buf.z_pow = a.z_pow - b.z_pow
    
    @inbounds for (d, e) in a.exps; buf.exps[d] += e; end
    @inbounds for (d, e) in b.exps; buf.exps[d] -= e; end
    
    return snapshot(buf)
end

Base.://(a::CycloMonomial, b::CycloMonomial) = a / b

function Base.inv(M::CycloMonomial)
    M.sign == 0 && throw(DivideError())
    inv_exps = [d => -e for (d, e) in M.exps] # keep it strictly sorted
    return CycloMonomial(M.sign, -M.z_pow, inv_exps)
end

# Check symbolic equality
function Base.:(==)(a::CycloMonomial, b::CycloMonomial)
    # Because the arrays are sequentially, they are strictly sorted.
    # We can just compare the Vectors directly!
    return a.sign == b.sign && a.z_pow == b.z_pow && a.exps == b.exps
end


 
# ----------------------------------------
# Human-Readable Output (Pretty Printing)
# ----------------------------------------

const SUBSCRIPTS = Dict('0'=>'₀', '1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄', 
                        '5'=>'₅', '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉')
const SUPERSCRIPTS = Dict('0'=>'⁰', '1'=>'¹', '2'=>'²', '3'=>'³', '4'=>'⁴', 
                          '5'=>'⁵', '6'=>'⁶', '7'=>'⁷', '8'=>'⁸', '9'=>'⁹', '-'=>'⁻')

to_subscript(n::Int) = map(c -> SUBSCRIPTS[c], string(n))
to_superscript(n::Int) = map(c -> SUPERSCRIPTS[c], string(n))

function Base.show(io::IO, M::CycloMonomial)
    M.sign == 0 && return print(io, "0")

    parts = String[]
    if M.z_pow != 0
        push!(parts, M.z_pow == 1 ? "z" : "z" * to_superscript(M.z_pow))
    end
    
    for (d, e) in M.exps
        base_str = "Φ" * to_subscript(d)
        push!(parts, e == 1 ? base_str : base_str * to_superscript(e))
    end
    
    if isempty(parts)
        print(io, M.sign == -1 ? "-1" : "1")
    else
        prefix = M.sign == -1 ? "-" : ""
        print(io, prefix * join(parts, " "))
    end
end


# ---------------------------------------
# -- Constructors for main results 
# (including 3j and 6j symbols) 
# ---------------------------------------


"""
    CycloResult

A high-performance, deferred-specialization representation of a quantum recoupling symbol.
Structures the evaluation as a hypergeometric ratio sequence to minimize algebraic divisions.
"""
struct CycloResult
    root::CycloMonomial       # Triangle coefficients (square-rooted part)
    radical::CycloMonomial        # Triangle coefficients Radical (the part inside the sqrt)
    base_term::CycloMonomial           # first term in Racah sum  
    ratios::Vector{CycloMonomial}  # ratios of hypergeometric steps
    z_range::UnitRange{Int}        # range of sum
    max_d::Int                     # maximum index d
end


function Base.show(io::IO, res::CycloResult)
    n_ratios = length(res.ratios)
    
    # Combine the rational prefactor with the base term for display purposes
    # The true mathematical form is: (Root * M0 * sqrt(Rad)) * [1 + R1 + R1*R2 + ...]
    overall_rat = res.root * res.base_term
    
    # Check if the radical is exactly 1 (sign=1, z_pow=0, no Φ polynomials)
    is_rad_one = res.radical.sign == 1 && res.radical.z_pow == 0 && isempty(res.radical.exps)
    
    println(io, "CycloResult (Hypergeometric Ratio Form)")
    println(io, "  ├─ Max Φ_d(q) required : d = ", res.max_d)
    
    if is_rad_one
        println(io, "  ├─ Overall Prefactor : ", overall_rat)
    else
        println(io, "  ├─ Prefactor (Rational) : ", overall_rat)
        println(io, "  ├─ Prefactor (Radical)  : √(", res.radical, ")")
    end
    
    if n_ratios == 0
        print(io,   "  └─ Ratios (R)     : [Empty Series]")
    elseif n_ratios <= 3
        println(io, "  └─ Ratios (R)     : ", n_ratios, " term(s)")
        for i in 1:n_ratios
            prefix = (i == n_ratios) ? "       └─ R" : "       ├─ R"
            print(io, prefix, to_subscript(i), " : ", res.ratios[i])
            i < n_ratios && println(io) 
        end
    else
        println(io, "  └─ Ratios (R)     : ", n_ratios, " term(s)")
        println(io, "       ├─ R", to_subscript(1), " : ", res.ratios[1])
        println(io, "       ├─ R", to_subscript(2), " : ", res.ratios[2])
        println(io, "       ├─ ... (", n_ratios - 3, " more) ...")
        print(io,   "       └─ R", to_subscript(n_ratios), " : ", res.ratios[end])
    end
end




# ------ Exact Algebraic Structures (using Nemo.jl)  ------------
# ==============================================================================
# Exact Algebraic Structures (Nemo.jl Engine)
# ==============================================================================

# """
#     ExactResult{T}

# A rigorously exact representation of a quantum symbol in a cyclotomic number field.
# Maintains the exact algebraic square-free remainder symbolically (`radical`) to 
# allow O(1) square-root extraction during multiplication, bypassing heavy CAS factorization.
# """
# struct ExactResult{T}
#     k::Int                  # Level k 
#     radical::CycloMonomial # Square free part inside sqrt: stored symbolically for fast multiplication!
#     factor::T             # The evaluated Nemo sum (Type T is a Nemo field element)
# end


"""
    HybridNemoResult{T}

The result of projecting a deferred CycloResult DAG into a exact cyclotomic field. 
To bypass dense polynomial square roots, the square-free radical geometry is 
preserved sparsely, while the rational remainder is evaluated into a dense field element.
"""
struct HybridNemoResult{T}
    k::Int                  # Topological level k
    radical::CycloMonomial  # The strictly square-free geometry (sparse)
    rational_factor::T      # The evaluated exact rational part (dense Nemo element)
end

# Helper function strictly for pretty-printing the square-free remainder
function _eval_nemo_print(m::CycloMonomial, k::Int, K, z)
    m.sign == 0 && return K(0)
    
    max_d = isempty(m.exps) ? 0 : m.exps[end].first
    V_exact, V_inv = get_phi_exact_table(max_d, k, z)
    
    A_val = K(1)
    @inbounds for (d, e) in m.exps
        if e > 0
            A_val *= (e == 1) ? V_exact[d] : V_exact[d]^e
        elseif e < 0
            A_val *= (e == -1) ? V_inv[d] : V_inv[d]^abs(e)
        end
    end
    
    A_val *= z^(m.z_pow)
    return m.sign == 1 ? A_val : -A_val
end

function Base.show(io::IO, res::HybridNemoResult)
    k_sub = to_subscript(res.k)
    print(io, "Exact SU(2)$k_sub Symbol:\n")
    
    if res.radical.sign == 0 || iszero(res.factor)
        print(io, "  0\n")
        return
    end

    # Dynamically extract the field and generator
    K = parent(res.factor)
    z = gen(K)
    
    A_val = _eval_nemo_print(res.radical, res.k, K, z)

    if isone(A_val)
        print(io, "  Value: ", res.factor)
    else
        print(io, "  Value: √(A) * B\n")
        print(io, "  ----------------\n")
        print(io, "  A (Radical) = ", A_val, "\n")
        print(io, "  B (Sum)     = ", res.factor)
    end
end

function Base.:(==)(a::HybridNemoResult, b::HybridNemoResult)
    a.k == b.k || return false
    iszero(a.factor) && return iszero(b.factor)
    iszero(b.factor) && return false
    return a.radical == b.radical && a.factor == b.factor
end

function Base.:+(a::HybridNemoResult, b::HybridNemoResult)
    @assert a.k == b.k "Cannot add results from different levels k"
    iszero(a.factor) && return b
    iszero(b.factor) && return a
    
    if a.radical == b.radical
        return HybridNemoResult(a.k, a.radical, a.factor + b.factor)
    else
        error("Cannot add ExactResults: They do not belong to the same topological square-class.")
    end
end

# ------------------------------------------------------
# Magic Multiplier (Extracts new squares dynamically!)
# -----------------------------------------------------
# ------------------------------------------------------------------------------
# Magic Multiplier (Extracts new squares dynamically!)
# ------------------------------------------------------------------------------
function Base.:*(a::HybridNemoResult, b::HybridNemoResult)
    @assert a.k == b.k "Cannot multiply different levels"
    
    iszero(a.factor) && return a
    iszero(b.factor) && return b
    
    m_prod = a.radical * b.radical
    
    sq_exps = Pair{Int,Int}[]
    rad_exps = Pair{Int,Int}[]
    for (d, e) in m_prod.exps
        q, r = divrem(e, 2)
        if q != 0; push!(sq_exps, d => q); end
        if r != 0; push!(rad_exps, d => r); end
    end
    q_z, r_z = divrem(m_prod.z_pow, 2)
    
    m_root = CycloMonomial(1, q_z, sq_exps)
    m_rad  = CycloMonomial(m_prod.sign, r_z, rad_exps)
    
    h = a.k + 2
    K = parent(a.factor)
    z = gen(K)
    
    max_d = isempty(m_root.exps) ? 0 : m_root.exps[end].first
    V_exact, V_inv = get_phi_exact_table(max_d, a.k, z)
    
    # Evaluate the newly formed root using the ultra-fast division-free projector
    exact_new_root = _project_ratio_nemo(m_root, V_exact, V_inv, z, h, K(0), K(1))
    
    new_sum = exact_new_root * a.factor * b.factor
    
    return HybridNemoResult(a.k, m_rad, new_sum)
end





