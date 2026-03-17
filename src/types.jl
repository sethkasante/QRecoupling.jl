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
It isolates the perfect square root from the strictly square-free remainder.
"""
function snapshot_square_root(buf::SymbolicBuffer)
    nnz_sq = 0
    nnz_rem = 0
    
    @inbounds for e in buf.exps
        if e != 0
            q, r = divrem(e, 2)
            if q != 0; nnz_sq += 1; end
            if r != 0; nnz_rem += 1; end
        end
    end
    
    sq_exps = Vector{Pair{Int, Int}}(undef, nnz_sq)
    rem_exps = Vector{Pair{Int, Int}}(undef, nnz_rem)
    
    idx_sq = 1
    idx_rem = 1
    @inbounds for d in 1:length(buf.exps)
        e = buf.exps[d]
        if e != 0
            q, r = divrem(e, 2)
            if q != 0
                sq_exps[idx_sq] = d => q
                idx_sq += 1
            end
            if r != 0
                rem_exps[idx_rem] = d => r
                idx_rem += 1
            end
        end
    end
    
    q_z, r_z = divrem(buf.z_pow, 2)
    
    m_root = CycloMonomial(1, q_z, sq_exps)
    m_rem  = CycloMonomial(buf.sign, r_z, rem_exps)
    
    return m_root, m_rem # returns root and remainder
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
    pref_root::CycloMonomial       # Triangle coefficients (square-rooted part)
    pref_rem::CycloMonomial        # Remainder (square-free part inside the sqrt)
    m_min::CycloMonomial           # first term in Racah sum  
    ratios::Vector{CycloMonomial}  # ratios of hypergeometric steps
    z_range::UnitRange{Int}        # range of sum
    max_d::Int                     # maximum index d
end

function Base.show(io::IO, ::MIME"text/plain", res::CycloResult)
    n_ratios = length(res.ratios)
    
    println(io, "CycloResult (Hypergeometric Ratio Form)")
    println(io, "  ├─ Max Φ_d(q) req : d = ", res.max_d)
    println(io, "  ├─ Δ (Root Part)  : ", res.pref_root)
    println(io, "  ├─ Δ (Rem Part)   : √(", res.pref_rem, ")")
    println(io, "  ├─ M₀ (Base Term) : ", res.m_min)
    
    if n_ratios == 0
        print(io,   "  └─ Ratios (R)     : [Empty Series]")
    elseif n_ratios <= 3
        println(io, "  └─ Ratios (R)     : ", n_ratios, " terms")
        for i in 1:n_ratios
            prefix = (i == n_ratios) ? "       └─ R" : "       ├─ R"
            print(io, prefix, to_subscript(i), " : ", res.ratios[i])
            i < n_ratios && println(io) 
        end
    else
        println(io, "  └─ Ratios (R)     : ", n_ratios, " terms")
        println(io, "       ├─ R", to_subscript(1), " : ", res.ratios[1])
        println(io, "       ├─ R", to_subscript(2), " : ", res.ratios[2])
        println(io, "       ├─ ... (", n_ratios - 3, " more) ...")
        print(io,   "       └─ R", to_subscript(n_ratios), " : ", res.ratios[end])
    end
end



# ------ Exact Algebraic Structures (using Nemo.jl)  ------------

"""
    ExactResult{T}
A rigorously exact representation of a quantum symbol. 
Maintains the exact algebraic square-free remainder separated from the evaluated sum.
"""
struct ExactResult{T}
    k::Int
    pref_rem::CycloMonomial # The purely square-free algebraic remainder!
    sum_part::T
end

function Base.show(io::IO, res::ExactResult)
    k_sub = to_subscript(res.k)
    print(io, "Exact SU(2)$k_sub Symbol:\n")
    print(io, "  Prefactor: √(", res.pref_rem, ")\n")
    print(io, "  Sum(Σ):    ", res.sum_part)
end