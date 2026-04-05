# ---------------------------------------------------------------------
# File: types.jl
# Core symbolic types for the Deferred Cyclotomic Representation (DCR).
# ---------------------------------------------------------------------

# --- Unicode Printing Utilities ---
const SUBSCRIPTS = Dict('0'=>'₀','1'=>'₁','2'=>'₂','3'=>'₃','4'=>'₄',
                        '5'=>'₅','6'=>'₆','7'=>'₇','8'=>'₈','9'=>'₉')

const SUPERSCRIPTS = Dict('0'=>'⁰','1'=>'¹','2'=>'²','3'=>'³','4'=>'⁴',
                          '5'=>'⁵','6'=>'⁶','7'=>'⁷','8'=>'⁸','9'=>'⁹','-'=>'⁻')

to_subscript(n::Integer)   = map(c -> SUBSCRIPTS[c], string(n))
to_superscript(n::Integer) = map(c -> SUPERSCRIPTS[c], string(n))

# ----------------------
# Cyclotomic Monomial
# ----------------------

"""
    CyclotomicMonomial
An unspecialized algebraic term: M(q) = σ * q^P * Π Φ_d(q)^e_d.
Stored as a sparse, sorted vector of (d => exponent) pairs.
"""
struct CyclotomicMonomial
    sign::Int                      # ∈ {-1, 0, 1}
    q_pow::Int                     # Power P of q
    phi_exps::Vector{Pair{Int,Int}} # Sparse (d => e_d)
    max_d::Int                     # Cached maximum index d
end

# Identity constructor (M = 1)
CyclotomicMonomial() = CyclotomicMonomial(1, 0, Pair{Int,Int}[], 0)

"""
    CycloBuffer
A dense mutable accumulator used during DCR compilation. 
Tracks 'max_d' to allow O(d) resets instead of O(N) scans.
"""
mutable struct CycloBuffer
    sign::Int
    q_pow::Int
    exps::Vector{Int}
    max_d::Int
end

CycloBuffer(cap::Int) = CycloBuffer(1, 0, zeros(Int, cap), 0)

"""
    DCR
The deferred representation: a combinatorial skeleton of a q-hypergeometric series.
Separates the combinatorial sum from geometrical preconditioning.
"""
struct DCR
    root::CyclotomicMonomial        # Rational prefactor (perfect squares)
    radical::CyclotomicMonomial     # Square-free radical (under the sqrt)
    base::CyclotomicMonomial        # Initial term at z_min
    ratios::Vector{CyclotomicMonomial} # Update sequence {R_z}
    z_range::UnitRange{Int}         # Summation bounds
    max_d::Int                      # Global max cyclotomic index
end

# --- Snapshot Logic (Dense -> Sparse) ---

function snapshot(buf::CycloBuffer)
    buf.sign == 0 && return CyclotomicMonomial(0, 0, Pair{Int,Int}[], 0)
    
    nnz = 0
    @inbounds for d in 1:buf.max_d
        if buf.exps[d] != 0; nnz += 1; end
    end
    
    sparse_exps = Vector{Pair{Int, Int}}(undef, nnz)
    idx = 1
    @inbounds for d in 1:buf.max_d
        e = buf.exps[d]
        if e != 0
            sparse_exps[idx] = d => e
            idx += 1
        end
    end
    return CyclotomicMonomial(buf.sign, buf.q_pow, sparse_exps, buf.max_d)
end

@inline function split_exp(e::Int)
    # Symmetric split: e = 2q + r where r ∈ {-1, 0, 1}
    # This prevents phase-flips in the complex plane compared to floor division
    r = (e & 1) == 0 ? 0 : (e > 0 ? 1 : -1)
    q = (e - r) >> 1
    return q, r
end

function snapshot_square_root(buf::CycloBuffer)
    nnz_sq, nnz_rad = 0, 0
    @inbounds for d in 1:buf.max_d
        e = buf.exps[d]
        if e != 0
            q, r = split_exp(e)
            q != 0 && (nnz_sq += 1)
            r != 0 && (nnz_rad += 1)
        end
    end

    sq_exps = Vector{Pair{Int, Int}}(undef, nnz_sq)
    rad_exps = Vector{Pair{Int, Int}}(undef, nnz_rad)
    
    idx_sq, idx_rad = 1, 1
    @inbounds for d in 1:buf.max_d
        e = buf.exps[d]
        if e != 0
            q, r = split_exp(e)
            if q != 0; sq_exps[idx_sq] = d => q; idx_sq += 1; end
            if r != 0; rad_exps[idx_rad] = d => r; idx_rad += 1; end
        end
    end
    
    q_q, q_r = split_exp(buf.q_pow)
    return (
        CyclotomicMonomial(buf.sign, q_q, sq_exps, buf.max_d), # Root gets the sign
        CyclotomicMonomial(1, q_r, rad_exps, buf.max_d)        # Radical is magnitude
    )
end

# --- REPL Display ---

function Base.show(io::IO, M::CyclotomicMonomial)
    M.sign == 0 && return print(io, "0")
    M.sign == -1 && print(io, "-")

    printed_lead = false
    if M.q_pow != 0
        print(io, M.q_pow == 1 ? "q" : "q" * to_superscript(M.q_pow))
        printed_lead = true
    end

    n = length(M.phi_exps)
    if n == 0
        !printed_lead && print(io, "1")
        return
    end

    # Balanced Truncation
    if n > 25
        head = 1:15
        tail = max(16, n-6):n
        indices = vcat(head, 0, tail)
    else
        indices = 1:n
    end

    for (i, idx) in enumerate(indices)
        (printed_lead || i > 1) && print(io, " ")
        if idx == 0
            print(io, "...")
            continue
        end
        d, e = M.phi_exps[idx]
        print(io, "Φ", to_subscript(d), e == 1 ? "" : to_superscript(e))
    end
end

@inline is_identity(m::CyclotomicMonomial) = m.sign == 1 && m.q_pow == 0 && isempty(m.phi_exps)

function Base.show(io::IO, dcr::DCR)
    println(io, "DCR (Deferred Cyclotomic Representation)")
    println(io, " ├─ Range    : ", dcr.z_range)
    println(io, " ├─ Max Index: d = ", dcr.max_d)

    print(io, " ├─ Radical  : ")
    if is_identity(dcr.radical)
        println(io, dcr.radical.sign == -1 ? "√(-1)" : "1")
    else
        println(io, "√(", dcr.radical, ")")
    end
    
    print(io, " ├─ Root     : ")
    println(io, is_identity(dcr.root) ? "1" : dcr.root)
    
    println(io, " ├─ Base Term: ", dcr.base)
    print(io,   " └─ Sequence : ", length(dcr.ratios), " update ratios {R_z}")
end




# --- Buffer Management ---

@inline function ensure_capacity!(buf::CycloBuffer, d_req::Int)
    curr = length(buf.exps)
    if d_req > curr
        new_cap = max(d_req + 10, curr * 2)
        resize!(buf.exps, new_cap)
        @inbounds fill!(view(buf.exps, curr+1:new_cap), 0)
    end
end

@inline function reset!(buf::CycloBuffer, sign::Int=1)
    buf.sign = sign
    buf.q_pow = 0
    @inbounds fill!(view(buf.exps, 1:buf.max_d), 0)
    buf.max_d = 0
end

@inline function update_exps!(buf::CycloBuffer, d::Int, p::Int)
    @inbounds buf.exps[d] += p
    d > buf.max_d && (buf.max_d = d)
end

# --- In-Place Arithmetic ---

function mul!(buf::CycloBuffer, a::CyclotomicMonomial, b::CyclotomicMonomial)
    (a.sign == 0 || b.sign == 0) && (buf.sign = 0; return buf)
    ensure_capacity!(buf, max(a.max_d, b.max_d))
    reset!(buf, a.sign * b.sign)
    buf.q_pow = a.q_pow + b.q_pow
    @inbounds for (d, e) in a.phi_exps; update_exps!(buf, d, e); end
    @inbounds for (d, e) in b.phi_exps; update_exps!(buf, d, e); end
    return buf
end

function div!(buf::CycloBuffer, a::CyclotomicMonomial, b::CyclotomicMonomial)
    b.sign == 0 && throw(DivideError())
    a.sign == 0 && (buf.sign = 0; return buf)
    ensure_capacity!(buf, max(a.max_d, b.max_d))
    reset!(buf, a.sign * b.sign)
    buf.q_pow = a.q_pow - b.q_pow
    @inbounds for (d, e) in a.phi_exps; update_exps!(buf, d, e); end
    @inbounds for (d, e) in b.phi_exps; update_exps!(buf, d, -e); end
    return buf
end

Base.:*(a::CyclotomicMonomial, b::CyclotomicMonomial) = snapshot(mul!(CycloBuffer(max(a.max_d, b.max_d)), a, b))
Base.:/(a::CyclotomicMonomial, b::CyclotomicMonomial) = snapshot(div!(CycloBuffer(max(a.max_d, b.max_d)), a, b))


# ---------------------------------------------------------------------
# Additional Utilities for types.jl
# ---------------------------------------------------------------------

# 1. Equality and Hashing (Critical for Caching)
function Base.:(==)(a::CyclotomicMonomial, b::CyclotomicMonomial)
    return a.sign == b.sign && a.q_pow == b.q_pow && a.phi_exps == b.phi_exps
end

function Base.hash(m::CyclotomicMonomial, h::UInt)
    h = hash(m.sign, h)
    h = hash(m.q_pow, h)
    return hash(m.phi_exps, h)
end

# 2. Inversion and Identity
Base.inv(m::CyclotomicMonomial) = CyclotomicMonomial(m.sign, -m.q_pow, [d => -e for (d, e) in m.phi_exps], m.max_d)

@inline Base.iszero(m::CyclotomicMonomial) = m.sign == 0
@inline Base.isone(m::CyclotomicMonomial)  = m.sign == 1 && m.q_pow == 0 && isempty(m.phi_exps)

# 3. Fused Multiply-Divide (The "Workhorse")
"""
    mul_div!(buf, a, b, c) -> buf
Calculates (a * b) / c in-place.
"""
function mul_div!(buf::CycloBuffer, a::CyclotomicMonomial, b::CyclotomicMonomial, c::CyclotomicMonomial)
    (a.sign == 0 || b.sign == 0) && (buf.sign = 0; return buf)
    c.sign == 0 && throw(DivideError())
    
    ensure_capacity!(buf, max(a.max_d, b.max_d, c.max_d))
    reset!(buf, a.sign * b.sign * c.sign)
    buf.q_pow = a.q_pow + b.q_pow - c.q_pow
    
    @inbounds for (d, e) in a.phi_exps; update_exps!(buf, d, e); end
    @inbounds for (d, e) in b.phi_exps; update_exps!(buf, d, e); end
    @inbounds for (d, e) in c.phi_exps; update_exps!(buf, d, -e); end
    return buf
end

# 4. Identity Monomial Constants
const ZERO_MONOMIAL = CyclotomicMonomial(0, 0, Pair{Int,Int}[], 0)
const ONE_MONOMIAL  = CyclotomicMonomial(1, 0, Pair{Int,Int}[], 0)





# types for projections 

# ---- exact projection ---- 

"""
    CycloExactResult{T}
The result of an exact projection. 
- `k`: The topological level.
- `radical`: The sparse square-free radical (preserved to avoid dense sqrts).
- `sum_factor`: The evaluated dense element in the Nemo field.
"""
struct CycloExactResult{T}
    k::Int
    radical::CyclotomicMonomial
    sum_factor::T
end

# --- REPL Display for Exact Results ---

function Base.show(io::IO, res::CycloExactResult)
    k_sub = to_subscript(res.k)
    println(io, "Exact SU(2)$k_sub Symbol")
    
    if iszero(res.sum_factor)
        println(io, " Value: 0")
        return
    end

    # Check if radical is 1
    if is_identity(res.radical)
        print(io, " Value: ")
        println(io, res.sum_factor)
    else
        println(io, " Value: √(A) * B")
        println(io, "  ├─ A (Radical): ", res.radical)
        print(io,   "  └─ B (Sum)    : ")
        # Truncate Nemo output if it's a massive polynomial
        str_val = string(res.sum_factor)
        if length(str_val) > 200
            println(io, str_val[1:100], " ... [Truncated] ... ", str_val[end-50:end])
        else
            println(io, str_val)
        end
    end
end




# ------ classical ----- 


struct ClassicalResult
    sign::Int
    sq_val::Rational{BigInt} # The squared symbol is always rational
end

function Base.show(io::IO, res::ClassicalResult)
    res.sign == 0 && return print(io, "0.0")
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    # If the square root is perfect, show the rational value
    if s_num^2 == num && s_den^2 == den
        print(io, res.sign < 0 ? "-" : "", Rational(s_num, s_den))
    else
        print(io, res.sign < 0 ? "-" : "", "√(", res.sq_val, ")")
    end
end



Base.Float64(res::ClassicalResult)  = res.sign * sqrt(Float64(res.sq_val))
Base.BigFloat(res::ClassicalResult) = res.sign * sqrt(BigFloat(res.sq_val))



# ------- buffers for projection to numbers ------ 


"""
    DiscreteBuffer{T}
Reusable buffer to eliminate heap allocations during the Log-Sum-Exp process.
"""
mutable struct DiscreteBuffer{T}
    log_mags::Vector{T}
    signs::Vector{Int8}
    DiscreteBuffer{T}(n) where T = new(Vector{T}(undef, n), Vector{Int8}(undef, n))
end

# Global workspace for the default Float64 path
const _WS_F64 = DiscreteBuffer{Float64}(4096)



mutable struct AnalyticBuffer{T}
    log_mags::Vector{T}
    phases::Vector{T}
    signs::Vector{Int8}
    AnalyticBuffer{T}(n) where T = new(Vector{T}(undef, n), Vector{T}(undef, n), Vector{Int8}(undef, n))
end

const _WS_C64 = AnalyticBuffer{Float64}(4096)