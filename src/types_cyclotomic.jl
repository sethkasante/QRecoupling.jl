
# -------------------------------------------------------------------------
#            --- Main symbolic types and structures for DRC ---
# Deferred cyclotomic representation (DRC) for q-hypergeometric series.
# -------------------------------------------------------------------------



# Printing utilities
const SUBSCRIPTS = Dict('0'=>'₀','1'=>'₁','2'=>'₂','3'=>'₃','4'=>'₄',
                        '5'=>'₅','6'=>'₆','7'=>'₇','8'=>'₈','9'=>'₉')

const SUPERSCRIPTS = Dict('0'=>'⁰','1'=>'¹','2'=>'²','3'=>'³','4'=>'⁴',
                          '5'=>'⁵','6'=>'⁶','7'=>'⁷','8'=>'⁸','9'=>'⁹','-'=>'⁻')

to_subscript(n::Integer)   = map(c -> SUBSCRIPTS[c], string(n))
to_superscript(n::Integer) = map(c -> SUPERSCRIPTS[c], string(n))




# ---- Struct Definitions ----- 

"""
    CyclotomicMonomial
An unspecialized algebraic term: M(q) = σ * qᴾ * Π Φ_d(q)ᵉᵈ.
Stored as a sparse, sorted vector of (d => exponent) pairs for efficiency.
Φ_d(q) are cyclotomic polynomials.
"""
struct CyclotomicMonomial
    sign::Int     # σ∈{-1,0,1}
    q_pow::Int    # Power P of q 
    phi_exps::Vector{Pair{Int,Int}}  # Sparse representation (d => exponent)
    max_d::Int    # maximum cyclotomic index d
end



# ---- Constants -----
const ZERO_MONOMIAL = CyclotomicMonomial(0, 0, Pair{Int,Int}[], 0)
const ONE_MONOMIAL  = CyclotomicMonomial(1, 0, Pair{Int,Int}[], 0)
# Identity cyclotomic monomial
CyclotomicMonomial() = ONE_MONOMIAL


"""
    CycloBuffer
A dense mutable accumulator used during DCR compilation to avoid 
intermediate sparse allocations. 
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
The Deferred Cyclotomic Representation: A combinatorial skeleton of a q-hypergeometric 
series. Separates combinatorial summation from numerical/geometrical projection.
"""
struct DCR
    root::CyclotomicMonomial     # perfect square-rooted prefactor
    radical::CyclotomicMonomial   # square-free radical (terms under the sqrt)
    base::CyclotomicMonomial      # initial term at z_min
    ratios::Vector{CyclotomicMonomial} # sequence of update ratios {R_z}
    z_range::UnitRange{Int}       # summation bounds
    max_d::Int                    # global maximum cyclotomic index
end

#zero DCR 
const ZERO_DCR = DCR(ZERO_MONOMIAL, ONE_MONOMIAL, ZERO_MONOMIAL, CyclotomicMonomial[], 0:-1, 0)


# ----- Buffer management & snapshots ------ 

@inline function ensure_capacity!(buf::CycloBuffer, d_req::Int)
    curr = length(buf.exps)
    #extend capacity if current size of exps < max d
    if d_req > curr
        new_cap = max(d_req + 10, curr * 2)
        resize!(buf.exps, new_cap)
        @inbounds fill!(view(buf.exps, curr+1:new_cap), 0)
    end
end

@inline function reset!(buf::CycloBuffer, sign::Int=1)
    #reset buf
    buf.sign = sign
    buf.q_pow = 0
    @inbounds fill!(view(buf.exps, 1:buf.max_d), 0)
    buf.max_d = 0
end

@inline function update_exps!(buf::CycloBuffer, d::Int, p::Int)
    #update exponents of Φ_d
    @inbounds buf.exps[d] += p
    d > buf.max_d && (buf.max_d = d)
end

#get sparse representation in CyclotomicMonomial form
function snapshot(buf::CycloBuffer)
    buf.sign == 0 && return ZERO_MONOMIAL
    
    #count non-zero exponents
    nnz = 0
    @inbounds for d in 1:buf.max_d
        buf.exps[d] != 0 && (nnz += 1)
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

#split exps as e = 2q + r, r∈{-1,0,1}
# @inline function split_exp(e::Int)
#     r = (e & 1) == 0 ? 0 : (e > 0 ? 1 : -1)
#     q = (e - r) >> 1
#     return q, r
# end

#Check: any branch choice related issue?
#split exps as e = 2q + r, r∈{0,1}
@inline function split_exp(e::Int)
    r = e & 1           
    q = (e - r) >> 1   
    return q, r
end
#snapshot for square root of cyclotomic monomials
#split into perfect square and radical
function snapshot_square_root(buf::CycloBuffer)
    buf.sign == 0 && return ZERO_MONOMIAL

    #count number of perfect squares and radicals 
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
    
    q_q, q_r = split_exp(buf.q_pow)
    return ( CyclotomicMonomial(buf.sign, q_q, sq_exps, buf.max_d), 
                CyclotomicMonomial(1, q_r, rad_exps, buf.max_d) )
end




#  --- Arithmetic operations on Cyclotomic monomials ----- 

function mul!(buf::CycloBuffer, a::CyclotomicMonomial, b::CyclotomicMonomial)
    (a.sign == 0 || b.sign == 0) && (buf.sign = 0; return buf)

    ensure_capacity!(buf, max(a.max_d, b.max_d))
    reset!(buf, a.sign * b.sign)
    buf.q_pow = a.q_pow + b.q_pow

    @inbounds for (d, e) in a.phi_exps 
        update_exps!(buf, d, e) 
    end

    @inbounds for (d, e) in b.phi_exps 
        update_exps!(buf, d, e) 
    end
    return buf
end

function div!(buf::CycloBuffer, a::CyclotomicMonomial, b::CyclotomicMonomial)
    b.sign == 0 && throw(DivideError())
    a.sign == 0 && (buf.sign = 0; return buf)

    ensure_capacity!(buf, max(a.max_d, b.max_d))
    reset!(buf, a.sign * b.sign)
    buf.q_pow = a.q_pow - b.q_pow

    @inbounds for (d, e) in a.phi_exps 
        update_exps!(buf, d, e)
    end

    @inbounds for (d, e) in b.phi_exps 
        update_exps!(buf, d, -e)
    end
    return buf
end

# (M_a * M_b )/ M_c
function mul_div!(buf::CycloBuffer, a::CyclotomicMonomial, b::CyclotomicMonomial, c::CyclotomicMonomial)
    (a.sign == 0 || b.sign == 0) && (buf.sign = 0; return buf)
    c.sign == 0 && throw(DivideError())

    ensure_capacity!(buf, max(a.max_d, b.max_d, c.max_d))
    reset!(buf, a.sign * b.sign * c.sign)
    buf.q_pow = a.q_pow + b.q_pow - c.q_pow

    @inbounds for (d, e) in a.phi_exps
        update_exps!(buf, d, e)
    end
    @inbounds for (d, e) in b.phi_exps
        update_exps!(buf, d, e)
    end
    @inbounds for (d, e) in c.phi_exps
        update_exps!(buf, d, -e)
    end
    return buf
end

Base.:*(a::CyclotomicMonomial, b::CyclotomicMonomial) = snapshot(mul!(CycloBuffer(max(a.max_d, b.max_d)), a, b))
Base.:/(a::CyclotomicMonomial, b::CyclotomicMonomial) = snapshot(div!(CycloBuffer(max(a.max_d, b.max_d)), a, b))
Base.inv(m::CyclotomicMonomial) = CyclotomicMonomial(m.sign, -m.q_pow, [d => -e for (d, e) in m.phi_exps], m.max_d)

Base.:-(m::CyclotomicMonomial) = CyclotomicMonomial(-m.sign, m.q_pow, m.phi_exps, m.max_d)

# Integer Powers
function Base.:^(m::CyclotomicMonomial, p::Int)
    m.sign == 0 && return p == 0 ? ONE_MONOMIAL : ZERO_MONOMIAL
    p == 0 && return ONE_MONOMIAL
    p == 1 && return m
    p == -1 && return inv(m)
    
    #sign: only flips if base is negative and power is odd
    new_sign = (m.sign == -1 && isodd(p)) ? -1 : 1
    new_q_pow = m.q_pow * p
    new_exps = Pair{Int,Int}[d => (e * p) for (d, e) in m.phi_exps]
    
    return CyclotomicMonomial(new_sign, new_q_pow, new_exps, m.max_d)
end


#  --- scalar Multiplication for Cyclotomic Monomials ----- 

function Base.:*(c::Int, m::CyclotomicMonomial)
    c == 0 && return ZERO_MONOMIAL
    c == 1 && return m
    c == -1 && return -m
    throw(ArgumentError("CyclotomicMonomials only support scalar coefficients of -1, 0, or 1. Got: $c"))
end

Base.:*(m::CyclotomicMonomial, c::Integer) = c * m

function Base.:/(m::CyclotomicMonomial, c::Integer)
    c == 0 && throw(DivideError())
    c == 1 && return m
    c == -1 && return -m
    throw(ArgumentError("CyclotomicMonomials only support scalar division by -1 or 1. Got: $c"))
end


Base.:/(c::Int, m::CyclotomicMonomial) = c * inv(m)



# ----- Comparisons ----- 
# utilities for comparing cyclotomic monomials 
@inline Base.iszero(m::CyclotomicMonomial) = m.sign == 0
@inline Base.isone(m::CyclotomicMonomial)  = m.sign == 1 && m.q_pow == 0 && isempty(m.phi_exps)

@inline function is_identity(m::CyclotomicMonomial)
    return m.sign == 1 && m.q_pow == 0 && isempty(m.phi_exps)
end

function Base.:(==)(a::CyclotomicMonomial, b::CyclotomicMonomial)
    return a.sign == b.sign && a.q_pow == b.q_pow && a.phi_exps == b.phi_exps
end

function Base.hash(m::CyclotomicMonomial, h::UInt)
    h = hash(m.sign, h)
    h = hash(m.q_pow, h)
    return hash(m.phi_exps, h)
end

"""
    _phi_exponent(m::CyclotomicMonomial, h::Int)
Returns the exponent of Φ_h in the monomial. Returns 0 if not present.
"""
@inline function _phi_exponent(m::CyclotomicMonomial, h::Int)
    for (d, e) in m.phi_exps
        d == h && return e
    end
    return 0
end

# ---- REPL ---- 

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

    # Balanced Truncation for readability
    if n > 25
        indices = vcat(1:15, 0, (n-6):n)
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