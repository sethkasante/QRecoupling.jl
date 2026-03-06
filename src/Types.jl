# src/Types.jl

"""
    CycloMonomial

Represents a product of cyclotomic polynomials evaluated at q:
Value = sign * z^(z_pow) * Π (Φ_d(q))^(exps[d])  where z = q^{1/2} = e^{i π / (k+2)}.

Uses a dense Vector{Int} for `exps` where the index `d` represents Φ_d.
"""
struct CycloMonomial
    sign::Int
    z_pow::Int # integer powers of z    
    exps::Vector{Int} # store exponents
end

# ----------------------
# Algebraic Operations
# ----------------------
#Multiplication
function Base.:*(a::CycloMonomial, b::CycloMonomial)
    len_a = length(a.exps)
    len_b = length(b.exps)
    max_len = max(len_a, len_b)
    
    new_exps = zeros(Int, max_len)
    
    # Bypass bounds checking and allow SIMD vectorization for raw speed
    @inbounds @simd for i in 1:len_a
        new_exps[i] += a.exps[i]
    end
    @inbounds @simd for i in 1:len_b
        new_exps[i] += b.exps[i]
    end
    
    # Strip trailing zeros to save memory and shorten future loops
    while !isempty(new_exps) && new_exps[end] == 0
        pop!(new_exps)
    end
    
    return CycloMonomial(a.sign * b.sign, a.z_pow + b.z_pow, new_exps)
end

#Division
function Base.:/(a::CycloMonomial, b::CycloMonomial)
    len_a = length(a.exps)
    len_b = length(b.exps)
    max_len = max(len_a, len_b)
    
    new_exps = zeros(Int, max_len)
    
    @inbounds @simd for i in 1:len_a
        new_exps[i] += a.exps[i]
    end
    @inbounds @simd for i in 1:len_b
        new_exps[i] -= b.exps[i]
    end
    
    while !isempty(new_exps) && new_exps[end] == 0
        pop!(new_exps)
    end
    
    return CycloMonomial(a.sign * b.sign, a.z_pow - b.z_pow, new_exps)
end

Base.://(a::CycloMonomial, b::CycloMonomial) = a / b

# Inversion to avoid division loops when numerator is 1
function Base.inv(M::CycloMonomial)
    M.sign == 0 && throw(DivideError())
    return CycloMonomial(M.sign, -M.z_pow, -M.exps)
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
    # Edge case for exact zero
    if M.sign == 0
        print(io, "0")
        return
    end

    parts = String[]
    
    # 1. Power of z = q^{1/2}
    if M.z_pow != 0
        push!(parts, M.z_pow == 1 ? "z" : "z" * to_superscript(M.z_pow))
    end
    
    # 2. Cyclotomic polynomials Φ_d (Index is the cyclotomic base)
    for d in 1:length(M.exps)
        e = M.exps[d]
        e == 0 && continue
        
        base_str = "Φ" * to_subscript(d)
        push!(parts, e == 1 ? base_str : base_str * to_superscript(e))
    end
    
    # 3. Final Assembly
    if isempty(parts)
        print(io, M.sign == -1 ? "-1" : "1")
    else
        prefix = M.sign == -1 ? "-" : ""
        print(io, prefix * join(parts, " "))
    end
end

# ----------------------------------------
# Result Structs
# ----------------------------------------

"""
    GenericResult

Holds the fully symbolic representation of a quantum symbol for ANY generic q.
"""
struct GenericResult
    pref_sq::CycloMonomial
    series::Vector{CycloMonomial}
end

function Base.show(io::IO, res::GenericResult)
    print(io, "√(", res.pref_sq, ") × (")
    
    n_terms = length(res.series)
    if n_terms <= 5
        join(io, res.series, "  +  ")
    else
        # Print first 2 terms, an ellipsis, and the last 2 terms
        print(io, res.series[1], "  +  ", res.series[2], "  +  ... (", n_terms - 4, " more terms) ...  +  ", res.series[end-1], "  +  ", res.series[end])
    end
    print(io, ")")
end

# function Base.show(io::IO, res::GenericResult)
#     print(io, "√(", res.pref_sq, ") × (")
#     join(io, res.series, "  +  ")
#     print(io, ")")
# end

"""
    ExactResult

Holds the exact algebraic evaluation of a quantum symbol at a specific root of unity.
"""
struct ExactResult
    k::Int
    pref_sq::nf_elem
    sum_cf::nf_elem
end

function Base.show(io::IO, res::ExactResult)
    k_sub = to_subscript(res.k)
    print(io, "Exact SU(2)$k_sub Symbol:\n")
    print(io, "  Prefactor(Δ²): ", res.pref_sq, "\n")
    print(io, "  Racah Sum(Σ):  ", res.sum_cf)
end

# ----------------------------------------
# Model Definitions
# ----------------------------------------

# struct ExactSU2kModel
#     k::Int
#     K::AnticNumberField            
#     z::nf_elem                     # Primitive 2N-th root of unity (z = q^{1/2})
#     Phi_eval::Vector{nf_elem}  # Cached evaluations of Φ_d(z^2)
# end

struct NumericSU2kModel
    k::Int
    logqnfact::Vector{BigFloat}
end