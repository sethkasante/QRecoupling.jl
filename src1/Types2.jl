# Types.jl

"""
    CycloMonomial

Represents exactly a product of cyclotomic polynomials evaluated at q:
Value = sign * z^(z_pow) * Π (Φ_d(z^2))^(exps[d])

Here, z = q^{1/2} = e^{i π / (k+2)}.
"""
struct CycloMonomial
    sign::Int
    z_pow::Int           # Handles the integer powers of z (which are fractional powers of q)
    exps::Dict{Int, Int} # Handles Π Φ_d(z^2)^{exps[d]}
end

# Exact Symbolic Multiplication
function Base.:*(a::CycloMonomial, b::CycloMonomial)
    small, large = length(a.exps) ≤ length(b.exps) ? (a,b) : (b,a)
    exps = copy(large.exps)
    for (d, e) in small.exps
        val = get(exps, d, 0) + e
        if val == 0
            delete!(exps, d)
        else
            exps[d] = val
        end
    end
    return CycloMonomial(a.sign * b.sign, a.z_pow + b.z_pow, exps)
end

# Exact Symbolic Division
function Base.:/(a::CycloMonomial, b::CycloMonomial)
    exps = copy(a.exps)
    for (d, e) in b.exps
        val = get(exps, d, 0) - e
        if val == 0
            delete!(exps, d)
        else
            exps[d] = val
        end
    end
    return CycloMonomial(a.sign * b.sign, a.z_pow - b.z_pow, exps)
end

Base.://(a::CycloMonomial, b::CycloMonomial) = a / b



# ----------------------------------------
# Human-Readable Output (Pretty Printing)
# ----------------------------------------

# Unicode mappings for clean subscripts and superscripts
const SUBSCRIPTS = Dict('0'=>'₀', '1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄', 
                        '5'=>'₅', '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉')
const SUPERSCRIPTS = Dict('0'=>'⁰', '1'=>'¹', '2'=>'²', '3'=>'³', '4'=>'⁴', 
                          '5'=>'⁵', '6'=>'⁶', '7'=>'⁷', '8'=>'⁸', '9'=>'⁹', '-'=>'⁻')

function to_subscript(n::Int)
    return map(c -> SUBSCRIPTS[c], string(n))
end

function to_superscript(n::Int)
    return map(c -> SUPERSCRIPTS[c], string(n))
end

# Overload Base.show to format the REPL output
function Base.show(io::IO, M::CycloMonomial)
    # Edge case for exact zero
    if M.sign == 0
        print(io, "0")
        return
    end

    parts = String[]
    
    # 1. z power (z = q^{1/2})
    if M.z_pow != 0
        if M.z_pow == 1
            push!(parts, "z")
        else
            push!(parts, "z" * to_superscript(M.z_pow))
        end
    end
    
    # 2. Cyclotomic polynomials Φ_d
    # Sort keys to ensure deterministic, mathematically neat ordering
    for d in sort(collect(keys(M.exps)))
        e = M.exps[d]
        e == 0 && continue
        
        base_str = "Φ" * to_subscript(d)
        if e == 1
            push!(parts, base_str)
        else
            push!(parts, base_str * to_superscript(e))
        end
    end
    
    # 3. Assembly
    if isempty(parts)
        # If there are no z powers or Φ terms, it's just the constant 1 or -1
        print(io, M.sign == -1 ? "-1" : "1")
    else
        # Prefix the sign, and join the terms with a space
        prefix = M.sign == -1 ? "-" : ""
        print(io, prefix * join(parts, " "))
    end
end



"""
    Generic6jResult

Holds the fully symbolic representation of a 6j symbol for ANY generic q.
"""
struct Generic6jResult
    pref_sq::CycloMonomial
    series::Vector{CycloMonomial}
end

function Base.show(io::IO, res::Generic6jResult)
    print(io, "√( ", res.pref_sq, " ) × ( ")
    join(io, res.series, "  +  ")
    print(io, " )")
end

"""
    Exact6jResult

Holds the exact algebraic evaluation of a 6j symbol at a specific root of unity.
"""
struct Exact6jResult
    k::Int
    pref_sq::nf_elem
    sum_cf::nf_elem
end

function Base.show(io::IO, res::Exact6jResult)
    print(io, "Exact SU(2)_$(res.k) Symbol:\n")
    print(io, "  Prefactor²: ", res.pref_sq, "\n")
    print(io, "  Racah Sum:  ", res.sum_cf)
end



# --- Exact and Numeric Models at level k ---

struct ExactSU2kModel
    k::Int
    K::AnticNumberField            
    z::nf_elem                     # Primitive 2N-th root of unity (z = q^{1/2})
    Phi_cache::Dict{Int, nf_elem}  # Cached evaluations of Φ_d(z^2)
end

struct NumericSU2kModel
    k::Int
    logqnfact::Vector{BigFloat}
end

#----   Objects to create ---- 
#TODO: Create an object for exact computations in cyclotomic fields Exact6j that has fields and returns a REPL of the form √{ζ₁₀³} × (ζ₁₀ - ζ₁₀³) ... 


#TODO: Create an object for Generic6j computations in cyclotomic monomials Generic? that has fields and returns a REPL with outputs of the form √{M} × [M₁, M₂, ...] OR of the form √{Φ₃⁻⁴ Φ₄⁻⁶} × [-z⁻⁶ Φ₂² Φ₃, z⁻¹⁰ Φ₂² Φ₃ Φ₄ , ...] 