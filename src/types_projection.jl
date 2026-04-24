
# -------------------------------------------------------------
#  -- Types and structures for DCR projections ---
# Discrete, analytical, exact and classical (q→1) projections
# -------------------------------------------------------------


# --- Exact projection types (Nemo/cyclotomic field Q(ζ)) ---- 

"""
    CompositeExactResult{T}
The unified result of an exact algebraic projection into a cyclotomic field Q(ζ).
Represents a formal sum of terms, automatically grouping those that share the same square-free radical.
A standard single TQFT symbol is just a composite with one dictionary entry.
"""
struct CompositeExactResult{T}
    k::Int
    terms::Dict{CyclotomicMonomial, T}
end

# Empty (or zero) constructor
function CompositeExactResult(k::Int, ::Type{T}) where T
    return CompositeExactResult{T}(k, Dict{CyclotomicMonomial, T}())
end

# single term constructor
function CompositeExactResult(k::Int, radical::CyclotomicMonomial, factor::T) where T
    dict = Dict{CyclotomicMonomial, T}()
    if !iszero(factor)
        dict[radical] = factor
    end
    return CompositeExactResult{T}(k, dict)
end

# --- Exact Arithmetic & Verification ---

function Base.:(==)(a::CompositeExactResult, b::CompositeExactResult)
    a.k != b.k && return false
    return a.terms == b.terms
end

Base.iszero(comp::CompositeExactResult) = isempty(comp.terms)
Base.zero(comp::CompositeExactResult{T}) where T = CompositeExactResult(comp.k, T)

# Scalar multiplication
function Base.:*(c::Number, comp::CompositeExactResult{T}) where T
    iszero(c) && return zero(comp)
    new_terms = Dict{CyclotomicMonomial, T}()
    for (rad, factor) in comp.terms
        new_factor = c * factor
        if !iszero(new_factor)
            new_terms[rad] = new_factor
        end
    end
    return CompositeExactResult{T}(comp.k, new_terms)
end
Base.:*(comp::CompositeExactResult, c::Number) = c * comp

# --- Addition: Composite + Composite ---
function Base.:+(a::CompositeExactResult{T}, b::CompositeExactResult{T}) where T
    a.k != b.k && error("Level mismatch: Cannot sum exact results from k=$(a.k) and k=$(b.k)")
    
    new_terms = copy(a.terms)
    for (rad, factor) in b.terms
        if haskey(new_terms, rad)
            new_terms[rad] += factor
            # Clean up exact zeros to keep the dictionary sparse
            iszero(new_terms[rad]) && delete!(new_terms, rad)
        else
            new_terms[rad] = factor
        end
    end
    
    return CompositeExactResult{T}(a.k, new_terms)
end

# --- Subtraction ---
Base.:-(a::CompositeExactResult{T}) where T = CompositeExactResult{T}(a.k, Dict(rad => -fac for (rad, fac) in a.terms))
Base.:-(a::CompositeExactResult, b::CompositeExactResult) = a + (-b)

# --- Multiplication: Composite * Composite ---
function Base.:*(a::CompositeExactResult{T}, b::CompositeExactResult{T}) where T
    a.k != b.k && error("Level k mismatch in multiplication.")
    
    result = zero(a)
    
    # allocate one buffer to reuse for all radical multiplications
    buf = CycloBuffer(1024) 
    
    for (rad_a, fac_a) in a.terms
        for (rad_b, fac_b) in b.terms
            empty!(buf)
            mul!(buf, rad_a, rad_b)
            
            # extract perfect squares from the combined radicals
            root_mono, new_rad = snapshot_square_root(buf)
            
            # project the perfect square root back into Nemo 
            root_val = project_exact(root_mono, a.k)
            new_factor = fac_a * fac_b * root_val
            
            # add the new term to our running result (uses + to group radicals automatically)
            term = CompositeExactResult(a.k, new_rad, new_factor)
            result = result + term
        end
    end
    
    return result
end

# --- collection API for CompositeExactResult ---

Base.length(comp::CompositeExactResult) = length(comp.terms)
Base.isempty(comp::CompositeExactResult) = isempty(comp.terms)

# Iterate exactly like a dictionary, yielding (CyclotomicMonomial, factor) pairs
Base.iterate(comp::CompositeExactResult, state...) = iterate(comp.terms, state...)

Base.keys(comp::CompositeExactResult) = keys(comp.terms)
Base.values(comp::CompositeExactResult) = values(comp.terms)



# --- Classical projection types (Rational BigInt) ---- 

"""
    ClassicalResult
Represents σ * √(sq_val) where sq_val is a Rational{BigInt}.
Preserves infinite precision for the classical Ponzano-Regge limit.
"""
struct ClassicalResult
    sign::Int
    sq_val::Rational{BigInt} 
end

function Base.:(==)(a::ClassicalResult, b::ClassicalResult)
    a.sign == b.sign && a.sq_val == b.sq_val
end

function Base.isapprox(a::ClassicalResult, b::ClassicalResult; atol=1e-15)
    return abs(Float64(a) - Float64(b)) < atol
end

function Base.:*(a::ClassicalResult, b::ClassicalResult)
    return ClassicalResult(a.sign * b.sign, a.sq_val * b.sq_val)
end

function Base.:/(a::ClassicalResult, b::ClassicalResult)
    return ClassicalResult(a.sign * b.sign, a.sq_val / b.sq_val)
end

Base.Float64(res::ClassicalResult) = res.sign * Float64(sqrt(BigFloat(res.sq_val)))
Base.BigFloat(res::ClassicalResult) = res.sign * sqrt(BigFloat(res.sq_val))

function Base.Rational{T}(res::ClassicalResult) where {T <: Integer}
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    if s_num^2 == num && s_den^2 == den
        return res.sign * Rational{T}(T(s_num), T(s_den))
    else
        throw(InexactError(:Rational, Rational{T}, res))
    end
end

Base.:+(a::ClassicalResult, b::Number) = Float64(a) + b
Base.:+(a::Number, b::ClassicalResult) = a + Float64(b)
Base.:*(res::ClassicalResult, x::Number) = Float64(res) * x
Base.:*(x::Number, res::ClassicalResult) = Float64(res) * x

# ---- REPL Printers ----


# Helper to safely truncate large Nemo strings
function _truncate_nemo_str(x, max_chunks=15)
    s = string(x)
    tokens = split(s, " ")
    
    # If the polynomial has too many terms, truncate the middle
    if length(tokens) > 2 * max_chunks + 1
        first_part = join(tokens[1:max_chunks], " ")
        last_part  = join(tokens[end-max_chunks+1:end], " ")
        return first_part * " ... " * last_part
    end
    
    return s
end

function Base.show(io::IO, comp::CompositeExactResult)
    k_sub = to_subscript(comp.k)
    print(io, "Exact SU(2)$k_sub Composite:\n  Value: ")
    
    if iszero(comp)
        print(io, "0")
        return
    end
    
    terms_strs = String[]
    for (rad, factor) in comp.terms
        fac_str = _truncate_nemo_str(factor)
        
        # If the radical is just 1, we don't need to print √(1)
        if is_identity(rad)
            push!(terms_strs, "($fac_str)")
        else
            # Project the radical to see what 'A' is
            rad_val = project_exact(rad, comp.k)
            rad_str = _truncate_nemo_str(rad_val, 30)
            push!(terms_strs, "[ √($rad_str) × ($fac_str) ]")
        end
    end
    
    print(io, join(terms_strs, " + "))
end

function Base.show(io::IO, res::ClassicalResult)
    res.sign == 0 && return print(io, "0.0")
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    if s_num^2 == num && s_den^2 == den
        print(io, res.sign < 0 ? "-" : "", Rational(s_num, s_den))
    else
        print(io, res.sign < 0 ? "-" : "", "√(", res.sq_val, ")")
    end
end