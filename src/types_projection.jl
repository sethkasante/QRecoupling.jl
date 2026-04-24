
# -------------------------------------------------------------
#  -- Types and structures for DCR projections ---
# Discrete, analytical, exact and classical (q→1) projections
# -------------------------------------------------------------


# --- Exact projection types (Nemo/cyclotomic field Q(ζ)) ---- 

"""
    CycloExactResult{T}
The result of an exact algebraic projection into a cyclotomic field Q(ζ).
- `k`: The level of the TQFT.
- `radical`: The algebraic square-free part (a cyclotomic monomial).
- `sum_factor`: The evaluated dense element in the Nemo field.
"""
struct CycloExactResult{T}
    k::Int
    radical::CyclotomicMonomial
    sum_factor::T
end

# --- Exact Arithmetic & Verification ---

function Base.:(==)(a::CycloExactResult, b::CycloExactResult)
    a.k != b.k && return false
    return a.radical == b.radical && a.sum_factor == b.sum_factor
end

# Add two exact results if they share the same radical.
# function Base.:+(a::CycloExactResult, b::CycloExactResult)
#     a.k != b.k && error("Level mismatch")
#     if a.radical == b.radical
#         return CycloExactResult(a.k, a.radical, a.sum_factor + b.sum_factor)
#     elseif iszero(a.sum_factor)
#         return b
#     elseif iszero(b.sum_factor)
#         return a
#     else
#         # If radicals differ, verification is done by checking if Sum(res_i) is zero.
#         # This usually requires moving radicals into the field (if they are perfect squares).
#         error("Cannot sum exact results with different radicals algebraically.")
#     end
# end

# Base.:-(a::CycloExactResult, b::CycloExactResult) = a + CycloExactResult(b.k, b.radical, -b.sum_factor)
Base.iszero(res::CycloExactResult) = iszero(res.sum_factor)

# Scalar multiplication
Base.:*(c::Number, res::CycloExactResult) = CycloExactResult(res.k, res.radical, c * res.sum_factor)

Base.zero(res::CycloExactResult) = CycloExactResult(res.k, ZERO_MONOMIAL, zero(res.sum_factor))
Base.one(res::CycloExactResult)  = CycloExactResult(res.k, ONE_MONOMIAL, one(res.sum_factor))

function Base.:*(a::CycloExactResult{T}, b::CycloExactResult{T}) where T
    a.k != b.k && error("Level k mismatch in multiplication.")
    
    buf = CycloBuffer(max(a.radical.max_d, b.radical.max_d))
    mul!(buf, a.radical, b.radical)
    
    # Extract perfect squares
    root_mono, new_rad = snapshot_square_root(buf)
    
    root_val = project_exact(root_mono, a.k)
    new_factor = a.sum_factor * b.sum_factor * root_val
    
    return CycloExactResult{T}(a.k, new_rad, new_factor)
end


#  --- Conversion to Float64 for sanity checks ---
# function Base.Complex{T}(res::CycloExactResult) where T <: AbstractFloat
#     # Project radical to numeric and multiply by the evaluated factor
#     rad_val = project_analytic(res.radical, exp(im*π/(res.k+2)))
#     return Complex{T}(rad_val * res.sum_factor)
# end


"""
    CompositeExactResult{T}
Represents a formal sum of `CycloExactResult`s with potentially different radicals.
Automatically groups terms that share the same square-free radical.
"""
struct CompositeExactResult{T}
    k::Int
    terms::Dict{CyclotomicMonomial, T}
end

# Constructor to build an empty composite result for a given level k
function CompositeExactResult(k::Int, ::Type{T}) where T
    return CompositeExactResult{T}(k, Dict{CyclotomicMonomial, T}())
end

# Constructor to upgrade a single CycloExactResult into a Composite
function CompositeExactResult(res::CycloExactResult{T}) where T
    dict = Dict{CyclotomicMonomial, T}()
    if !iszero(res)
        dict[res.radical] = res.sum_factor
    end
    return CompositeExactResult{T}(res.k, dict)
end

# --- Addition: Cyclo + Cyclo ---
function Base.:+(a::CycloExactResult, b::CycloExactResult)
    a.k != b.k && error("Cannot sum exact results from different levels (k=$(a.k) vs k=$(b.k))")
    
    # If radicals match, keep it as a simple CycloExactResult
    if a.radical == b.radical
        return CycloExactResult(a.k, a.radical, a.sum_factor + b.sum_factor)
    end
    
    # If radicals differ, promote to CompositeExactResult
    comp = CompositeExactResult(a)
    return comp + b
end

# --- Addition: Composite + Cyclo ---
function Base.:+(comp::CompositeExactResult{T}, b::CycloExactResult{T}) where T
    comp.k != b.k && error("Level mismatch")
    iszero(b) && return comp
    
    new_terms = copy(comp.terms)
    if haskey(new_terms, b.radical)
        new_terms[b.radical] += b.sum_factor
        # Clean up exact zeros to keep the dictionary sparse
        iszero(new_terms[b.radical]) && delete!(new_terms, b.radical)
    else
        new_terms[b.radical] = b.sum_factor
    end
    
    return CompositeExactResult{T}(comp.k, new_terms)
end

# Allow commutative addition
Base.:+(a::CycloExactResult, b::CompositeExactResult) = b + a

# --- Addition: Composite + Composite ---
function Base.:+(a::CompositeExactResult{T}, b::CompositeExactResult{T}) where T
    a.k != b.k && error("Level mismatch")
    
    new_terms = copy(a.terms)
    for (rad, factor) in b.terms
        if haskey(new_terms, rad)
            new_terms[rad] += factor
            iszero(new_terms[rad]) && delete!(new_terms, rad)
        else
            new_terms[rad] = factor
        end
    end
    
    return CompositeExactResult{T}(a.k, new_terms)
end

# --- Subtraction ---
Base.:-(a::CycloExactResult) = CycloExactResult(a.k, a.radical, -a.sum_factor)
Base.:-(a::CompositeExactResult{T}) where T = CompositeExactResult{T}(a.k, Dict(rad => -fac for (rad, factor) in a.terms))

Base.:-(a::CycloExactResult, b::CycloExactResult) = a + (-b)
Base.:-(a::CompositeExactResult, b::CycloExactResult) = a + (-b)
Base.:-(a::CycloExactResult, b::CompositeExactResult) = a + (-b)
Base.:-(a::CompositeExactResult, b::CompositeExactResult) = a + (-b)






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

# --- Classical arithmetic (for identity verification) ---

function Base.:(==)(a::ClassicalResult, b::ClassicalResult)
    a.sign == b.sign && a.sq_val == b.sq_val
end

function Base.isapprox(a::ClassicalResult, b::ClassicalResult; atol=1e-15)
    return abs(Float64(a) - Float64(b)) < atol
end

# Multiplication is exact: (s1√v1) * (s2√v2) = (s1s2)√(v1v2)
function Base.:*(a::ClassicalResult, b::ClassicalResult)
    return ClassicalResult(a.sign * b.sign, a.sq_val * b.sq_val)
end

function Base.:/(a::ClassicalResult, b::ClassicalResult)
    return ClassicalResult(a.sign * b.sign, a.sq_val / b.sq_val)
end

#convert to Float
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

# addition is trickier: we usually convert to Float/BigFloat for sums of symbols
Base.:+(a::ClassicalResult, b::Number) = Float64(a) + b
Base.:+(a::Number, b::ClassicalResult) = a + Float64(b)


Base.:*(res::ClassicalResult, x::Number) = Float64(res) * x
Base.:*(x::Number, res::ClassicalResult) = Float64(res) * x


# ---- REPL ----

function Base.show(io::IO, res::CycloExactResult)
    k_sub = to_subscript(res.k)
    print(io, "Exact SU(2)$k_sub Symbol: \n")
    if iszero(res.sum_factor)
        print(io, "0")
    elseif is_identity(res.radical)
        print(io, res.sum_factor)
    else
        rad = project_exact(res.radical,res.k)
        print(io, "  Value: √(A) * B\n")
        print(io, "  ----------------\n")
        print(io, "  A (Radical) = ", rad, "\n")
        # print(io, "√(", res.radical, ") × ")
        print(io, "  B (Sum)     = ", res.sum_factor)
    end
end

function Base.show(io::IO, res::ClassicalResult)
    res.sign == 0 && return print(io, "0.0")
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    # Check for perfect square to simplify display
    if s_num^2 == num && s_den^2 == den
        print(io, res.sign < 0 ? "-" : "", Rational(s_num, s_den))
    else
        print(io, res.sign < 0 ? "-" : "", "√(", res.sq_val, ")")
    end
end