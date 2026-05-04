
# -------------------------------------------------------------
#  -- Types and structures for DCR projections ---
# Discrete, analytical, exact and classical (q→1) projections
# -------------------------------------------------------------


# --- Exact projection types (Nemo/cyclotomic field Q(ζ)) ---- 

"""
    CompositeExactResult{T}
The unified result of an exact algebraic projection into a cyclotomic field Q(ζ).
Represents a formal sum of terms, automatically grouping those that share the same square-free radical.
"""
struct CompositeExactResult{T}
    k::Int
    terms::Dict{CyclotomicMonomial, T}
end

#empty 
function CompositeExactResult(k::Int, ::Type{T}) where T
    return CompositeExactResult{T}(k, Dict{CyclotomicMonomial, T}())
end

#single term 
function CompositeExactResult(k::Int, radical::CyclotomicMonomial, factor::T) where T
    dict = Dict{CyclotomicMonomial, T}()
    if !iszero(factor)
        dict[radical] = factor
    end
    return CompositeExactResult{T}(k, dict)
end

# --- Base & Properties ---

Base.:(==)(a::CompositeExactResult, b::CompositeExactResult) = (a.k == b.k) && (a.terms == b.terms)
Base.iszero(comp::CompositeExactResult) = isempty(comp.terms)
Base.zero(comp::CompositeExactResult{T}) where T = CompositeExactResult(comp.k, T)

Base.length(comp::CompositeExactResult) = length(comp.terms)
Base.isempty(comp::CompositeExactResult) = isempty(comp.terms)
Base.iterate(comp::CompositeExactResult, state...) = iterate(comp.terms, state...)
Base.keys(comp::CompositeExactResult) = keys(comp.terms)
Base.values(comp::CompositeExactResult) = values(comp.terms)

# --- basic arithmetics ---

function Base.:+(a::CompositeExactResult{T}, b::CompositeExactResult{T}) where T
    a.k != b.k && error("Level mismatch: Cannot sum exact results from k=$(a.k) and k=$(b.k)")
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

Base.:-(a::CompositeExactResult{T}) where T = CompositeExactResult{T}(a.k, Dict(rad => -fac for (rad, fac) in a.terms))
Base.:-(a::CompositeExactResult, b::CompositeExactResult) = a + (-b)

function Base.:*(a::CompositeExactResult{T}, b::CompositeExactResult{T}) where T
    a.k != b.k && error("Level k mismatch in multiplication.")
    result = zero(a)
    buf = CycloBuffer(1024) # one buffer for all multiplications
    for (rad_a, fac_a) in a.terms
        for (rad_b, fac_b) in b.terms
            reset!(buf)
            mul!(buf, rad_a, rad_b)
            root_mono, new_rad = snapshot_square_root(buf)
            root_val = project_exact(root_mono, a.k) #project perfect square
            new_factor = fac_a * fac_b * root_val
            term = CompositeExactResult(a.k, new_rad, new_factor)
            result = result + term
        end
    end
    return result
end

function Base.:*(c, comp::CompositeExactResult{T}) where T
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
Base.:*(comp::CompositeExactResult, c) = c * comp

function Base.:+(c, comp::CompositeExactResult{T}) where T
    iszero(c) && return comp
    c_nemo = if c isa T
        c
    elseif !isempty(comp.terms)
        parent(first(values(comp.terms)))(c)
    else
        K, _ = cyclotomic_field(2 * (comp.k + 2), "ζ")
        K(c)
    end
    return CompositeExactResult(comp.k, ONE_MONOMIAL, c_nemo) + comp
end
Base.:+(comp::CompositeExactResult, c) = c + comp
Base.:-(c, comp::CompositeExactResult) = c + (-comp)
Base.:-(comp::CompositeExactResult, c) = comp + (-c)

function Base.:/(comp::CompositeExactResult{T}, c) where T
    iszero(c) && throw(DivideError())
    new_terms = Dict{CyclotomicMonomial, T}()
    for (rad, factor) in comp.terms
        new_terms[rad] = factor / c
    end
    return CompositeExactResult{T}(comp.k, new_terms)
end

function Base.:/(c, comp::CompositeExactResult)
    if length(comp.terms) == 1
        rad = first(keys(comp.terms))
        fac = first(values(comp.terms))
        if is_identity(rad)
            return c / fac
        end
    end
    error("Algebraic division by a CompositeExactResult with unresolved square roots is not implemented.")
end

function Base.:/(a::CompositeExactResult, b::CompositeExactResult)
    a.k != b.k && error("Level k mismatch in division")
    if length(b.terms) == 1
        rad = first(keys(b.terms))
        fac = first(values(b.terms))
        if is_identity(rad)
            return a / fac 
        end
    end
    error("Algebraic division by a CompositeExactResult with unresolved square roots is not implemented.")
end


# --- Classical projection types (Rational BigInt) ---- 

"""
    ClassicalResult
Represents σ * √(sq_val) where sq_val is a Rational{BigInt}.
Preserves precision for the classical Ponzano-Regge limit.
"""
struct ClassicalResult <: Real
    sign::Int
    sq_val::Rational{BigInt} 
end

# --- ClassicalResult: Base & Properties ---

Base.:(==)(a::ClassicalResult, b::ClassicalResult) = (a.sign == b.sign) && (a.sq_val == b.sq_val)
Base.isapprox(a::ClassicalResult, b::ClassicalResult; atol=1e-15) = abs(Float64(a) - Float64(b)) < atol

Base.zero(::Type{ClassicalResult}) = ClassicalResult(0, 0//1)
Base.zero(res::ClassicalResult)    = ClassicalResult(0, 0//1)
Base.iszero(res::ClassicalResult)  = res.sign == 0 || res.sq_val == 0//1

Base.one(::Type{ClassicalResult})  = ClassicalResult(1, 1//1)
Base.one(res::ClassicalResult)     = ClassicalResult(1, 1//1)
Base.isone(res::ClassicalResult)   = res.sign == 1 && res.sq_val == 1//1

Base.sign(res::ClassicalResult)    = res.sign
Base.abs(res::ClassicalResult)     = ClassicalResult(abs(res.sign), res.sq_val)


# --- ClassicalResult: math operations ---

Base.:-(a::ClassicalResult) = ClassicalResult(-a.sign, a.sq_val)
Base.inv(a::ClassicalResult) = iszero(a) ? throw(DivideError()) : ClassicalResult(a.sign, 1 // a.sq_val)

Base.:*(a::ClassicalResult, b::ClassicalResult) = ClassicalResult(a.sign * b.sign, a.sq_val * b.sq_val)
Base.:/(a::ClassicalResult, b::ClassicalResult) = a * inv(b)

function Base.:^(res::ClassicalResult, p::Int)
    p == 0 && return one(ClassicalResult)
    iszero(res) && return zero(ClassicalResult)
    # also handles negative powers (inverts the fraction)
    new_sign = iseven(p) ? 1 : res.sign
    return ClassicalResult(new_sign, res.sq_val^p)
end

function Base.:*(c::Union{Int, Rational}, res::ClassicalResult)
    iszero(c) && return zero(ClassicalResult)
    iszero(res) && return res
    return ClassicalResult(res.sign * sign(c), res.sq_val * c^2)
end
Base.:*(res::ClassicalResult, c::Union{Int, Rational}) = c * res

function Base.:/(res::ClassicalResult, c::Union{Int, Rational})
    iszero(c) && throw(DivideError())
    iszero(res) && return res
    return ClassicalResult(res.sign * sign(c), res.sq_val / c^2)
end

function _is_perfect_square(r::Rational)
    num, den = numerator(r), denominator(r)
    return isqrt(num)^2 == num && isqrt(den)^2 == den
end

function Base.:+(a::ClassicalResult, b::ClassicalResult)
    iszero(a) && return b
    iszero(b) && return a
    
    # merge exact perfect squares exactly
    if _is_perfect_square(a.sq_val) && _is_perfect_square(b.sq_val)
        val = Rational{BigInt}(a) + Rational{BigInt}(b)
        return ClassicalResult(sign(val), val^2)
    end
    
    # merge identical radicals exactly
    if a.sq_val == b.sq_val
        if a.sign == b.sign
            return ClassicalResult(a.sign, a.sq_val * 4) 
        else
            return zero(ClassicalResult)
        end
    end
    
    # fallback to BigFloat 
    return BigFloat(a) + BigFloat(b)
end

Base.:-(a::ClassicalResult, b::ClassicalResult) = a + (-b)

# --- ClassicalResult: mixed type fallbacks ---

function Base.:+(a::ClassicalResult, c::Union{Integer, Rational})
    iszero(c) && return a
    return a + ClassicalResult(sign(c), Rational{BigInt}(c^2))
end
Base.:+(c::Union{Integer, Rational}, a::ClassicalResult) = a + c
Base.:-(a::ClassicalResult, c::Union{Integer, Rational}) = a + (-c)
Base.:-(c::Union{Integer, Rational}, a::ClassicalResult) = c + (-a)

Base.:+(a::ClassicalResult, b::AbstractFloat) = BigFloat(a) + BigFloat(b)
Base.:+(a::AbstractFloat, b::ClassicalResult) = BigFloat(a) + BigFloat(b)
Base.:-(a::ClassicalResult, b::AbstractFloat) = BigFloat(a) - BigFloat(b)
Base.:-(a::AbstractFloat, b::ClassicalResult) = BigFloat(a) - BigFloat(b)

Base.:*(res::ClassicalResult, x::AbstractFloat) = BigFloat(res) * BigFloat(x)
Base.:*(x::AbstractFloat, res::ClassicalResult) = BigFloat(x) * BigFloat(res)
Base.:/(res::ClassicalResult, x::AbstractFloat) = BigFloat(res) / BigFloat(x)


# --- Conversions ---

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

# ---- REPL print ----

function _truncate_nemo_str(x, max_chunks=20)
    s = string(x)
    tokens = split(s, " ")
    #truncate at middle 
    if length(tokens) > 2 * max_chunks + 1
        first_part = join(tokens[1:max_chunks], " ")
        last_part  = join(tokens[end-max_chunks+1:end], " ")
        return first_part * " ... " * last_part
    end
    return s
end

function Base.show(io::IO, comp::CompositeExactResult)
    h = comp.k + 2
    h_sub = to_subscript(2h)
    print(io, "Exact Algebraic Result in ℚ(ζ$h_sub):\n  Value: ")
    
    if iszero(comp)
        print(io, "0")
        return
    end
    
    terms_strs = String[]
    for (rad, factor) in comp.terms
        fac_str = _truncate_nemo_str(factor)
        #don't print √(1)
        if is_identity(rad)
            push!(terms_strs, "($fac_str)")
        else
            #project radical to exact
            rad_val = project_exact(rad, comp.k)
            rad_str = _truncate_nemo_str(rad_val, 30)
            push!(terms_strs, "[ √($rad_str) × ($fac_str) ]")
        end
    end
    print(io, join(terms_strs, " + "))
end

function Base.show(io::IO, res::ClassicalResult)
    iszero(res) && return print(io, "0//1")
    num, den = numerator(res.sq_val), denominator(res.sq_val)
    s_num, s_den = isqrt(num), isqrt(den)
    
    if s_num^2 == num && s_den^2 == den
        print(io, res.sign < 0 ? "-" : "", Rational(s_num, s_den))
    else
        print(io, res.sign < 0 ? "-" : "", "√(", res.sq_val, ")")
    end
end