
# -------------------------------------------------------------
#  -- Types and structures for DCR projections ---
# Discrete, analytical, exact and classical (q→1) projections
# -------------------------------------------------------------



#  ---- Buffers for discrete (level k) and analytical computations ------  

"""
    DiscreteBuffer{T}
A workspace for Log-Sum-Exp summation at roots of unity (Real-valued).
"""
mutable struct DiscreteBuffer{T}
    log_mags::Vector{T}
    signs::Vector{Int8}
    capacity::Int
end

DiscreteBuffer{T}(n) where T = DiscreteBuffer{T}(Vector{T}(undef, n), Vector{Int8}(undef, n), n)

"""
    AnalyticBuffer{T}
Workspace for Log-Sum-Exp summation in the generic complex plane (Complex-valued).
"""
mutable struct AnalyticBuffer{T}
    log_mags::Vector{T}
    phases::Vector{T}
    signs::Vector{Int8}
    capacity::Int
end

AnalyticBuffer{T}(n) where T = AnalyticBuffer{T}(Vector{T}(undef, n), Vector{T}(undef, n), Vector{Int8}(undef, n), n)

# global pre-allocated workspaces
const _WS_F64 = DiscreteBuffer{Float64}(8192)
const _WS_C64 = AnalyticBuffer{Float64}(8192)




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

#  --- Conversion to Float64 for sanity checks ---
# function Base.Complex{T}(res::CycloExactResult) where T <: AbstractFloat
#     # Project radical to numeric and multiply by the evaluated factor
#     rad_val = project_analytic(res.radical, exp(im*π/(res.k+2)))
#     return Complex{T}(rad_val * res.sum_factor)
# end



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
    print(io, "Exact SU(2)$k_sub Symbol: ")
    if iszero(res.sum_factor)
        print(io, "0")
    elseif is_identity(res.radical)
        print(io, res.sum_factor)
    else
        print(io, "√(", res.radical, ") × [Nemo Element]")
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