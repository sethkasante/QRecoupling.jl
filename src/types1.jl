# ==============================================================================
# File: types_projections.jl
# Result types and execution buffers for DCR projections.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Exact Projection Types (Nemo/Cyclotomic Field)
# ------------------------------------------------------------------------------

"""
    CycloExactResult{T}
The result of an exact algebraic projection into a cyclotomic field.
- `k`: The level of the TQFT.
- `radical`: The algebraic square-free part (monomial).
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

# Conversion to Float64 for sanity checks
function Base.Complex{T}(res::CycloExactResult) where T <: AbstractFloat
    # Project radical to numeric and multiply by the evaluated factor
    rad_val = project_monomial_analytic(res.radical, exp(im*π/(res.k+2)))
    return Complex{T}(rad_val * res.sum_factor)
end

# ------------------------------------------------------------------------------
# 2. Classical Projection Types (Rational BigInt)
# ------------------------------------------------------------------------------

"""
    ClassicalResult
Represents σ * √(sq_val) where sq_val is a Rational{BigInt}.
Preserves infinite precision for the classical Ponzano-Regge limit.
"""
struct ClassicalResult
    sign::Int
    sq_val::Rational{BigInt} 
end

# --- Classical Arithmetic (For Identity Verification) ---

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

# Addition is trickier; we usually convert to Float/BigFloat for sums of symbols
Base.:+(a::ClassicalResult, b::Number) = Float64(a) + b
Base.:+(a::Number, b::ClassicalResult) = a + Float64(b)

# ------------------------------------------------------------------------------
# 3. Projection Buffers (Zero-Allocation Architecture)
# ------------------------------------------------------------------------------

"""
    DiscreteBuffer{T}
Workspace for Log-Sum-Exp summation at roots of unity (Real-valued).
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

# Global Pre-allocated Workspaces
const _WS_F64 = DiscreteBuffer{Float64}(8192)
const _WS_C64 = AnalyticBuffer{Float64}(8192)

# ------------------------------------------------------------------------------
# 4. Display Logic
# ------------------------------------------------------------------------------

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