# ---------------------------
#   -- QPhase Struct -- 
# to handle phase functions 
# ---------------------------

"""
    QPhase

Represents an exact fractional q-phase: `sign * q^(q_pow)`.
"""
struct QPhase
    sign::Int8
    q_pow::Rational{Int}
end

# --- Basic Identities ---
Base.iszero(p::QPhase) = p.sign == 0
Base.one(::Type{QPhase}) = QPhase(Int8(1), 0//1)
Base.zero(::Type{QPhase}) = QPhase(Int8(0), 0//1)
Base.sign(p::QPhase) = p.sign
Base.copy(p::QPhase) = QPhase(p.sign, p.q_pow)

Base.:(==)(a::QPhase, b::QPhase) = (iszero(a) && iszero(b)) || (a.sign == b.sign && a.q_pow == b.q_pow)

# --- Unary Operators ---
Base.:-(p::QPhase) = QPhase(Int8(-p.sign), p.q_pow)
Base.:+(p::QPhase) = p

# --- Multiplicative Group ---
function Base.:*(a::QPhase, b::QPhase)
    (iszero(a) || iszero(b)) && return zero(QPhase)
    return QPhase(a.sign * b.sign, a.q_pow + b.q_pow)
end

function Base.inv(p::QPhase)
    iszero(p) && throw(DivideError())
    return QPhase(p.sign, -p.q_pow)
end

Base.:/(a::QPhase, b::QPhase) = a * inv(b)

# --- Exponentiation ---
function Base.:^(p::QPhase, n::Integer)
    iszero(p) && return n == 0 ? one(QPhase) : zero(QPhase)
    new_sign = iseven(n) ? Int8(1) : p.sign
    return QPhase(new_sign, p.q_pow * n)
end

function Base.:^(p::QPhase, r::Rational)
    iszero(p) && return r == 0 ? one(QPhase) : zero(QPhase)
    # Prevent complex numbers natively appearing from fractional powers of negative signs
    if p.sign == -1 && iseven(denominator(r))
        throw(DomainError(r, "Cannot take fractional power with an even denominator of a negative QPhase exactly within real signs."))
    end
    new_sign = (p.sign == -1 && isodd(numerator(r))) ? Int8(-1) : Int8(1)
    return QPhase(new_sign, p.q_pow * r)
end

# --- Printing ---
function Base.show(io::IO, p::QPhase)
    iszero(p) && return print(io, "0")
    s_str = p.sign == -1 ? "-" : ""
    pow = p.q_pow
    if pow == 0
        print(io, p.sign == -1 ? "-1" : "1")
    elseif pow == 1
        print(io, "$(s_str)q")
    else
        print(io, "$(s_str)q^($pow)")
    end
end

# --- QPhase * CyclotomicMonomial ---
function Base.:*(phase::QPhase, m::CyclotomicMonomial)
    # Assuming ZERO_MONOMIAL is defined in your constants
    (iszero(phase) || iszero(m)) && return ZERO_MONOMIAL 
    
    if denominator(phase.q_pow) == 1
        # It's a clean integer! Safe to absorb.
        return CyclotomicMonomial(
            phase.sign * m.sign,
            m.q_pow + Int(numerator(phase.q_pow)),
            m.phi_exps,
            m.max_d     # <-- Passed safely through!
        )
    else
        throw(ArgumentError("Cannot absorb fractional QPhase (q^$(phase.q_pow)) into a discrete CyclotomicMonomial. Multiply fractional phases against the parent CompositeExactResult instead."))
    end
end

Base.:*(m::CyclotomicMonomial, phase::QPhase) = phase * m
Base.:/(m::CyclotomicMonomial, phase::QPhase) = m * inv(phase)

# --- QPhase * CompositeExactResult ---
function Base.:*(phase::QPhase, comp::CompositeExactResult{T}) where T
    # If phase is zero, return an empty dictionary composite
    iszero(phase) && return CompositeExactResult{T}(comp.k, comp.global_phase, comp.radical, Dict{CyclotomicMonomial, T}())
    
    # Shift the global phase tracker
    new_phase = comp.global_phase + phase.q_pow
    
    # Flip the signs in the evaluated Nemo dictionaries if necessary
    new_terms = Dict{CyclotomicMonomial, T}()
    for (rad, val) in comp.terms
        new_terms[rad] = phase.sign == 1 ? val : -val
    end
    
    return CompositeExactResult{T}(comp.k, new_phase, comp.radical, new_terms)
end

Base.:*(comp::CompositeExactResult, phase::QPhase) = phase * comp
Base.:/(comp::CompositeExactResult, phase::QPhase) = comp * inv(phase)


# --- Helpers ---- 

"""
    q_phase(pow; sign=1)

Creates a pure `QPhase` object with a fractional/rational power of q.
"""
q_phase(pow; sign=1) = QPhase(Int8(sign), Rational{Int}(pow))

"""
    q_phi(d::Int, e::Int=1)

Creates an isolated `CyclotomicMonomial` for the polynomial Φ_d(q^2)^e.
"""
q_phi(d::Int, e::Int=1) = CyclotomicMonomial(Int8(1), 0, [d => e], d)

"""
    q_mono(q_pow::Int; sign=1, phi_exps=Pair{Int,Int}[])

Creates a full `CyclotomicMonomial` requiring a strict integer power for q.
Automatically calculates `max_d` from the provided polynomial exponents.
"""
function q_mono(q_pow::Int; sign=1, phi_exps=Pair{Int,Int}[])
    # Compute max_d dynamically so the CycloBuffer doesn't overflow later
    d_max = isempty(phi_exps) ? 1 : maximum(first, phi_exps)
    return CyclotomicMonomial(Int8(sign), q_pow, phi_exps, d_max)
end



# apis 

"""
    rmatrix(j1::Spin, j2::Spin, j3::Spin; k=nothing, q=nothing, exact::Bool=false, T::Type=ComplexF64)

Returns the R-matrix phase. 
Formula: R = (-1)^{j_1 + j_2 - j_3} q^{j_3(j_3+1) - j_1(j_1+1) - j_2(j_2+1)}. 
If no evaluation target is provided, returns the exact `QPhase`.
"""
function rmatrix(j1::Spin, j2::Spin, j3::Spin; 
                 k=nothing, q=nothing, exact::Bool=false, T::Type=ComplexF64)
    
    J1, J2, J3 = doubled(j1, j2, j3)
    
    # check admissibility 
    if !_δ(J1, J2, J3)
        return (exact || (isnothing(k) && isnothing(q))) ? zero(QPhase) : T(0)
    end
    
    # if k is provided
    if !isnothing(k) && !_qδ(J1, J2, J3, k)
        return (exact || isnothing(q)) ? zero(QPhase) : T(0)
    end

    # Phase formula variables
    p = (J3*(J3+2) - J1*(J1+2) - J2*(J2+2)) ÷ 2
    s = iseven((J1 + J2 - J3) ÷ 2) ? Int8(1) : Int8(-1)
    
    # --- exact q phase ---
    if (isnothing(k) && isnothing(q)) || exact
        return QPhase(s, p // 2)
    end
    
    # --- classical limit ---
    if !isnothing(q) && (q == 1 || q == 1.0)
        return T(s)
    end
    
    # --- Numeric (level k) ---
    if !isnothing(k)
        h = k + 2
        phase_angle = (pi * p) / (2 * h)
        return T(s * cis(phase_angle))
    end
    
    # --- generic q ---
    if !isnothing(q)
        q_C = complex(float(q))
        return T(s * exp((p / 2) * log(q_C)))
    end
end