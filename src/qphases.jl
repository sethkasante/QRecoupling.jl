# ---------------------------
#   -- QPhase Struct -- 
# to handle phase functions 
# ---------------------------

"""
    QPhase

Represents the formal monomial `sign * q^(q_pow)`, with `sign` in `(-1, 0, 1)`.
It retains the chosen q branch, not a level. Rational powers act on the formal
exponent and the real sign; they need not equal a principal power after evaluation.
"""
struct QPhase
    sign::Int8
    q_pow::Rational{Int}
    function QPhase(sign, q_pow)
        sign in (-1,0,1) || throw(ArgumentError("QPhase sign must be -1, 0 or 1"))
        new(Int8(sign), iszero(sign) ? 0//1 : Rational{Int}(q_pow))
    end
end

# --- Basic identities ---
Base.iszero(p::QPhase) = p.sign == 0
Base.one(::Type{QPhase}) = QPhase(Int8(1), 0//1)
Base.zero(::Type{QPhase}) = QPhase(Int8(0), 0//1)
Base.one(::QPhase) = one(QPhase)
Base.zero(::QPhase) = zero(QPhase)
Base.isone(p::QPhase) = p.sign == 1 && iszero(p.q_pow)
Base.sign(p::QPhase) = p.sign
Base.copy(p::QPhase) = QPhase(p.sign, p.q_pow)

Base.:(==)(a::QPhase, b::QPhase) = (iszero(a) && iszero(b)) || (a.sign == b.sign && a.q_pow == b.q_pow)
Base.hash(p::QPhase, h::UInt) = hash((p.sign,p.q_pow),h)

# --- unary operators ---
Base.:-(p::QPhase) = QPhase(Int8(-p.sign), p.q_pow)
Base.:+(p::QPhase) = p

# --- multiplication ---
function Base.:*(a::QPhase, b::QPhase)
    (iszero(a) || iszero(b)) && return zero(QPhase)
    return QPhase(a.sign * b.sign, a.q_pow + b.q_pow)
end

function Base.inv(p::QPhase)
    iszero(p) && throw(DivideError())
    return QPhase(p.sign, -p.q_pow)
end

Base.:/(a::QPhase, b::QPhase) = a * inv(b)

function Base.:^(p::QPhase, n::Integer)
    if iszero(p)
        n < 0 && throw(DivideError())
        return n == 0 ? one(QPhase) : zero(QPhase)
    end
    new_sign = iseven(n) ? Int8(1) : p.sign
    return QPhase(new_sign, p.q_pow * n)
end

function Base.:^(p::QPhase, r::Rational)
    if iszero(p)
        r < 0 && throw(DivideError())
        return r == 0 ? one(QPhase) : zero(QPhase)
    end
    # Prevent complex numbers natively appearing from fractional powers of negative signs
    if p.sign == -1 && iseven(denominator(r))
        throw(DomainError(r, "Cannot take fractional power with an even denominator of a negative QPhase exactly within real signs."))
    end
    new_sign = (p.sign == -1 && isodd(numerator(r))) ? Int8(-1) : Int8(1)
    return QPhase(new_sign, p.q_pow * r)
end

# --- print repl ---
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

"Display a non-reciprocal phase over x without discarding its choice of q."
function Base.show(io::IO, ::MIME"text/plain", p::QPhase)
    println(io,"Exact phase over x = q + q⁻¹")
    if iszero(p) || iszero(p.q_pow)
        print(io,"  = "); show(io,p)
        return
    end
    println(io,"  = ",p.sign<0 ? "−" : "","u^(",p.q_pow,")")
    print(io,"  u² − x·u + 1 = 0, with u = q; fractional powers retain the q-phase branch")
end

# --- QPhase * CyclotomicMonomial ---
function Base.:*(phase::QPhase, m::CyclotomicMonomial)
    (iszero(phase) || iszero(m)) && return ZERO_MONOMIAL 
    
    if denominator(phase.q_pow) == 1
        # It's a clean integer! Safe to absorb.
        return CyclotomicMonomial(phase.sign * m.sign,
            m.q_pow + Int(numerator(phase.q_pow)),
            m.phi_exps,
            m.max_d     
        )
    else
        throw(ArgumentError("Cannot absorb fractional QPhase (q^$(phase.q_pow)) into a CyclotomicMonomial; retain the phase separately."))
    end
end

Base.:*(m::CyclotomicMonomial, phase::QPhase) = phase * m
Base.:/(m::CyclotomicMonomial, phase::QPhase) = m * inv(phase)

# --- QPhase * CompositeExactResult ---
# Only integer powers of q = ζ lie in ℚ(ζ_{2(k+2)}); they multiply every factor exactly.
function Base.:*(phase::QPhase, comp::CompositeExactResult{T}) where T
    (iszero(phase) || isempty(comp.terms)) && return zero(comp)
    denominator(phase.q_pow) == 1 || throw(ArgumentError(
        "q^($(phase.q_pow)) is not an element of ℚ(ζ$(to_subscript(2 * (comp.k + 2)))); only integer powers of q can multiply an exact result."))
    ζ = gen(parent(first(values(comp.terms))))
    f = Int(phase.sign) * ζ^Int(numerator(phase.q_pow))
    return CompositeExactResult{T}(comp.k, Dict{CyclotomicMonomial, T}(rad => f * val for (rad, val) in comp.terms))
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
    # Compute max_d 
    d_max = isempty(phi_exps) ? 1 : maximum(first, phi_exps)
    return CyclotomicMonomial(Int8(sign), q_pow, phi_exps, d_max)
end



# apis 

"""
    rmatrix(j1::Spin, j2::Spin, j3::Spin; k=nothing, q=nothing, exact::Bool=false, T::Type=ComplexF64)

Returns the R-matrix phase. 
Formula: R = (-1)^{j_1 + j_2 - j_3} q^{j_3(j_3+1) - j_1(j_1+1) - j_2(j_2+1)}. 
The default is the classical value. `rmatrix(Symbolic(), ...)` returns a `QPhase`.
At a level, `exact=true` retains the exact `QPhase` representation.
Numerical level phases are complex even when `T` or `Level(k; T=...)` names a real type.
"""
function rmatrix(j1::Spin, j2::Spin, j3::Spin; 
                 k=nothing, q=nothing, exact::Bool=false, T::Type=ComplexF64)
    
    q = _evaluation_q(k,q,exact)
    k isa AbstractVector && return [rmatrix(j1,j2,j3;k=kk,exact=exact,T=T) for kk in k]
    _check_phase_q(q)
    E = _phase_eltype(T,k,q)
    J1, J2, J3 = doubled(j1, j2, j3)
    
    # check admissibility 
    if !_δ(J1, J2, J3)
        return exact ? (isnothing(k) ? 0 : zero(QPhase)) : E(0)
    end
    
    if !isnothing(k) && !_qδ(J1, J2, J3, Int(k))
        return exact ? zero(QPhase) : E(0)
    end

    p = (J3*(J3+2) - J1*(J1+2) - J2*(J2+2)) ÷ 2
    s = iseven((J1 + J2 - J3) ÷ 2) ? Int8(1) : Int8(-1)
    
    _is_classical(q) && return exact ? Int(s) : T(s)

    # exact q-phase at a level
    if exact
        return QPhase(s, p // 2)
    end
    
    # level k
    if !isnothing(k)
        return E(s * _level_phase(p,Int(k),real(E)))
    end
    
    # generic q
    if !isnothing(q)
        v = s * qhalfpow(_phase_parameter(q,T), p)
        return E(v)
    end
end

@inline function _check_phase_q(q)
    q === nothing || (q isa Number && isfinite(q) && !iszero(q)) ||
        throw(DomainError(q,"q must be finite and nonzero"))
    return nothing
end

@inline function _phase_eltype(::Type{T},k,q) where {T}
    k === nothing || return Complex{real(float(T))}
    (q === nothing || _is_classical(q)) && return T
    Q = typeof(float(q))
    return promote_type(T, q isa Real && q < 0 ? Complex{real(Q)} : Q)
end

"q^(p/2) at a level, reducing large integer exponents before floating-point division."
@inline function _level_phase(p::Int,k::Int,::Type{R}) where {R}
    h = k+2
    r = -2h <= p <= 2h ? p : rem(p,4h)
    return cispi(R(r)/R(2h))
end

"Preserve the input precision, widening the parameter before a BigFloat phase calculation."
function _phase_parameter(q, ::Type{T}) where {T}
    real(T) === BigFloat && return q isa Real ? BigFloat(q) : Complex{BigFloat}(q)
    return q
end

"""
    qhalfpow(q, p::Integer)

`q^{p/2}` for an integer `p` — the half-integral power every braiding phase in the package ends in.

Positive real q uses a real power, avoiding unnecessary complex logarithms.
On the negative real axis, including an imaginary part of either signed zero,
the package takes arg(q)=π: the magnitude is a real power and i^p is applied
by exact sign changes and component swaps. Away from that axis the principal
complex power is used.
"""
function qhalfpow(q::Number, p::Integer)
    pp = Int(p)
    if real(q) < 0 && iszero(imag(q))
        iseven(pp) && return complex(float(real(q))^(pp ÷ 2))
        r = float(-real(q))
        v = r^(pp / 2)
        z = zero(v)
        return mod(pp,4) == 1 ? complex(z,v) : complex(z,-v)
    end
    if (q isa Real || iszero(imag(q))) && real(q) > 0
        r = float(real(q))
        return iseven(pp) ? r^(pp ÷ 2) : r^(pp / 2)
    end
    return complex(float(q))^(pp / 2)
end
