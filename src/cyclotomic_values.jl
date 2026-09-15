
# ---------------------------------------------------------------------------------
#   Cyclotomic values at q = e^{iπ/h}, h = k + 2
#
#   For x = q² (a primitive h-th root of unity) and d ≠ h:
#     Φ_d(x)   = s_d · |Φ_d(x)| · e^{iπ ω_d / 2h}
#     |Φ_d(x)| = Λ̃(d/h) · Π_{e | d, h ∤ e} (2|sin(πe/h)|)^{μ(d/e)}
#     s_d      = (-1)^{Σ_{e | d, h ∤ e} ⌊e/h⌋ + [h | d] φ(d/h)}
#     ω_d      = 2φ(d) + h·[d = 1]
#   with Λ̃(m) = ℓ if m is a power of the prime ℓ and 1 otherwise (Λ̃ = 1 when h ∤ d).
#   Φ_h(x) = 0: its exponent is tracked as a valuation and never evaluated.
#
#   Balanced convention: for d ≥ 2, Ψ_d(q) = q^{-φ(d)} Φ_d(q²) = s_d |Φ_d(x)| is real on |q| = 1,
#   so a monomial is σ q^P Π Ψ_d^{e_d} with P = balanced_phase(m), and
#   √(q^P Π Ψ_d^{e_d}) := q^{P/2} √(Π Ψ_d^{e_d}), the branch continuous from q = 1.
# ---------------------------------------------------------------------------------


"Möbius function and Euler totient on 1..D."
function arith_sieve(D::Int)
    μ = zeros(Int, D)
    D >= 1 && (μ[1] = 1)
    for i in 1:D, j in 2i:i:D
        μ[j] -= μ[i]
    end
    φ = collect(1:D)
    for p in 2:D
        φ[p] == p || continue      # untouched so far, so p is prime
        for m in p:p:D
            φ[m] -= φ[m] ÷ p
        end
    end
    return μ, φ
end

"Divisor lists for 1..D."
function divisor_lists(D::Int)
    divs = [Int[] for _ in 1:D]
    for d in 1:D, n in d:d:D
        push!(divs[n], d)
    end
    return divs
end

"Λ̃(m) = ℓ if m = ℓ^a for a prime ℓ and a ≥ 1, otherwise 1."
function lambda_tilde(m::Int)
    m < 2 && return 1
    p = 2
    while p * p <= m && m % p != 0
        p += 1
    end
    m % p == 0 || return m          # no factor up to √m: m is prime
    while m % p == 0
        m ÷= p
    end
    return m == 1 ? p : 1
end

"Euler totient of a single n ≥ 1."
function _totient(n::Int)
    n < 1 && throw(DomainError(n, "totient requires n ≥ 1"))
    r = n; m = n; p = 2
    while p * p <= m
        if m % p == 0
            while m % p == 0
                m ÷= p
            end
            r -= r ÷ p
        end
        p += 1
    end
    m > 1 && (r -= r ÷ m)
    return r
end

"""
    balanced_phase(m::CyclotomicMonomial) -> Int

The exponent P with m(q) = σ q^P Π Ψ_d(q)^{e_d}, where Ψ_d(q) = q^{-φ(d)} Φ_d(q²).
"""
function balanced_phase(m::CyclotomicMonomial)
    P = m.q_pow
    for (d, e) in m.phi_exps
        P += e * _totient(d)
    end
    return P
end


# ---- Numeric table at q = e^{iπ/h} ----

"""
    RootOfUnityTable{T}

Φ_d(q²) at q = e^{iπ/h} for d = 1..D as a log-magnitude and a phase in units of π/2h (mod 4h).
The phase includes the sign s_d (a sign of -1 adds 2h). The entry at d = h is a placeholder;
Φ_h is handled as a valuation. `φ` is the Euler totient, used for balanced phases.
"""
struct RootOfUnityTable{T}
    h::Int
    logmag::Vector{T}
    ω::Vector{Int}
    φ::Vector{Int}
end

function build_rou_table(D::Int, k::Int, ::Type{T}) where {T}
    h = k + 2
    D = max(D, h)
    μ, φ = arith_sieve(D)
    logmag = zeros(T, D)
    par = zeros(Int, D)
    for e in 1:D
        e % h == 0 && continue
        le = log(2 * abs(sinpi(T(e % h) / h)))
        fl = e ÷ h
        for n in e:e:D
            mu = μ[n÷e]
            mu == 0 && continue
            logmag[n] += mu * le
            par[n] += fl
        end
    end
    ω = Vector{Int}(undef, D)
    for d in 1:D
        if d == h
            logmag[d] = zero(T)
            par[d] = 0
        elseif d % h == 0
            logmag[d] += log(T(lambda_tilde(d ÷ h)))
            par[d] += φ[d÷h]
        end
        ω[d] = mod(2φ[d] + (d == 1 ? h : 0) + (isodd(par[d]) ? 2h : 0), 4h)
    end
    return RootOfUnityTable{T}(h, logmag, ω, φ)
end

"Balanced phase using the totient table (all d must be covered by `tab`)."
function balanced_phase(m::CyclotomicMonomial, tab::RootOfUnityTable)
    P = m.q_pow
    @inbounds for (d, e) in m.phi_exps
        P += e * tab.φ[d]
    end
    return P
end

"""
    mono_at_root(m, tab) -> (logabs, σ, w, v)

m(q) = σ · exp(logabs) · e^{iπ w / 2h} · Φ_h(q²)^v at q = e^{iπ/h}, with w reduced mod 4h and
σ = m.sign (the signs of the cyclotomic factors are part of w).
"""
@inline function mono_at_root(m::CyclotomicMonomial, tab::RootOfUnityTable{T}) where {T}
    h = tab.h
    lm = zero(T)
    w = 2 * m.q_pow
    v = 0
    @inbounds for (d, e) in m.phi_exps
        if d == h
            v += e
        else
            lm += e * tab.logmag[d]
            w += e * tab.ω[d]
        end
    end
    return lm, Int(m.sign), mod(w, 4h), v
end

"x · e^{iπ w / 2h}, kept real when the phase is 0 or π."
@inline function apply_phase(x::Number, w::Int, h::Int)
    w = mod(w, 4h)
    w == 0 && return x
    w == 2h && return -x
    return x * cispi(real(typeof(float(x)))(w) / (2h))
end
