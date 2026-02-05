module Quantum6jSU2k

using LRUCache

export SU2kModel, q6j

# ===============================
# Model
# ===============================

struct SU2kModel
    k::Int
    prec::Int
end

const LOGQ_CACHE = LRU{Tuple{Int,Int}, Vector{BigFloat}}(maxsize=16)
const SIXJ_CACHE = LRU{Tuple{Int,Int,NTuple{6,Int}}, BigFloat}(maxsize=512)

# ===============================
# Utilities
# ===============================

@inline spin2(j::Real) = begin
    J = round(Int, 2j)
    abs(2j - J) < 1e-12 || throw(ArgumentError("spin $j not half-integer"))
    J
end

@inline admissible(Ja,Jb,Jc,k) =
    (Ja + Jb + Jc ≤ k) &&
    ((Ja + Jb + Jc) % 2 == 0) &&
    (Ja ≤ Jb + Jc) && (Jb ≤ Ja + Jc) && (Jc ≤ Ja + Jb)

# ===============================
# Log q-factorials
# ===============================

function logqfact_table(k::Int, prec::Int)
    key = (k,prec)
    haskey(LOGQ_CACHE,key) && return LOGQ_CACHE[key]

    setprecision(prec) do
        N = k + 3
        tab = Vector{BigFloat}(undef, N+1)
        tab[1] = 0  # log([0]_q!)

        denom = sin(big(pi)/(k+2))
        for n in 1:N
            qn = sin(big(pi)*n/(k+2)) / denom
            tab[n+1] = tab[n] + log(abs(qn))
        end

        LOGQ_CACHE[key] = tab
        return tab
    end
end

@inline logqfact(tab,n) = tab[n+1]

# ===============================
# Triangle prefactor
# ===============================

function logΔ(Ja,Jb,Jc,k,tab)
    admissible(Ja,Jb,Jc,k) || return -Inf
    a = (Ja+Jb-Jc)÷2
    b = (Ja-Jb+Jc)÷2
    c = (-Ja+Jb+Jc)÷2
    s = (Ja+Jb+Jc)÷2
    (logqfact(tab,a) + logqfact(tab,b) + logqfact(tab,c) - logqfact(tab,s+1)) / 2
end

# ===============================
# Quantum 6j
# ===============================

function q6j(j1,j2,j3,j4,j5,j6; k::Int, prec::Int=256)
    model = SU2kModel(k,prec)

    J = ntuple(i->spin2((j1,j2,j3,j4,j5,j6)[i]), 6)
    key = (k,prec,J)
    haskey(SIXJ_CACHE,key) && return SIXJ_CACHE[key]

    J1,J2,J3,J4,J5,J6 = J

    if !( admissible(J1,J2,J3,k) &&
          admissible(J1,J5,J6,k) &&
          admissible(J2,J4,J6,k) &&
          admissible(J3,J4,J5,k) )
        return zero(BigFloat)
    end

    tab = logqfact_table(k,prec)

    A = (
        (J1+J2+J3)÷2,
        (J1+J5+J6)÷2,
        (J2+J4+J6)÷2,
        (J3+J4+J5)÷2
    )
    B = (
        (J1+J2+J4+J5)÷2,
        (J1+J3+J4+J6)÷2,
        (J2+J3+J5+J6)÷2
    )

    zmin = maximum(A)
    zmax = minimum(B)
    zmin > zmax && return zero(BigFloat)

    log_pref =
        logΔ(J1,J2,J3,k,tab) +
        logΔ(J1,J5,J6,k,tab) +
        logΔ(J2,J4,J6,k,tab) +
        logΔ(J3,J4,J5,k,tab)

    setprecision(prec) do
        sumz = zero(BigFloat)
        for z in zmin:zmax
            lognum = logqfact(tab,z+1)
            logden = sum(logqfact(tab,z-a) for a in A) +
                     sum(logqfact(tab,b-z) for b in B)
            sumz += (-1)^z * exp(lognum - logden)
        end

        val = exp(log_pref) * sumz
        SIXJ_CACHE[key] = val
        return val
    end
end

end # module