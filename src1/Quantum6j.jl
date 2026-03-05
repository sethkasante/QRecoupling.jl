module Quantum6j

#cyclotomicfield.jl

using Hecke

export QSU2kField, 
        qroot, 
        qinteger, qfactorial, 
        qδ, qΔ2, 
        q6jseries

"""
    QSU2kField(k)

Cyclotomic field for SU(2)_k with ζ = exp(iπ/(k+2))
"""
struct QSU2kField
    k::Int
    n::Int
    K::AbsSimpleNumField 
    q::AbsSimpleNumFieldElem
end

function QSU2kField(k::Int)
    n = 2*(k + 2)
    K, ζ = cyclotomic_field(n, "ζ")
    return QSU2kField(k, n, K, ζ)
end

function cyclo_poly_value(Q::QSU2kField, d::Int)
    ζ = Q.q
    N = Q.n

    # ζ^(N/d) is a primitive d-th root if d | N
    if N % d == 0
        ζd = ζ^(div(N,d))
    else
        # embed Φ_d(q) symbolically anyway
        ζd = ζ^(div(N, gcd(N,d)))
    end

    prod(ζd^a - 1 for a in 1:d if gcd(a,d) == 1)
end

"""
    qroot(Q)

Return q = ζ = exp(iπ/(k+2)) : a primitive root of unity 
"""
qroot(Q::QSU2kField) = Q.q 


"""
    qinteger(Q, n)

Quantum integer [n]_q exactly in ℚ(ζ)
"""
function qinteger(Q::QSU2kField, n::Int)
    q = qroot(Q)
    return n == 0 ? 1 : (q^n - q^(-n)) // (q - q^(-1))
end

qfactorial(qf::QSU2kField, n::Int) = n == 0 ? 1 : prod(qinteger(qf,i) for i in 1:n)

#store qfactorials in a table. For a given k, store up to [k]! 


# q triangle condition: checks admissible triple 
"""
    qδ(qf, j1, j2, j3) -> ::Bool

Checks the triangle conditions `j3 ≤ j1 + j2`, `j1 ≤ j2 + j3`, `j2 ≤ j3 + j1` and j1+j2+j3 ≤ k.
"""
qδ(qf::QSU2kField, j1, j2, j3) = isinteger(j1+j2+j3) && (j3 <= j1 + j2) && (j1 <= j2 + j3) && (j2 <= j3 + j1) && (j1+j2+j3) <= qf.k 


# squared triangle coefficient
function qΔ2(qf::QSU2kField, j1, j2, j3)
    # check triangle conditions 
    qδ(qf, j1, j2, j3) || throw(DomainError("spins are not admissible"))
    n1 = qfactorial(qf, j1 + j2 - j3) 
    n2 = qfactorial(qf, j1 - j2 + j3) 
    n3 = qfactorial(qf, - j1 + j2 + j3) 
    num = n1*n2*n3
    den = qfactorial(qf, j1 + j2 + j3 + 1) 
    # result
    return num // den
end



function q6jseries(qf::QSU2kField, j1, j2, j3, j4, j5, j6) 
    α1, α2, α3, α4 = j1 + j2 + j3, j1 + j5 + j6, j2 + j4 + j6, j3 + j4 + j5
    β1, β2, β3 = j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3)
    
    sol = 0 
    for z in zrange
        num = iseven(z) ? qfactorial(qf,z+1) : -qfactorial(qf,z+1)
        den = qfactorial(qf,z-α1) * qfactorial(qf,z-α2) * 
                qfactorial(qf,z-α3) * qfactorial(qf,z-α4) * 
                qfactorial(qf,β1- z) * qfactorial(qf,β2- z) * qfactorial(qf,β3- z)
        sol += (num // den)
    end
    return sol 
end


end # module Quantum6j
