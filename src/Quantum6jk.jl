

#qsu2k_v2.jl
module QuantumSU2k


export qinteger, qfactorial, 
        log_qfactorial, log_qΔ, qδtet, 
        qΔcoeffs, q6jsummand, qRacahsum, _qracah6j




function qinteger(n::Int,k::Int)
    (sinpi(n/(k+2))) / sinpi(1/(k+2))
end 

function log_qfactorial(n::Int,k::Int)
    n > k +1 && return throw(ArgumentError("Poles"))
    (n == 0 || n == 1) && return 0.0
    return sum(log(qinteger(i,k)) for i in 1:n)
end 

function log_qΔ2(j1,j2,j3,k)
    log_qfactorial(j1 + j2 - j3,k) + log_qfactorial(j1 - j2 + j3,k) + 
        log_qfactorial(-j1 + j2 + j3,k) - log_qfactorial(j1 + j2 + j3 + 1,k)
end

function log_qΔ(j1,j2,j3,k)
    !qδ(j1, j2, j3, k) && return 0.0
    (log_qfactorial(j1 + j2 - j3,k) + log_qfactorial(j1 - j2 + j3,k) + 
        log_qfactorial(-j1 + j2 + j3,k) - log_qfactorial(j1 + j2 + j3 + 1,k))/2
end

function qΔcoeffs(j1,j2,j3,j4,j5,j6,k)
    log_terms = log_qΔ(j1, j2, j3, k) + log_qΔ(j1, j5, j6, k) + log_qΔ(j2, j4, j6, k) +
        log_qΔ(j3, j4, j5, k)
    return exp(log_terms)
end

function q6jsummand(z,α1,α2,α3,α4,β1,β2,β3,k;prec=256)
    setprecision(prec) do
        log_term = log_qfactorial(z+1,k) - ( log_qfactorial(z-α1,k) + 
                log_qfactorial(z-α2,k) + log_qfactorial(z-α3,k) +
                log_qfactorial(z-α4,k) + log_qfactorial(β1-z,k) + 
                log_qfactorial(β2-z,k) + log_qfactorial(β3-z,k) )
        # println(exp(log_num - log_den))
        return exp(log_term)
    end
end

function qRacahsum(α1,α2,α3,α4,β1,β2,β3,k;prec=256)
    #range of z values
    zrange = max(α1,α2,α3,α4):min(β1,β2,β3)
    # n = length(zrange)
    setprecision(prec) do
        sumz = big(0)
        for z in zrange
            sz = q6jsummand(z,α1,α2,α3,α4,β1,β2,β3,k) 
            sumz += iseven(z) ? sz : -sz 
        end
    return sumz
    end
end


function _qracah6j(j1, j2, j3, j4, j5, j6, k;prec=256)::Float64
    
    T = qΔcoeffs(j1,j2,j3,j4,j5,j6,k)
    
    #WignerSymbols reorders this for storage purposes 
    α1, α2, α3, α4 = j1 + j2 + j3, j1 + j5 + j6, j2 + j4 + j6, j3 + j4 + j5
    β1, β2, β3 = j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6

    S = qRacahsum(α1,α2,α3,α4,β1,β2,β3, k;prec=256)
    
    return T*S 
end


# classical and quantum triangle conditions at level k: 
# checks admissible triple 
"""
    δ(j1, j2, j3,k) -> ::Bool
    qδ(j1, j2, j3,k) -> ::Bool

Checks the triangle conditions `j3 ≤ j1 + j2`, `j1 ≤ j2 + j3`, `j2 ≤ j3 + j1` and extra condition `j1 + j2 + j3 ≤ k` for qδ.
"""
δ(j1, j2, j3) = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)

# (j3<= j1+j2) && (j1 <= j2 + j3) && (j2 <= j3 + j1) 

qδ(j1, j2, j3, k) = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 


δtet(j1, j2, j3, j4, j5, j6) = δ(j1, j2, j3) && δ(j1, j5, j6) && δ(j2, j4, j6) && δ(j3, j4, j5)

qδtet(j1, j2, j3, j4, j5, j6, k) = qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)

#

end # module Quantum6j
