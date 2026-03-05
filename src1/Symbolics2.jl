#Symbolics.jl

# ============================================================
# Symbolic Engine in Cyclo Monomials
# ============================================================

#quantum factorials 
function qfactorial_symb(n::Int)
    if n == 0
        return CycloMonomial(1, 0, Dict{Int, Int}())
    end

    exps = Dict{Int, Int}()
    for d in 2:n
        power = div(n, d) 
        if power > 0
            exps[d] = power
        end
    end
    
    # Power of z = q^{1/2}. q^{-n(n-1)/4} = z^{-n(n-1)/2}
    z_pow = - (n * (n - 1)) ÷ 2
    
    return CycloMonomial(1, z_pow, exps)
end

# triangle coefficient 
function qdelta2_symb(j1, j2, j3)
    num = qfactorial_symb(Int(j1+j2-j3)) * qfactorial_symb(Int(j1-j2+j3)) * qfactorial_symb(Int(-j1+j2+j3))
    den = qfactorial_symb(Int(j1+j2+j3+1))
    return num / den
end

function qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    return qdelta2_symb(j1, j2, j3) * qdelta2_symb(j1, j5, j6) * qdelta2_symb(j2, j4, j6) * qdelta2_symb(j3, j4, j5)
end

# racah summand terms  
function q6jsummand_symb(z, α1, α2, α3, α4, β1, β2, β3)
    num = qfactorial_symb(z+1)
    den = qfactorial_symb(z-α1) * qfactorial_symb(z-α2) * qfactorial_symb(z-α3) *
          qfactorial_symb(z-α4) * qfactorial_symb(β1-z) * qfactorial_symb(β2-z) *
          qfactorial_symb(β3-z) 
    
    res = num / den
    if isodd(z)
        return CycloMonomial(-res.sign, res.z_pow, res.exps)
    else
        return res
    end
end

# vector of racah summand  
function q6jseries_symb(j1, j2, j3, j4, j5, j6)::Vector{CycloMonomial}
    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    S_z = CycloMonomial[]
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3)
    @inbounds for z in zrange
        push!(S_z, q6jsummand_symb(z, α1, α2, α3, α4, β1, β2, β3)) 
    end
    return S_z
end


function qracah6j_generic(j1, j2, j3, j4, j5, j6)
    if !δtet(j1, j2, j3, j4, j5, j6)
        return CycloMonomial(0,0,Dict{Int,Int}())
    end

    Tc2 = qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    series = q6jseries_symb(j1, j2, j3, j4, j5, j6)
    return Generic6jResult(Tc2,series)
end



#------- quantum 3j symbols -----

function q3jsummand_symb(z,α1, α2, β1, β2, β3)
    num = CycloMonomial(1,0,Dict{Int,Int}())
    den = qfactorial_symb(z) * qfactorial_symb(α1 + z) * qfactorial_symb(α2 + z) *
        qfactorial_symb(β1 - z) * qfactorial_symb(β2 - z) * qfactorial_symb(β3 - z)
    res = num / den #faster way?
    if isodd(Int(z + α1 - α2)) # include extra sign: (-1)^{j1-j2-m3}
        return CycloMonomial(-res.sign, res.z_pow, res.exps) 
    else
        return res
    end
end


function q3jseries_symb(j1,j2,j3,m1,m2)

    α1 = Int(j3 - j2 + m1) 
    α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 -j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    S_z = CycloMonomial[]
    zrange = max(α1,α2,zero(α1)):min(β1,β2,β3)
    @inbounds for z in zrange
        push!(S_z, q3jsummand_symb(z, α1, α2, β1, β2, β3)) 
    end

    return S_z
end

function qracah3j_generic(j1, j2, j3, m1,m2, m3 = -m1-m2)

    if !δ(j1, j2, j3) || !iszero(m1+m2+m3)
        return CycloMonomial(0,0,Dict{Int,Int}())
    end
    #prefactor_square
    pref_sq = qdelta2_symb(j1,j2,j3) * qfactorial_symb(j1+m1) * 
        qfactorial_symb(j1-m1) * qfactorial_symb(j2+m2) * 
        qfactorial_symb(j2-m2) * qfactorial_symb(j3-m1-m2) * 
        qfactorial_symb(j3+m1+m2)
    series = q3jseries_symb(j1,j2,j3,m1,m2)
    return Generic6jResult(pref_sq,series)
end






#====== TOFIX Below: ================#

# evaluation at q=1 

function phi_at_one(d::Int)
    factors = collect(keys(factor(d)))
    if length(factors) == 1
        return factors[1] # d is a power of a single prime p
    else
        return 1 # d has multiple distinct prime factors
    end
end

"""
    evaluate_classical(M::CycloMonomial) -> Rational{BigInt}

Evaluates a CycloMonomial at q=1. 
Uses the property: Φ_d(1) = p if d = p^k, else 1 (for d > 1).
"""
function evaluate_classical(M::CycloMonomial)
    M.sign == 0 && return big(0)
    
    res = BigInt(M.sign)
    # At q=1, z = 1, so z^z_pow = 1
    
    for (d, e) in M.exps
        d <= 1 && continue # Φ_1(1) is 0, handled by admissibility
        
        # Calculate Φ_d(1)
        val = phi_at_one(d)
        if e > 0
            res *= BigInt(val)^e
        elseif e < 0
            # Use Rational for division to maintain exactness
            return Rational(res) / (Rational(phi_at_one(d))^abs(e))
        end
    end
    return res
end

function qracah6j_classical(j1, j2, j3, j4, j5, j6)
    # if !δtet(j1, j2, j3, j4, j5, j6)
    #     return BigInt(0)
    # end
    q6j = qracah6j_generic(j1, j2, j3, j4, j5, j6)
    sumz = sum(evaluate_classical.(q6j.series) )
    return sqrt(evaluate_classical(q6j.pref_sq)) * sumz
end