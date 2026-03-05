# src/Symbolics.jl

# ============================================================
# Symbolic Engine in Cyclo Monomials
# ============================================================

# Quantum factorials 
function qfactorial_symb(n::Int)
    if n <= 1
        return CycloMonomial(1, 0, Int[])
    end

    # Pre-allocate exactly n elements. exps[d] is the exponent for Φ_d
    exps = zeros(Int, n)
    @inbounds for d in 2:n
        power = div(n, d) 
        if power > 0
            exps[d] = power
        end
    end
    
    # Power of z = q^{1/2}. q^{-n(n-1)/4} = z^{-n(n-1)/2}
    z_pow = - (n * (n - 1)) ÷ 2
    
    return CycloMonomial(1, z_pow, exps)
end

# Triangle coefficient 
function qdelta2_symb(j1, j2, j3)
    num = qfactorial_symb(Int(j1+j2-j3)) * qfactorial_symb(Int(j1-j2+j3)) * qfactorial_symb(Int(-j1+j2+j3))
    den = qfactorial_symb(Int(j1+j2+j3+1))
    return num / den
end

function qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    return qdelta2_symb(j1, j2, j3) * qdelta2_symb(j1, j5, j6) * qdelta2_symb(j2, j4, j6) * qdelta2_symb(j3, j4, j5)
end

# Racah summand terms  
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

# Vector of Racah summand  
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
    if !δtet(j1, j2, j3, j4, j5, j6) # Assuming δtet is classical admissibility
        return GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[])
    end

    Tc2 = qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    series = q6jseries_symb(j1, j2, j3, j4, j5, j6)
    return GenericResult(Tc2, series)
end

# ======================
# Quantum 3j symbols 
# ======================

function q3jsummand_symb(z, α1, α2, β1, β2, β3)
    den = qfactorial_symb(z) * qfactorial_symb(α1 + z) * qfactorial_symb(α2 + z) *
          qfactorial_symb(β1 - z) * qfactorial_symb(β2 - z) * qfactorial_symb(β3 - z)
          
    res = inv(den) # inversion 
    
    if isodd(Int(z + α1 - α2)) # include extra sign: (-1)^{j1-j2-m3}
        return CycloMonomial(-res.sign, res.z_pow, res.exps) 
    else
        return res
    end
end

function q3jseries_symb(j1, j2, j3, m1, m2)
    α1 = Int(j3 - j2 + m1) 
    α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 - j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    S_z = CycloMonomial[]
    zrange = max(α1, α2, 0):min(β1, β2, β3)
    @inbounds for z in zrange
        push!(S_z, q3jsummand_symb(z, α1, α2, β1, β2, β3)) 
    end

    return S_z
end

function qracah3j_generic(j1, j2, j3, m1, m2, m3 = -m1-m2)
    if !δ(j1, j2, j3) || !iszero(m1 + m2 + m3)
        return GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[])
    end
    
    pref_sq = qdelta2_symb(j1, j2, j3) * qfactorial_symb(Int(j1+m1)) * qfactorial_symb(Int(j1-m1)) * qfactorial_symb(Int(j2+m2)) * qfactorial_symb(Int(j2-m2)) * qfactorial_symb(Int(j3-m1-m2)) * qfactorial_symb(Int(j3+m1+m2))
              
    series = q3jseries_symb(j1, j2, j3, m1, m2)
    return GenericResult(pref_sq, series)
end



# ============================================================
# Classical Evaluation at q = 1
# ============================================================

"""
    phi_at_one(d::Int) -> Int

Fast, allocation-free check. Φ_d(1) = p if d = p^k (a prime power), else 1.
"""
function phi_at_one(d::Int)
    d <= 1 && return 1
    p = 2
    # Find smallest prime factor
    while d % p != 0
        p += 1
        p * p > d && (p = d; break)
    end
    
    # Check if d is purely a power of p
    temp = d
    while temp % p == 0
        temp ÷= p
    end
    
    return temp == 1 ? p : 1
end

"""
    evaluate_classical(M::CycloMonomial) -> Rational{BigInt}

Evaluates a CycloMonomial at q=1. Returns an exact fraction to prevent precision loss.
"""
function evaluate_classical(M::CycloMonomial)
    M.sign == 0 && return big(0) // 1
    
    num = BigInt(1)
    den = BigInt(1)
    
    # Loop over the Vector indices (which represent the cyclotomic base d)
    @inbounds for d in 2:length(M.exps)
        e = M.exps[d]
        e == 0 && continue
        
        val = phi_at_one(d)
        val == 1 && continue # 1^e is 1, so skip
        
        if e > 0
            num *= BigInt(val)^e
        elseif e < 0
            den *= BigInt(val)^abs(e)
        end
    end
    
    final_num = M.sign == -1 ? -num : num
    return final_num // den
end

"""
    qracah6j_classical(j1, j2, j3, j4, j5, j6) -> Float64

Evaluates the classical Ponzano-Regge 6j symbol.
"""
function qracah6j_classical(j1, j2, j3, j4, j5, j6)
    q6j = qracah6j_generic(j1, j2, j3, j4, j5, j6)
    
    if q6j.pref_sq.sign == 0
        return 0.0
    end
    
    sumz = sum(evaluate_classical.(q6j.series))
    pref_sq_val = evaluate_classical(q6j.pref_sq)
    
    # Square root forces us into floats, but we deferred it to the absolute last step
    return Float64(sqrt(BigFloat(pref_sq_val)) * BigFloat(sumz))
end