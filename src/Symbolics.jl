# src/Symbolics.jl

# ============================================================
# Core Allocation-Free In-Place Operations
# ============================================================

@inline function _ensure_capacity!(exps::Vector{Int}, n::Int)
    old_len = length(exps)
    if n > old_len
        resize!(exps, n)
        @inbounds for i in (old_len+1):n
            exps[i] = 0
        end
    end
    return nothing
end

function mul_qfact(M::CycloMonomial, n::Int, power::Int=1)
    n <= 1 && return M
    new_z_pow = M.z_pow + power * (-(n * (n - 1)) ÷ 2)
    _ensure_capacity!(M.exps, n)
    @inbounds for d in 2:n
        M.exps[d] += power * div(n, d)
    end
    return CycloMonomial(M.sign, new_z_pow, M.exps)
end

function mul_qint(M::CycloMonomial, n::Int, power::Int=1)
    n <= 1 && return M
    new_z_pow = M.z_pow + power * (1 - n)
    _ensure_capacity!(M.exps, n)
    @inbounds for d in 2:n
        if n % d == 0
            M.exps[d] += power
        end
    end
    return CycloMonomial(M.sign, new_z_pow, M.exps)
end

# function Base.inv(M::CycloMonomial)
#     M.sign == 0 && throw(DivideError())
#     return CycloMonomial(M.sign, -M.z_pow, -M.exps)
# end

# ============================================================
# Symbolic Engine (Single-Allocation Setup)
# ============================================================

function qfactorial_symb(n::Int)
    res = CycloMonomial(1, 0, Int[])
    return mul_qfact(res, n, 1)
end

function qdelta2_symb(j1::Spin, j2::Spin, j3::Spin)
    res = CycloMonomial(1, 0, Int[])
    res = mul_qfact(res, Int(j1+j2-j3), 1)
    res = mul_qfact(res, Int(j1-j2+j3), 1)
    res = mul_qfact(res, Int(-j1+j2+j3), 1)
    res = mul_qfact(res, Int(j1+j2+j3+1), -1)
    return res
end

function qtricoeff2_symb(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res = CycloMonomial(1, 0, Int[])
    res = mul_qfact(res, Int(j1+j2-j3), 1); res = mul_qfact(res, Int(j1-j2+j3), 1); res = mul_qfact(res, Int(-j1+j2+j3), 1); res = mul_qfact(res, Int(j1+j2+j3+1), -1)
    res = mul_qfact(res, Int(j1+j5-j6), 1); res = mul_qfact(res, Int(j1-j5+j6), 1); res = mul_qfact(res, Int(-j1+j5+j6), 1); res = mul_qfact(res, Int(j1+j5+j6+1), -1)
    res = mul_qfact(res, Int(j2+j4-j6), 1); res = mul_qfact(res, Int(j2-j4+j6), 1); res = mul_qfact(res, Int(-j2+j4+j6), 1); res = mul_qfact(res, Int(j2+j4+j6+1), -1)
    res = mul_qfact(res, Int(j3+j4-j5), 1); res = mul_qfact(res, Int(j3-j4+j5), 1); res = mul_qfact(res, Int(-j3+j4+j5), 1); res = mul_qfact(res, Int(j3+j4+j5+1), -1)
    return res
end

# ============================================================
# Racah Series (Hypergeometric Ratio Method)
# ============================================================

function q6jsummand_symb_first(z::Int, α1::Int, α2::Int, α3::Int, α4::Int, β1::Int, β2::Int, β3::Int)
    res = CycloMonomial(iseven(z) ? 1 : -1, 0, Int[])
    res = mul_qfact(res, z+1, 1)
    res = mul_qfact(res, z-α1, -1); res = mul_qfact(res, z-α2, -1); res = mul_qfact(res, z-α3, -1); res = mul_qfact(res, z-α4, -1)
    res = mul_qfact(res, β1-z, -1); res = mul_qfact(res, β2-z, -1); res = mul_qfact(res, β3-z, -1)
    return res
end

function q6jseries_symb(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::Vector{CycloMonomial}
    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    S_z = CycloMonomial[]
    z_min = max(α1, α2, α3, α4)
    z_max = min(β1, β2, β3)
    
    if z_min <= z_max
        term = q6jsummand_symb_first(z_min, α1, α2, α3, α4, β1, β2, β3)
        push!(S_z, term)

        for z in z_min+1 : z_max
            next_term = CycloMonomial(-term.sign, term.z_pow, copy(term.exps))
            
            next_term = mul_qint(next_term, z + 1, 1)
            next_term = mul_qint(next_term, β1 - z + 1, 1)
            next_term = mul_qint(next_term, β2 - z + 1, 1)
            next_term = mul_qint(next_term, β3 - z + 1, 1)
            
            next_term = mul_qint(next_term, z - α1, -1)
            next_term = mul_qint(next_term, z - α2, -1)
            next_term = mul_qint(next_term, z - α3, -1)
            next_term = mul_qint(next_term, z - α4, -1)
            
            push!(S_z, next_term)
            term = next_term
        end
    end
    return S_z
end

function qracah6j_generic(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    # if !δtet(j1, j2, j3, j4, j5, j6)
    #     return GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[])
    # end
    #admissibility checks in the main file

    Tc2 = qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    series = q6jseries_symb(j1, j2, j3, j4, j5, j6)
    return GenericResult(Tc2, series)
end

# ============================================================
# Quantum 3j symbols (Optimized)
# ============================================================

function q3jsummand_symb_first(z::Int, α1::Int, α2::Int, β1::Int, β2::Int, β3::Int)
    sign_val = isodd(Int(z + α1 - α2)) ? -1 : 1
    res = CycloMonomial(sign_val, 0, Int[])
    res = mul_qfact(res, z, -1)
    res = mul_qfact(res, α1 + z, -1); res = mul_qfact(res, α2 + z, -1)
    res = mul_qfact(res, β1 - z, -1); res = mul_qfact(res, β2 - z, -1); res = mul_qfact(res, β3 - z, -1)
    return res
end

function q3jseries_symb(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin)
    α1 = Int(j3 - j2 + m1) 
    α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 - j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    S_z = CycloMonomial[]
    z_min = max(α1, α2, 0)
    z_max = min(β1, β2, β3)

    if z_min <= z_max
        term = q3jsummand_symb_first(z_min, α1, α2, β1, β2, β3)
        push!(S_z, term)

        for z in z_min+1 : z_max
            next_term = CycloMonomial(-term.sign, term.z_pow, copy(term.exps))
            
            next_term = mul_qint(next_term, β1 - z + 1, 1)
            next_term = mul_qint(next_term, β2 - z + 1, 1)
            next_term = mul_qint(next_term, β3 - z + 1, 1)
            
            next_term = mul_qint(next_term, z, -1)
            next_term = mul_qint(next_term, α1 + z, -1)
            next_term = mul_qint(next_term, α2 + z, -1)
            
            push!(S_z, next_term)
            term = next_term
        end
    end
    return S_z
end

function qracah3j_generic(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    # if !δ(j1, j2, j3) || !iszero(m1 + m2 + m3)
    #     return GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[])
    # end
    #admissibility checks in the main file!
    
    pref_sq = qdelta2_symb(j1, j2, j3) 
    pref_sq = mul_qfact(pref_sq, Int(j1+m1), 1); pref_sq = mul_qfact(pref_sq, Int(j1-m1), 1)
    pref_sq = mul_qfact(pref_sq, Int(j2+m2), 1); pref_sq = mul_qfact(pref_sq, Int(j2-m2), 1)
    pref_sq = mul_qfact(pref_sq, Int(j3-m1-m2), 1); pref_sq = mul_qfact(pref_sq, Int(j3+m1+m2), 1)
              
    series = q3jseries_symb(j1, j2, j3, m1, m2)
    return GenericResult(pref_sq, series)
end

# ============================================================
# Classical Evaluation (q = 1)
# ============================================================

function phi_at_one(d::Int)
    d <= 1 && return 1
    p = 2
    while d % p != 0
        p += 1
        p * p > d && (p = d; break)
    end
    temp = d
    while temp % p == 0
        temp ÷= p
    end
    return temp == 1 ? p : 1
end

function evaluate_classical(M::CycloMonomial)
    M.sign == 0 && return big(0) // 1
    num = BigInt(1)
    den = BigInt(1)
    
    @inbounds for d in 2:length(M.exps)
        e = M.exps[d]
        e == 0 && continue
        val = phi_at_one(d)
        val == 1 && continue 
        
        if e > 0
            num *= BigInt(val)^e
        elseif e < 0
            den *= BigInt(val)^abs(e)
        end
    end
    
    final_num = M.sign == -1 ? -num : num
    return final_num // den
end

function qracah6j_classical(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    q6j = qracah6j_generic(j1, j2, j3, j4, j5, j6)
    q6j.pref_sq.sign == 0 && return 0.0
    
    sumz = sum(evaluate_classical.(q6j.series))
    pref_sq_val = evaluate_classical(q6j.pref_sq)
    
    return Float64(sqrt(BigFloat(pref_sq_val)) * BigFloat(sumz))
end

# Inside src/Symbolics.jl

function qracah3j_classical(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin = -m1-m2)
    q3j = qracah3j_generic(j1, j2, j3, m1, m2, m3)
    
    if q3j.pref_sq.sign == 0
        return 0.0
    end
    
    sumz = sum(evaluate_classical.(q3j.series))
    pref_sq_val = evaluate_classical(q3j.pref_sq)
    
    return Float64(sqrt(BigFloat(pref_sq_val)) * BigFloat(sumz))
end