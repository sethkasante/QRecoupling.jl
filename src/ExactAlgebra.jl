# src/ExactAlgebra.jl

struct ExactSU2kModel
    k::Int
    K::AnticNumberField            
    z::nf_elem                     
    q_facts::Vector{nf_elem}       
end

# ============================================================
# Exact Algebraic Model using (Nemo.jl)
# ============================================================

function ExactSU2kModel(k::Int)
    N = k + 2
    K, z = cyclotomic_field(2N, "ζ") 
    
    q_facts = Vector{nf_elem}(undef, k + 3)
    q_facts[1] = K(1) 
    
    z_inv = inv(z)
    den_inv = inv(z - z_inv)
    z_n = z
    z_inv_n = z_inv
    
    @inbounds for n in 1:(k+2)
        q_int = (z_n - z_inv_n) * den_inv
        q_facts[n+1] = q_facts[n] * q_int
        z_n *= z
        z_inv_n *= z_inv
    end
    return ExactSU2kModel(k, K, z, q_facts)
end

# ============================================================
# Exact Cyclo Field Evaluations 
# ============================================================

@inline function q_delta2_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin)
    a = Int(j1 + j2 - j3); b = Int(j1 - j2 + j3)
    c = Int(-j1 + j2 + j3); d = Int(j1 + j2 + j3)
    
    num = model.q_facts[a+1] * model.q_facts[b+1] * model.q_facts[c+1]
    den = model.q_facts[d+2]
    return num * inv(den)
end

@inline function qtricoeff2_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    return q_delta2_exact(model, j1, j2, j3) * q_delta2_exact(model, j1, j5, j6) * q_delta2_exact(model, j2, j4, j6) * q_delta2_exact(model, j3, j4, j5)
end

function q6jseries_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3, model.k) 
    sum_cf = model.K(0)

    @inbounds for z in zrange
        num = model.q_facts[z+2]
        den = model.q_facts[z-α1+1] * model.q_facts[z-α2+1] * model.q_facts[z-α3+1] * model.q_facts[z-α4+1] * model.q_facts[β1-z+1] * model.q_facts[β2-z+1] * model.q_facts[β3-z+1]
        
        term = num * inv(den)
        sum_cf += iseven(z) ? term : -term
    end
    return sum_cf
end


function qracah6j_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    Tc2 = qtricoeff2_exact(model, j1, j2, j3, j4, j5, j6)
    Sum_cf = q6jseries_exact(model, j1, j2, j3, j4, j5, j6)
    return ExactResult(model.k, Tc2, Sum_cf)
end


# ============================================================
# Exact Quantum 3j Symbol
# ============================================================

@inline function q3j_pref_sq_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin)
    delta2 = q_delta2_exact(model, j1, j2, j3)
    
    # Exact memory lookups for the 6 factorials
    facts = model.q_facts[Int(j1+m1)+1] * model.q_facts[Int(j1-m1)+1] * model.q_facts[Int(j2+m2)+1] * model.q_facts[Int(j2-m2)+1] * model.q_facts[Int(j3-m1-m2)+1] * model.q_facts[Int(j3+m1+m2)+1]
            
    return delta2 * facts
end

function q3jseries_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin)
    α1 = Int(j3 - j2 + m1) 
    α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 - j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    # Corrected z_min bound to prevent negative factorials
    zrange = max(-α1, -α2, 0):min(β1, β2, β3, model.k)
    sum_cf = model.K(0)
    
    # Pre-calculate the static part of the phase exponent
    phase_offset = α1 - α2

    @inbounds for z in zrange
        den = model.q_facts[z+1] * model.q_facts[α1+z+1] * model.q_facts[α2+z+1] * model.q_facts[β1-z+1] * model.q_facts[β2-z+1] * model.q_facts[β3-z+1]
        
        term = inv(den)
        # Phase is (-1)^{z + α1 - α2}
        sum_cf += isodd(z + phase_offset) ? -term : term
    end
    
    return sum_cf
end

function qracah3j_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin)
    # The m-admissibility validated in the main file
    pref_sq = q3j_pref_sq_exact(model, j1, j2, j3, m1, m2)
    Sum_cf = q3jseries_exact(model, j1, j2, j3, m1, m2)
    
    return ExactResult(model.k, pref_sq, Sum_cf)
end



# ============================================================
# Float Projection (Horner's Method)
# ============================================================

function evaluate_exact(res::ExactResult, prec::Int=256)
    if res.pref_sq == 0 && res.sum_cf == 0
        return zero(BigFloat)
    end
    
    setprecision(BigFloat, prec) do
        N = res.k + 2
        target_z = cispi(big"1.0" / N)
        
        sum_bf = horner_eval(res.sum_cf, target_z)
        pref2_bf = horner_eval(res.pref_sq, target_z)
        
        return real(sqrt(pref2_bf) * sum_bf)
    end
end

function horner_eval(ev::nf_elem, root_val::Complex{BigFloat})
    d = degree(parent(ev)) - 1
    d < 0 && return zero(Complex{BigFloat})
    
    res = Complex{BigFloat}(coeff(ev, d))
    for i in (d-1):-1:0
        c = coeff(ev, i)
        val = BigFloat(Nemo.numerator(c)) / BigFloat(Nemo.denominator(c))
        res = res * root_val + val
    end
    return res
end