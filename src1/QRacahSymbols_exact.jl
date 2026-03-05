module ExactQRacah

using Nemo

export SU2kExactModel, exact_qracah6j, exact_racah_sum

# ============================================================
# Exact Algebraic Model
# ============================================================

struct SU2kExactModel
    k::Int
    K::AnticNumberField  # The Cyclotomic Field
    z::nf_elem           # The primitive root of unity 'q'
    qnfact::Vector{nf_elem} # Exact quantum factorials
end

function SU2kExactModel(k::Int)
    N = k + 2
    # Define the cyclotomic field Q(zeta_2N)
    # Nemo will automatically reduce all math modulo the cyclotomic polynomial!
    K, z = cyclotomic_field(2 * N, "z")
    
    # Calculate exact quantum integers: [n]_q = (z^n - z^-n) / (z - z^-1)
    den = z - inv(z)
    
    qnfact = Vector{elem_type(K)}(undef, N + 1)
    qnfact[1] = K(1) # [0]! = 1
    
    curr_fact = K(1)
    for n in 1:N
        if n == 1
            qn = K(1)
        else
            qn = (z^n - inv(z)^n) // den
        end
        curr_fact *= qn
        qnfact[n+1] = curr_fact
    end
    
    return SU2kExactModel(k, K, z, qnfact)
end

# (Assume your standard admissibility functions qδtet, ishalfInt etc. are included here)
@inline ishalfInt(j)::Bool = isinteger(2*j) && j ≥ 0
@inline δ(j1, j2, j3)::Bool = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)
@inline qδ(j1, j2, j3, k)::Bool = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 
@inline qδtet(j1, j2, j3, j4, j5, j6, k)::Bool = 
    ishalfInt(j1) && ishalfInt(j2) && ishalfInt(j3) && ishalfInt(j4) && 
    ishalfInt(j5) && ishalfInt(j6) && 
    qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && 
    qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)

# ============================================================
# Exact Racah Sum Evaluation
# ============================================================

"""
Computes the alternating Racah sum EXACTLY using cyclotomic field arithmetic.
This completely prevents catastrophic cancellation and expression swell.
"""
function exact_racah_sum(model::SU2kExactModel, j1, j2, j3, j4, j5, j6)
    tab = model.qnfact
    
    α1 = Int(j1 + j2 + j3) 
    α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6) 
    α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5)
    β2 = Int(j1 + j3 + j4 + j6) 
    β3 = Int(j2 + j3 + j5 + j6)

    zrange = max(α1, α2, α3, α4):min(β1, β2, β3)
    
    # Initialize sum as 0 inside the Number Field
    S = model.K(0) 
    
    for z in zrange
        num = tab[z+2]
        den = tab[z-α1+1] * tab[z-α2+1] * tab[z-α3+1] * tab[z-α4+1] *
              tab[β1-z+1] * tab[β2-z+1] * tab[β3-z+1]
        
        # Exact algebraic division
        term = num // den
        
        if iseven(z)
            S += term
        else
            S -= term
        end
    end
    return S
end

# ============================================================
# Bridge back to BigFloat
# ============================================================

"""
Evaluates the exact number field element to a BigFloat.
"""
function evaluate_nf(elem::nf_elem, k::Int)
    N = k + 2
    z_bf = exp(im * big(pi) / N) 
    
    val = zero(Complex{BigFloat})
    poly = Nemo.polynomial(elem)
    
    # Evaluate the bounded polynomial using exact rational coefficients
    for i in 0:degree(poly)
        c = coeff(poly, i)
        num = BigFloat(Nemo.numerator(c))
        den = BigFloat(Nemo.denominator(c))
        val += (num / den) * (z_bf^i)
    end
    
    # The Racah sum is purely real, drop the 0im
    return real(val) 
end

"""
Computes the purely numerical Δ triangle coefficients (safe from cancellation)
"""
function qdim_bf(j, k::Int)
    θ = big(pi) / (k + 2)
    return sin(Int(2j + 1) * θ) / sin(θ)
end

function qnfact_bf(n::Int, k::Int)
    θ = big(pi) / (k + 2)
    res = big(1.0)
    for i in 1:n
        res *= sin(i * θ) / sin(θ)
    end
    return res
end

function delta_bf(j1, j2, j3, k::Int)
    a, b, c, d = Int(j1+j2-j3), Int(j1-j2+j3), Int(-j1+j2+j3), Int(j1+j2+j3)
    num = qnfact_bf(a, k) * qnfact_bf(b, k) * qnfact_bf(c, k)
    den = qnfact_bf(d + 1, k)
    return sqrt(num / den)
end

"""
The final exact-to-float hybrid 6j evaluator.
"""
function exact_qracah6j(model::SU2kExactModel, j1, j2, j3, j4, j5, j6)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return big(0.0)
    end
    
    # 1. Compute the unstable sum EXACTLY in the cyclotomic field
    exact_S = exact_racah_sum(model, j1, j2, j3, j4, j5, j6)
    
    # 2. Evaluate to a standard BigFloat
    S_bf = evaluate_nf(exact_S, model.k)
    
    # 3. Compute the perfectly stable square roots numerically
    T = delta_bf(j1, j2, j3, model.k) * delta_bf(j1, j5, j6, model.k) *
        delta_bf(j2, j4, j6, model.k) * delta_bf(j3, j4, j5, model.k)
        
    return T * S_bf
end

end # end module