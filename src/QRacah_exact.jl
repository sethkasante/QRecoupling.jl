module ExactSymbolicQRacah

using Nemo

export CycloMonomial, ExactSU2kModel, exact_qracah6j

# ============================================================
# 1. Symbolic Data Structures
# ============================================================

"""
    CycloMonomial

Represents exactly a product of cyclotomic polynomials evaluated at q:
Value = sign * q^(q_pow) * Π (Φ_d(q))^(exps[d])

This allows exact multiplication, division, and square roots via simple integer arithmetic.
"""
struct CycloMonomial
    sign::Int
    q_pow::Int
    exps::Dict{Int, Int}
end

# Exact Symbolic Multiplication
function Base.:*(a::CycloMonomial, b::CycloMonomial)
    exps = copy(a.exps)
    for (d, e) in b.exps
        val = get(exps, d, 0) + e
        if val == 0
            delete!(exps, d)
        else
            exps[d] = val
        end
    end
    CycloMonomial(a.sign * b.sign, a.q_pow + b.q_pow, exps)
end

# Exact Symbolic Division
function Base.:/(a::CycloMonomial, b::CycloMonomial)
    exps = copy(a.exps)
    for (d, e) in b.exps
        val = get(exps, d, 0) - e
        if val == 0
            delete!(exps, d)
        else
            exps[d] = val
        end
    end
    CycloMonomial(a.sign * b.sign, a.q_pow - b.q_pow, exps)
end

"""
    exact_sqrt(M::CycloMonomial) -> (CycloMonomial, CycloMonomial)

Splits a CycloMonomial M into M_out and M_in such that M = (M_out^2) * M_in.
M_out is pulled outside the square root, M_in remains inside.
"""
function exact_sqrt(M::CycloMonomial)
    if M.sign < 0
        error("Cannot take exact real sqrt of negative CycloMonomial")
    end
    
    out_exps = Dict{Int, Int}()
    in_exps  = Dict{Int, Int}()
    
    for (d, e) in M.exps
        out_e = div(e, 2)
        in_e  = e % 2
        if out_e != 0; out_exps[d] = out_e; end
        if in_e != 0;  in_exps[d]  = in_e;  end
    end
    
    out_q = div(M.q_pow, 2)
    in_q  = M.q_pow % 2
    
    M_out = CycloMonomial(1, out_q, out_exps)
    M_in  = CycloMonomial(1, in_q, in_exps)
    
    return M_out, M_in
end

# ============================================================
# 2. The Model and Precomputation
# ============================================================

struct ExactSU2kModel
    k::Int
    K::AnticNumberField            # The Q(q) Cyclotomic Number Field
    z::nf_elem                     # The primitive 2N-th root of unity (q)
    Phi_cache::Dict{Int, nf_elem}  # Pre-evaluated Φ_d(q) to avoid recomputing
    qfacts::Vector{CycloMonomial}  # Pre-factorized symbolic quantum factorials
end

function ExactSU2kModel(k::Int)
    N = k + 2
    K, z = cyclotomic_field(2N, "z")
    
    # 1. Precompute Cyclotomic Polynomials evaluated at z
    Zx, x = polynomial_ring(ZZ, "x")
    max_d = 4 * (k + 2) # Upper bound for divisor indices
    Phi_cache = Dict{Int, nf_elem}()
    for d in 3:max_d
        poly = cyclotomic(d, x)
        Phi_cache[d] = evaluate(poly, z) 
    end
    
    # 2. Precompute Symbolic Quantum Factorials [n]_q!
    # identity: [n]_q = q^{1-n} * Π_{d|2n, d>2} Φ_d(q)
    max_n = 2 * k + 4 
    qfacts = Vector{CycloMonomial}(undef, max_n + 1)
    qfacts[1] = CycloMonomial(1, 0, Dict{Int,Int}()) # [0]! = 1
    
    for n in 1:max_n
        exps = Dict{Int,Int}()
        for d in 3:2n
            if 2n % d == 0
                exps[d] = 1
            end
        end
        qn = CycloMonomial(1, 1 - n, exps)
        qfacts[n+1] = qfacts[n] * qn
    end
    
    return ExactSU2kModel(k, K, z, Phi_cache, qfacts)
end

# ============================================================
# 3. Field Evaluation & BigFloat Bridge
# ============================================================

"""
    evaluate_nf(M::CycloMonomial, model::ExactSU2kModel) -> nf_elem

Takes a symbolic dictionary of quantum primes and evaluates it exactly
into the Nemo cyclotomic field.
"""
function evaluate_nf(M::CycloMonomial, model::ExactSU2kModel)
    K = model.K
    res = K(M.sign)
    
    if M.q_pow > 0
        res *= model.z^M.q_pow
    elseif M.q_pow < 0
        res *= inv(model.z)^(-M.q_pow)
    end
    
    for (d, e) in M.exps
        if e > 0
            res *= model.Phi_cache[d]^e
        elseif e < 0
            res *= inv(model.Phi_cache[d])^(-e)
        end
    end
    return res
end

"""
    to_bigfloat(elem::nf_elem, k::Int) -> BigFloat

Safely projects the exact algebraic field element back to a real BigFloat.
"""
function to_bigfloat(elem::nf_elem, k::Int)
    N = k + 2
    z_bf = exp(im * big(pi) / N) 
    
    val = zero(Complex{BigFloat})
    poly = Nemo.polynomial(elem)
    
    for i in 0:degree(poly)
        c = coeff(poly, i)
        num = BigFloat(Nemo.numerator(c))
        den = BigFloat(Nemo.denominator(c))
        val += (num / den) * (z_bf^i)
    end
    
    return real(val) 
end

# (Admissibility conditions omitted for brevity, assume qδtet is standard)
@inline ishalfInt(j)::Bool = isinteger(2*j) && j ≥ 0
@inline δ(j1, j2, j3)::Bool = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)
@inline qδ(j1, j2, j3, k)::Bool = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 
@inline qδtet(j1, j2, j3, j4, j5, j6, k)::Bool = 
    ishalfInt(j1) && ishalfInt(j2) && ishalfInt(j3) && ishalfInt(j4) && 
    ishalfInt(j5) && ishalfInt(j6) && 
    qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)

# ============================================================
# 4. The Final 6j Evaluator
# ============================================================

function exact_qracah6j(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return big(0.0)
    end
    
    tab = model.qfacts
    
    # 1. Compute symbolic Delta^2 prefactor
    function delta2(a, b, c)
        x, y, z, w = Int(a+b-c), Int(a-b+c), Int(-a+b+c), Int(a+b+c)
        return (tab[x+1] * tab[y+1] * tab[z+1]) / tab[w+2]
    end
    
    Pref_sq = delta2(j1,j2,j3) * delta2(j1,j5,j6) * delta2(j2,j4,j6) * delta2(j3,j4,j5)
    
    # Extract exact symbolic square root
    Pref_out, Pref_in = exact_sqrt(Pref_sq)
    
    # 2. Run the exact Racah Sum
    α1 = Int(j1 + j2 + j3) 
    α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6) 
    α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5)
    β2 = Int(j1 + j3 + j4 + j6) 
    β3 = Int(j2 + j3 + j5 + j6)

    zrange = max(α1, α2, α3, α4):min(β1, β2, β3)
    
    S_sum = model.K(0) # Initialize sum exactly in the number field
    
    for z in zrange
        # Pure symbolic construction of the term
        num = tab[z+2]
        den = tab[z-α1+1] * tab[z-α2+1] * tab[z-α3+1] * tab[z-α4+1] *
              tab[β1-z+1] * tab[β2-z+1] * tab[β3-z+1]
              
        term_symb = num / den
        
        # Apply parity sign
        sign_z = iseven(z) ? term_symb.sign : -term_symb.sign
        term_symb = CycloMonomial(sign_z, term_symb.q_pow, term_symb.exps)
        
        # Evaluate to number field and add (no expression swell here!)
        S_sum += evaluate_nf(term_symb, model)
    end
    
    # 3. Final Assembly
    # Evaluate exact components back to BigFloat
    out_bf = to_bigfloat(evaluate_nf(Pref_out, model), model.k)
    in_bf  = to_bigfloat(evaluate_nf(Pref_in, model), model.k)
    S_bf   = to_bigfloat(S_sum, model.k)
    
    return out_bf * sqrt(abs(in_bf)) * S_bf
end

end # module