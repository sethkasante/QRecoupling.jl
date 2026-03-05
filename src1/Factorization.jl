# src/Factorization.jl

"""
    qfactorial_symb(n::Int) -> CycloMonomial

Computes the exact symbolic representation of [n]_q! using the identity:
[n]_q! = q^{-n(n-1)/4} * Π_{d=2}^n Φ_d(q)^floor(n/d)
"""
function qfactorial_symb(n::Int)
    if n == 0
        return CycloMonomial(1, 0//1, Dict{Int, Int}())
    end

    exps = Dict{Int, Int}()
    for d in 2:n
        power = div(n, d) # Integer division matches floor(n/d)
        if power > 0
            exps[d] = power
        end
    end
    
    # Calculate the fractional q power
    q_pow = - (n * (n - 1)) // 4
    
    return CycloMonomial(1, q_pow, exps)
end

"""
    exact_sqrt(M::CycloMonomial) -> (CycloMonomial, CycloMonomial)

Splits M into (M_out^2) * M_in for exact square root extraction.
"""
function exact_sqrt(M::CycloMonomial)
    M.sign < 0 && error("Cannot take exact real sqrt of negative CycloMonomial")
    
    out_exps = Dict{Int, Int}()
    in_exps  = Dict{Int, Int}()
    
    for (d, e) in M.exps
        out_e = div(e, 2)
        in_e  = e % 2
        if out_e != 0; out_exps[d] = out_e; end
        if in_e != 0;  in_exps[d]  = in_e;  end
    end
    
    # Rational division for the q_pow
    out_q = M.q_pow / 2
    in_q  = 0//1 # We absorb the remainder entirely into the out_q if we define v = q^(1/2) 
                 # Wait, if q_pow is a fraction, we must be careful.
                 # Let's keep it rigorous:
    
    # Actually, the simplest approach is just halving the rational exponent directly.
    # Since we evaluate at q eventually, fractional powers of q are just powers of v.
    M_out = CycloMonomial(1, M.q_pow / 2, out_exps)
    M_in  = CycloMonomial(1, 0//1, in_exps) # Only the cyclotomic polynomials leave a remainder
    
    return M_out, M_in
end