#TQFT.jl

"""
    qdim(j, k) -> CycloMonomial
    qdim(j, model::ExactSU2kModel) -> nf_elem

Returns the quantum dimension of a spin j representation: [2j+1]_q.
[cite: 11, 13, 26, 27]
"""
function qdim_symb(j)
    n = Int(2j + 1)
    # [n]_q = [n]_q! / [n-1]_q!
    return qfactorial_symb(n) / qfactorial_symb(n-1) [cite: 11, 13]
end

function qdim(j, model::ExactSU2kModel)
    return evaluate_cyclof(qdim_symb(j), model) [cite: 2, 3]
end

"""
    rmatrix(j1, j2, j3) -> CycloMonomial

Returns the braiding R-symbol (eigenvalue) for the process j1 ⊗ j2 -> j3.
R^{j1,j2}_{j3} = (-1)^{j1+j2-j3} q^{(j3(j3+1) - j1(j1+1) - j2(j2+1))/2}
In our z-basis (z = q^{1/2}), the exponent is j3(j3+1) - j1(j1+1) - j2(j2+1).
[cite: 13, 18, 26]
"""
function rmatrix_symb(j1, j2, j3)
    if !δ(j1, j2, j3) return CycloMonomial(0, 0, Dict{Int,Int}()) end
    
    # Phase sign (-1)^{j1+j2-j3}
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1 #[cite: 18]
    
    # Casimir power: j(j+1). We use 2j to keep everything as integers.
    # Exponent = [2j3(2j3+2) - 2j1(2j1+2) - 2j2(2j2+2)] / 4
    # Since z = q^{1/2}, the z_pow is exactly this numerator / 2 if we use q^{1/2} logic.
    function casimir_pow(j)
        return Int(2j * (2j + 2)) ÷ 4
    end
    
    z_pow = casimir_pow(j3) - casimir_pow(j1) - casimir_pow(j2) [cite: 12, 13]
    
    return CycloMonomial(s, z_pow, Dict{Int,Int}()) [cite: 18]
end

"""
    fsymbol(j1, j2, j3, j4, j5, j6, model::ExactSU2kModel)

Returns the F-matrix element (fusion/recoupling coefficient).
Often normalized as: F = (-1)^{j1+j2+j3+j4} √([2j3+1][2j6+1]) {j1 j2 j3; j4 j5 j6}_q
[cite: 5, 15, 16, 17, 26]
"""
function fsymbol_exact(j1, j2, j3, j4, j5, j6, model::ExactSU2kModel)
    # Admissibility check [cite: 5]
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) return model.K(0) end
    
    # Get exact 6j components [cite: 5, 6]
    pref2_nf, sum_nf = qracah6j_symb(model, j1, j2, j3, j4, j5, j6)
    
    # Fusion normalization dimensions [cite: 11]
    d3_sq = evaluate_cyclof(qdim_symb(j3)^2, model)
    d6_sq = evaluate_cyclof(qdim_symb(j6)^2, model)
    
    # Return as a tuple (Sign, Square_of_Prefactor, Racah_Sum) 
    # to allow the user to handle the square root numerically or symbolically.
    sign = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    return (sign, pref2_nf * d3_sq * d6_sq, sum_nf)
end


