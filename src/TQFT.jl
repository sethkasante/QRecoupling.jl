# src/TQFT.jl

# ============================================================
# 1. Quantum Dimensions [2j+1]_q
# ============================================================

"""
    qdim_symb(j) -> CycloMonomial

Returns the exact symbolic quantum dimension.
"""
function qdim_symb(j)
    n = Int(2j + 1)
    return qfactorial_symb(n) / qfactorial_symb(n - 1)
end

"""
    qdim_exact(model::ExactSU2kModel, j) -> nf_elem

Returns the exact quantum dimension in the cyclotomic number field.
"""
function qdim_exact(model::ExactSU2kModel, j)
    # Map the symbolic representation directly into the field
    return evaluate_cyclofield(qdim_symb(j), model)
end

# ============================================================
# 2. R-Matrix (Braiding)
# ============================================================

"""
    rmatrix_symb(j1, j2, j3) -> CycloMonomial

Returns the exact symbolic R-matrix eigenvalue.
"""
function rmatrix_symb(j1, j2, j3)
    if !δ(j1, j2, j3)
        return CycloMonomial(0, 0, Int[])
    end
    
    phase_exp = Int(j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1
    
    # Int[] means there are no cyclotomic polynomial factors (Φ_d), just the phase and sign
    return CycloMonomial(s, phase_exp, Int[])
end

"""
    rmatrix_exact(model::ExactSU2kModel, j1, j2, j3) -> nf_elem

Returns the exact R-matrix eigenvalue in the cyclotomic number field.
"""
function rmatrix_exact(model::ExactSU2kModel, j1, j2, j3)
    return evaluate_cyclofield(rmatrix_symb(j1, j2, j3), model)
end

# ============================================================
# 3. F-Symbol (Fusion / Normalized 6j)
# ============================================================
# Formula: F = (-1)^{j1+j2+j3+j4} * √([2j3+1][2j6+1]) * {6j}

"""
    fsymbol_generic(j1, j2, j3, j4, j5, j6) -> GenericResult

Returns the fully symbolic F-symbol. 
The dimension prefactors are absorbed into the Triangle squared prefactor, 
and the global sign phase is distributed into the Racah series.
"""
function fsymbol_generic(j1, j2, j3, j4, j5, j6)
    res_6j = qracah6j_generic(j1, j2, j3, j4, j5, j6)
    
    if res_6j.pref_sq.sign == 0
        return res_6j
    end
    
    # 1. Scale the squared prefactor by [2j3+1] * [2j6+1]
    d3_symb = qdim_symb(j3)
    d6_symb = qdim_symb(j6)
    new_pref_sq = res_6j.pref_sq * d3_symb * d6_symb
    
    # 2. Distribute the global phase into the summand series
    phase_s = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    
    if phase_s == -1
        # Flip the sign of every monomial in the series
        new_series = [CycloMonomial(-term.sign, term.z_pow, term.exps) for term in res_6j.series]
    else
        new_series = res_6j.series
    end
    
    return GenericResult(new_pref_sq, new_series)
end

"""
    fsymbol_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6) -> ExactResult

Returns the exact F-symbol evaluated in the cyclotomic number field.
"""
function fsymbol_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    res_6j = qracah6j_exact(model, j1, j2, j3, j4, j5, j6)
    
    # Check if the 6j symbol evaluated to exactly 0 (admissibility failure)
    if res_6j.pref_sq == 0
        return res_6j
    end
    
    # 1. Scale the squared field prefactor by exact dimensions
    d3_exact = qdim_exact(model, j3)
    d6_exact = qdim_exact(model, j6)
    new_pref_sq = res_6j.pref_sq * d3_exact * d6_exact
    
    # 2. Apply the global phase to the evaluated field sum
    phase_s = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    new_sum = phase_s == 1 ? res_6j.sum_cf : -res_6j.sum_cf
    
    return ExactResult(model.k, new_pref_sq, new_sum)
end

# ============================================================
# 4. G-Symbol (Tetrahedral Weight)
# ============================================================
# Formula: G = √(Π [2j_i + 1]) * {6j}

function gsymbol_generic(j1, j2, j3, j4, j5, j6)
    res_6j = qracah6j_generic(j1, j2, j3, j4, j5, j6)
    if res_6j.pref_sq.sign == 0 
        return res_6j 
    end
    
    # Multiply the squared prefactor by the product of all 6 edge dimensions
    dim_prod = qdim_symb(j1) * qdim_symb(j2) * qdim_symb(j3) * qdim_symb(j4) * qdim_symb(j5) * qdim_symb(j6)
               
    new_pref_sq = res_6j.pref_sq * dim_prod
    return GenericResult(new_pref_sq, res_6j.series)
end

function gsymbol_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    res_6j = qracah6j_exact(model, j1, j2, j3, j4, j5, j6)
    if res_6j.pref_sq == 0 return res_6j end
    
    dim_prod = qdim_exact(model, j1) * qdim_exact(model, j2) * qdim_exact(model, j3) * qdim_exact(model, j4) * qdim_exact(model, j5) * qdim_exact(model, j6)
               
    new_pref_sq = res_6j.pref_sq * dim_prod
    return ExactResult(model.k, new_pref_sq, res_6j.sum_cf)
end

# ============================================================
# Extension Suite: TQFT Functions (Numeric)
# ============================================================

"""
    qdim_numeric(j, model::NumericSU2kModel)

Quantum dimension [2j+1]_q using the log-tables.
"""
function qdim_numeric(j, model::NumericSU2kModel)
    n = Int(2j + 1)
    # [n]_q = [n]_q! / [n-1]_q!
    # exp(log([n]!) - log([n-1]!))
    return exp(model.logqnfact[n+1] - model.logqnfact[n])
end

"""
    fsymbol_numeric(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)

Returns the normalized F-matrix element.
F = (-1)^{j1+j2+j3+j4} √([2j3+1][2j6+1]) {j1 j2 j3; j4 j5 j6}_q
"""
function fsymbol_numeric(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)
    val_6j = _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
    
    d3 = qdim_numeric(j3, model)
    d6 = qdim_numeric(j6, model)
    
    phase = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    return phase * sqrt(d3 * d6) * val_6j
end

"""
    rmatrix_numeric(j1, j2, j3, k::Int)

Returns the braiding R-symbol. Since this is purely a phase, it is evaluated 
directly without the log-tables.
"""
function rmatrix_numeric(j1, j2, j3, k::Int)
    phase_exp = Int(j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1
    return s * cispi(phase_exp / (k + 2))
end

"""
    gsymbol_numeric(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)

The G-symbol (Tetrahedral Weight) used in state sums.
G = {6j} * √(Π [2j_i + 1]_q)
"""
function gsymbol_numeric(model::NumericSU2kModel, j1, j2, j3, j4, j5, j6)
    val_6j = _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
    
    dims_prod = qdim_numeric(j1, model) * qdim_numeric(j2, model) * qdim_numeric(j3, model) * qdim_numeric(j4, model) * qdim_numeric(j5, model) * qdim_numeric(j6, model)
                
    return val_6j * sqrt(dims_prod)
end