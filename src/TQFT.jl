# ============================================================
# Quantum Integers [n] and Dimensions [2j+1]_q
# ============================================================

"""
    qint(n::Int; mode=:generic)
    qint(n::Int, k::Int; mode=:exact, T=Float64)
    
Public API for quantum integers.
"""
function qint(n::Int; mode=:generic)
    @assert mode == :generic "Numeric and exact modes require level k. Use qint(n, k; mode=...)"
    n == 0 && return CycloMonomial(0, 0, Int[]) 
    n == 1 && return CycloMonomial(1, 0, Int[]) 
    
    # [n]_q factors into Φ_d(q) for all divisors d of n (d > 1)
    exps = zeros(Int, n)
    for d in 2:n
        if n % d == 0
            exps[d] = 1
        end
    end
    return CycloMonomial(1, 0, exps)
end

function qint(n::Int, k::Int; mode=:numeric, T::Type{<:AbstractFloat}=Float64)
    if mode == :exact
        n == 0 && return Nemo.QQFieldElem(0)
        n == 1 && return Nemo.QQFieldElem(1)
        
        # Performance: Only use the full model cache if it's already there.
        # Otherwise, compute a single element to avoid building the whole factorial table.
        if haskey(EXACT_MODEL_CACHE, k)
            model = EXACT_MODEL_CACHE[k]
            return model.q_facts[n+1] * inv(model.q_facts[n])
        else
            K, z = Nemo.cyclotomic_field(2 * (k + 2), "ζ")
            z_inv = inv(z)
            return (z^n - z_inv^n) * inv(z - z_inv)
        end
        
    elseif mode == :numeric
        n == 0 && return zero(T)
        n == 1 && return one(T)
        θ = T(π) / T(k + 2)
        return sin(T(n) * θ) / sin(θ)
    else
        throw(ArgumentError("Unknown mode: $mode"))
    end
end


#-------- Internal functions ------- 

function qdim_symb(j::Spin)
    n = Int(2j + 1)
    #optimized for dim computations 
    buf = SymbolicBuffer(n + 1)
    add_qint!(buf, n, 1)
    return snapshot(buf)
end

function qdim_exact(model::ExactSU2kModel, j::Spin)
    n = Int(2j + 1)
    # [n] = [n]! / [n-1]!
    return model.q_facts[n+1] * inv(model.q_facts[n])
end

function qdim_numeric(j::Spin, model::NumericSU2kModel{T})::T where {T}
    n = Int(2j + 1)
    return exp(model.logqnfact[n+1] - model.logqnfact[n])
end

#--------- Public API for quantum dimensions  ----------------

function qdim(j::Spin; mode=:generic)
    if !ishalfInt(j)
        if mode == :generic return CycloMonomial(0, 0, Int[]) end
        if mode == :classical return 0.0 end
    end
    
    mode == :generic && return qdim_symb(j) 
    error("Mode :$mode requires a level `k`.")
end

# mainly for exact and numeric 
function qdim(j::Spin, k::Int; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    if mode == :generic || mode == :classical return qdim(j; mode=mode) end
    
    # j must be a valid spin, and 2j <= k
    if !ishalfInt(j) || (2j > k)
        if mode == :exact
            # check cache first to avoid instantiation
            if haskey(EXACT_MODEL_CACHE, k)
                return EXACT_MODEL_CACHE[k].K(0)
            else
                K, _ = Nemo.cyclotomic_field(2*(k+2), "ζ")
                return K(0)
            end
        else
            return zero(T)
        end
    end

    if mode == :exact
        #TODO:Change to haskey otherwise compute qint 
        model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
        return qdim_exact(model, j)
    elseif mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return qdim_numeric(j, model)
    end
    error("Unknown mode: $mode")
end


# ============================================================
# 2. R-Matrix (Braiding)
# ============================================================

# R^{j3}_{j1,j2} phase calculation
function rmatrix_symb(j1::Spin, j2::Spin, j3::Spin)
    phase_exp = Int(j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1
    return CycloMonomial(s, phase_exp, Int[])
end

function rmatrix_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin)
    phase_exp = Int(j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1
    res = model.z^phase_exp
    return s == 1 ? res : -res
end

function rmatrix_numeric(j1::Spin, j2::Spin, j3::Spin, k::Int; T::Type{<:AbstractFloat}=Float64)
    phase_exp = Int(j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1
    return s * cispi(T(phase_exp) / T(k + 2))
end

# ------  Public API for Rmatrix Braiding  ------
function rmatrix(j1::Spin, j2::Spin, j3::Spin, k::Union{Int, Nothing}=nothing; mode=:generic, T=Float64)
    # gatekeeper
    if k === nothing
        if !δ(j1, j2, j3) 
            return (mode == :generic ? CycloMonomial(0,0,Int[]) : 0.0)
        end
        
        mode == :generic && return rmatrix_symb(j1, j2, j3)
        error("Mode $mode requires level k")
    else
        if !qδ(j1, j2, j3, k) 
            return (mode == :exact ? qint(0, k; mode=:exact) : zero(T))
        end

        if mode == :exact
            model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
            return rmatrix_exact(model, j1, j2, j3)
        elseif mode == :numeric
            return rmatrix_numeric(j1, j2, j3, k; T=T)
        end
    end
end


# ============================================================
# 3. F-Symbol (Fusion / Normalized 6j)
# ============================================================

function fsymbol_generic(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res_6j = qracah6j_generic(j1, j2, j3, j4, j5, j6)
    res_6j.pref_sq.sign == 0 && return res_6j
    
    # 1. Update the prefactor using the SymbolicBuffer to avoid CycloMonomial multiplications
    cap = length(res_6j.pref_sq.exps) + Int(2 * max(j3, j6) + 2)
    buf = SymbolicBuffer(cap)
    
    # Initialize buffer with the existing 6j prefactor
    buf.sign = res_6j.pref_sq.sign
    buf.z_pow = res_6j.pref_sq.z_pow
    copyto!(buf.exps, 1, res_6j.pref_sq.exps, 1, length(res_6j.pref_sq.exps))
    
    # Add dimensions [2j3+1]_q and [2j6+1]_q
    add_qint!(buf, Int(2j3 + 1), 1)
    add_qint!(buf, Int(2j6 + 1), 1)
    
    new_pref_sq = snapshot(buf)
    
    # 2. Apply overall phase shift to the Racah sum series
    phase_s = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    if phase_s == -1
        new_series = [CycloMonomial(-term.sign, term.z_pow, copy(term.exps)) for term in res_6j.series]
        return GenericResult(new_pref_sq, new_series)
    end
    
    return GenericResult(new_pref_sq, res_6j.series)
end

function fsymbol_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res_6j = qracah6j_exact(model, j1, j2, j3, j4, j5, j6)
    res_6j.pref_sq == 0 && return res_6j
    
    new_pref_sq = res_6j.pref_sq * qdim_exact(model, j3) * qdim_exact(model, j6)
    phase_s = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    new_sum = phase_s == 1 ? res_6j.sum_cf : -res_6j.sum_cf
    
    return ExactResult(model.k, new_pref_sq, new_sum)
end

function fsymbol_numeric(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::T where {T}
    val_6j = _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
    phase = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    return phase * sqrt(qdim_numeric(j3, model) * qdim_numeric(j6, model)) * val_6j
end

#------ Public API for F-symbols ----------

function fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::OptInt=nothing; 
                 mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    
    # Admissibility checks
    if k === nothing
        if !δtet(j1, j2, j3, j4, j5, j6) 
            return GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[])
        else
            return fsymbol_generic(j1, j2, j3, j4, j5, j6)
        end
        error("Mode $mode requires level k")
    else
        # level k admissibility conditions
        if !qδtet(j1, j2, j3, j4, j5, j6, k)
            return (mode == :exact ? ExactResult(k, 0, 0) : zero(T))
        end
        
        if mode == :exact
            model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
            return fsymbol_exact(model, j1, j2, j3, j4, j5, j6)
        elseif mode == :numeric
            model = NumericSU2kModel(k; T=T, prec=prec)
            return fsymbol_numeric(model, j1, j2, j3, j4, j5, j6)
        end
        throw(ArgumentError("Unknown mode: $mode"))
    end
end


# ============================================================
# 4. G-Symbol (Tetrahedral Weight)
# ============================================================

function gsymbol_generic(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res_6j = qracah6j_generic(j1, j2, j3, j4, j5, j6)
    res_6j.pref_sq.sign == 0 && return res_6j
    
    cap = length(res_6j.pref_sq.exps) + Int(2 * max(j1, j2, j3, j4, j5, j6) + 2)
    buf = SymbolicBuffer(cap)
    
    buf.sign = res_6j.pref_sq.sign
    buf.z_pow = res_6j.pref_sq.z_pow
    copyto!(buf.exps, 1, res_6j.pref_sq.exps, 1, length(res_6j.pref_sq.exps))
    
    add_qint!(buf, Int(2j1 + 1), 1); add_qint!(buf, Int(2j2 + 1), 1); add_qint!(buf, Int(2j3 + 1), 1)
    add_qint!(buf, Int(2j4 + 1), 1); add_qint!(buf, Int(2j5 + 1), 1); add_qint!(buf, Int(2j6 + 1), 1)
    
    return GenericResult(snapshot(buf), res_6j.series)
end

function gsymbol_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res_6j = qracah6j_exact(model, j1, j2, j3, j4, j5, j6)
    res_6j.pref_sq == 0 && return res_6j
    
    dim_prod = qdim_exact(model, j1) * qdim_exact(model, j2) * qdim_exact(model, j3) * qdim_exact(model, j4) * qdim_exact(model, j5) * qdim_exact(model, j6)
    return ExactResult(model.k, res_6j.pref_sq * dim_prod, res_6j.sum_cf)
end

function gsymbol_numeric(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::T where {T}
    val_6j = _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
    dims_prod = qdim_numeric(j1, model) * qdim_numeric(j2, model) * qdim_numeric(j3, model) * qdim_numeric(j4, model) * qdim_numeric(j5, model) * qdim_numeric(j6, model)
    return val_6j * sqrt(dims_prod)
end


# ----- Public API for G-symbols ----------- 

function gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::OptInt=nothing; 
                 mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    
    # Admissibility checks
    if k === nothing
        if !δtet(j1, j2, j3, j4, j5, j6) 
            return GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[])
        else
            return gsymbol_generic(j1, j2, j3, j4, j5, j6)
        end
        error("Mode $mode requires level k")
    else
        # level k admissibility conditions
        if !qδtet(j1, j2, j3, j4, j5, j6, k)
            return (mode == :exact ? ExactResult(k, 0, 0) : zero(T))
        end
        
        if mode == :exact
            model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
            return gsymbol_exact(model, j1, j2, j3, j4, j5, j6)
        elseif mode == :numeric
            model = NumericSU2kModel(k; T=T, prec=prec)
            return gsymbol_numeric(model, j1, j2, j3, j4, j5, j6)
        end
        throw(ArgumentError("Unknown mode: $mode"))
    end
end

