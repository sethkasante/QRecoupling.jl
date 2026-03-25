
# -----------------------------------------------------------
#              ----  Main APIs ----
# Master public-facing API for QRacahSymbols.jl.
# Handles mode routing, precision management, and caches.
# -----------------------------------------------------------


"""
    clear_numeric_caches!()
Empties the `LOGQFACT_CACHE` used for direct Log-Sum-Exp numeric evaluations.
"""
clear_numeric_caches!() = (empty!(LOGQFACT_CACHE); nothing)

"""
    clear_sieve_caches!()
Empties the cyclotomic sieves and analytic tracking arrays.
"""
clear_sieve_caches!() = (empty!(GLOBAL_SIEVE_CACHE); empty!(UNIT_CIRCLE_CACHE); empty!(ANALYTIC_CACHE); nothing)

"""
    clear_exact_caches!()
Empties the dense cyclotomic polynomial evaluations used by the exact Nemo engine.
Highly recommended when switching topological levels `k` to free RAM.
"""
clear_exact_caches!() = (empty!(EXACT_PHI_CACHE); nothing)

"""
    clear_caches!()
Aggressively clears all internal caches.
"""
clear_caches!() = (clear_numeric_caches!(); clear_sieve_caches!(); clear_exact_caches!(); nothing)


# Note: base_term has sign = 0 to trigger the fast-fail in downstream evaluators!
const EMPTY_CYCLO_RESULT = CycloResult(
    EMPTY_MONOMIAL, EMPTY_MONOMIAL, ZERO_MONOMIAL, CycloMonomial[], 0:-1, 0
)

# -------------------------------------
#  ---- Core 6j and 3j Symbols ---
# -------------------------------------

"""
    q6j(j1, j2, j3, j4, j5, j6, [k]; mode=:cyclo, T=Float64, prec=256)

Evaluates the quantum 6j-symbol for the SU(2)_k fusion category.
Routes computation to the algebraic (`:cyclo`), classical (`:classical`), 
rigorous cyclotomic (`:exact`), or dense numerical (`:numeric`) engines.
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::OptInt=nothing; 
             mode=nothing, T::Type{<:AbstractFloat}=Float64, prec=256)
    
    mode = isnothing(mode) ? (isnothing(k) ? :cyclo : :numeric) : mode

    # 1. k-Independent Topologies
    if mode == :cyclo || mode == :classical || mode == :classical_exact
        if !δtet(j1, j2, j3, j4, j5, j6)
            mode == :cyclo && return EMPTY_CYCLO_RESULT
            mode == :classical_exact && return ClassicalResult(0, 0//1)
            return 0.0
        end
        
        mode == :cyclo && return q6j_cyclo(j1, j2, j3, j4, j5, j6)
        mode == :classical && return q6j_classical(j1, j2, j3, j4, j5, j6)
        mode == :classical_exact && return q6j_classical_exact(j1, j2, j3, j4, j5, j6)
    end
    
    # 2. Level k Required
    isnothing(k) && throw(ArgumentError("Mode :$mode requires a topological level k."))

    # 3. Quantum Admissibility
    if !qδtet(j1, j2, j3, j4, j5, j6, k) 
        if mode == :exact
            _, z = Nemo.cyclotomic_field(2 * (k + 2), "ζ")
            return CycloExactResult(k, EMPTY_MONOMIAL, zero(z))
        else
            return zero(T)
        end
    end

    c_spins = canonical_spins(j1, j2, j3, j4, j5, j6)
    
    # 4. Engine Dispatch
    if mode == :exact
        res_cyclo = q6j_cyclo(c_spins[1], c_spins[2], c_spins[3], c_spins[4], c_spins[5], c_spins[6])
        return cyclo_to_exact(res_cyclo, k)
        
    elseif mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return _q6j_stable(model, c_spins[1], c_spins[2], c_spins[3], c_spins[4], c_spins[5], c_spins[6])
    else
        throw(ArgumentError("Unknown mode: $mode"))
    end
end


"""
    q3j(j1, j2, j3, m1, m2, [m3], [k]; mode=:cyclo, T=Float64, prec=256)

Evaluates the quantum 3j-symbol (Wigner symbol). 
If `m3` is omitted, it is automatically assumed to be `-m1 - m2`.
"""
q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin; kwargs...) = 
    q3j(j1, j2, j3, m1, m2, -m1-m2, nothing; kwargs...)

function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin, k::OptInt=nothing; 
             mode=nothing, T::Type{<:AbstractFloat}=Float64, prec=256)
    
    mode = isnothing(mode) ? (isnothing(k) ? :cyclo : :numeric) : mode

    if mode == :cyclo || mode == :classical || mode == :classical_exact
        if !δ(j1, j2, j3) || !iszero(m1 + m2 + m3)
            mode == :cyclo && return EMPTY_CYCLO_RESULT
            mode == :classical_exact && return ClassicalResult(0, 0//1)
            return 0.0
        end
        
        mode == :cyclo && return q3j_cyclo(j1, j2, j3, m1, m2, m3)
        mode == :classical && return q3j_classical(j1, j2, j3, m1, m2, m3)
        mode == :classical_exact && return q3j_classical_exact(j1, j2, j3, m1, m2, m3)
    end
    
    isnothing(k) && throw(ArgumentError("Mode :$mode requires a level k."))

    if !qδ(j1, j2, j3, k) || !iszero(m1 + m2 + m3) 
        if mode == :exact
            _, z = Nemo.cyclotomic_field(2 * (k + 2), "ζ")
            return CycloExactResult(k, EMPTY_MONOMIAL, zero(z))
        else
            return zero(T)
        end
    end

    if mode == :exact
        res_cyclo = q3j_cyclo(j1, j2, j3, m1, m2, m3)
        return cyclo_to_exact(res_cyclo, k)
        
    elseif mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return _q3j_stable(model, j1, j2, j3, m1, m2)
    else
        throw(ArgumentError("Unknown mode: $mode"))
    end
end


"""
    evaluate_cyclo(res::CycloResult, k::OptInt, [T=Float64]; prec=512, kwargs...)

Projects a deferred `CycloResult` into a concrete numerical value or exact 
algebraic field element for a given discrete level `k`.

### Arguments
- `res`: The `CycloResult` DAG to project.
- `k`: The topological level.
- `T`: The target evaluation type (defaults to `Float64`). 
  - Use `Float64`, `BigFloat`, or `Complex` for standard numeric evaluation.
  - Use `ExactResult` or `CycloExactResult` for exact algebraic evaluation in Nemo.

### Examples
    evaluate_cyclo(res, 5)                # Numeric projection at SU(2)_5 (Float64)
    evaluate_cyclo(res, 5, BigFloat)      # High-precision numeric at SU(2)_5
    evaluate_cyclo(res, 5, ExactResult)   # Exact algebraic evaluation
"""
function evaluate_cyclo(res::CycloResult, k::OptInt=nothing, ::Type{T}=Float64; q=nothing, theta=nothing, prec=512) where {T}
    #we don't check the quantum level bound here j1+j2+j3 <k  
    # ---- Exact Algebraic Evaluation -----

    if T == ExactResult || T == CycloExactResult || T == Nemo.nf_elem
        isnothing(k) && throw(ArgumentError("Exact projection requires a level k."))
        return cyclo_to_exact(res, k, T)
    end

    # ---- Numeric Evaluation -----
    if T <: Number
        if !isnothing(k)
            return evaluate_level(res, k, T; prec=prec)
        elseif !isnothing(theta)
            # q on unit circle 
            return evaluate_unit_circle(res, Float64(theta), T; prec=prec)
        elseif !isnothing(q)
            # generic q (analytic continuation)
            return evaluate_analytic(res, q, T; prec=prec)
        else
            throw(ArgumentError("Must provide a valid level k, or use kwargs q/theta."))
        end
    end

    throw(ArgumentError("Unsupported evaluation target type: $T"))
end

function evaluate_cyclo(res::CycloResult, ::Type{T}=Float64; k=nothing, q=nothing, theta=nothing, prec=512) where {T}
    return evaluate_cyclo(res, k, T; q=q, theta=theta, prec=prec)
end

# ----------------------------------------------------------------------
# Topological Tensors API (Integer, dimensions, braiding, F/G Symbols)
# ----------------------------------------------------------------------



"""
    qint(n::Int, k::Int=0; mode=:numeric, T=Float64, prec=256)

Evaluates the symmetric quantum integer [n]_q. 
Available modes: 
- `:numeric` (Hardware or arbitrary precision float)
- `:cyclo`   (Deferred sparse directed acyclic graph)
- `:exact`   (Hybrid exact Nemo cyclotomic field element)
"""
function qint(n::Int, k::OptInt=nothing; mode=nothing, T::Type{<:AbstractFloat}=Float64, prec=128)
    mode = isnothing(mode) ? (isnothing(k) ? :cyclo : :numeric) : mode

    if mode == :cyclo || mode == :classical
        mode == :cyclo && return _qint_cyclo(n) 
        return n
    end
    isnothing(k) && throw(ArgumentError("Mode :$mode requires a level k."))

    if mode == :exact
        k <= 0 && throw(ArgumentError("Mode :exact requires a valid level k > 0."))
        val = _qint_exact(n, k)
        return CycloExactResult(k, EMPTY_MONOMIAL, val)
        
    elseif mode == :numeric
        k <= 0 && throw(ArgumentError("Mode :numeric requires a valid level k > 0."))
        return _qint_numeric(n, k, T, prec)
        
    else
        throw(ArgumentError("Unknown mode: $mode. Choose :numeric, :cyclo, :exact."))
    end
end

"""
    qdim(j::Spin, [k::Int]; mode=:cyclo, T=Float64, prec=256)

Evaluates the quantum dimension `[2j+1]_q` of a spin `j` representation.
"""
function qdim(j::Spin, k::OptInt=nothing; mode=nothing, T::Type{<:AbstractFloat}=Float64, prec=128)
    mode = isnothing(mode) ? (isnothing(k) ? :cyclo : :numeric) : mode

    n = round(Int, 2j + 1)

    if mode == :cyclo || mode == :classical
        if !ishalfInt(j)
            mode == :cyclo && return ZERO_MONOMIAL
            return zero(T)
        end
        mode == :cyclo && return _qint_cyclo(n) 
        return T(n)
    end

    isnothing(k) && throw(ArgumentError("Mode :$mode requires a level k."))
    
    if !ishalfInt(j) || (2j > k)
        mode == :exact && return CycloExactResult(k, EMPTY_MONOMIAL, Nemo.QQFieldElem(0))
        return zero(T)
    end

    if mode == :exact
        return qint(n, k; mode=:exact)
    elseif mode == :numeric
        return qint(n, k; mode=:numeric, T=T,prec=prec)
    end
    error("Unknown mode: $mode")
end

"""
    rmatrix(j1, j2, j3, [k]; mode=:cyclo, T=Float64)

Evaluates the framing phase (R-Matrix) for the fusion vertex `j1 ⊗ j2 → j3`.
"""
function rmatrix(j1::Spin, j2::Spin, j3::Spin, k::OptInt=nothing; 
                 mode=nothing, T::Type{<:AbstractFloat}=Float64)
    mode = isnothing(mode) ? (isnothing(k) ? :cyclo : :numeric) : mode

    if mode == :cyclo || mode == :classical
        !(abs(j1-j2) <= j3 <= j1+j2 && isinteger(j1+j2+j3)) && return (mode == :cyclo ? ZERO_MONOMIAL : 0.0)
        mode == :cyclo && return rmatrix_cyclo(j1, j2, j3)
        return iseven(round(Int, j1 + j2 - j3)) ? 1.0 : -1.0
    end

    isnothing(k) && throw(ArgumentError("Mode :$mode requires level k."))
    
    !(abs(j1-j2) <= j3 <= min(j1+j2, k - (j1+j2))) && return (mode == :exact ? CycloExactResult(k, EMPTY_MONOMIAL, Nemo.QQFieldElem(0)) : zero(T))

    mode == :exact && return rmatrix_exact(j1, j2, j3, k)
    mode == :numeric && return rmatrix_numeric(j1, j2, j3, k; T=T)
    error("Unknown mode: $mode")
end

"""
    fsymbol(j1, j2, j3, j4, j5, j6, [k]; mode=:cyclo, T=Float64, prec=256)

Evaluates the Unitary F-Matrix element representing a change of basis 
in the fusion tree. It is derived by scaling the 6j-symbol by the 
quantum dimensions of the internal channels.
"""
function fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::OptInt=nothing; 
                 mode=nothing, T::Type{<:AbstractFloat}=Float64, prec=256)
    mode = isnothing(mode) ? (isnothing(k) ? :cyclo : :numeric) : mode

    if mode == :cyclo
        return fsymbol_cyclo(j1, j2, j3, j4, j5, j6)
    elseif mode == :exact
        isnothing(k) && throw(ArgumentError("Mode :exact requires level k."))
        return cyclo_to_exact(fsymbol_cyclo(j1, j2, j3, j4, j5, j6), k)
    elseif mode == :numeric
        isnothing(k) && throw(ArgumentError("Mode :numeric requires level k."))
        model = NumericSU2kModel(k; T=T, prec=prec)
        return fsymbol_numeric(model, j1, j2, j3, j4, j5, j6)
    elseif mode == :classical
        res_6j = q6j_classical_exact(j1, j2, j3, j4, j5, j6)
        dim_mult = (2j3 + 1) * (2j6 + 1)
        phase = iseven(round(Int, j1 + j2 + j4 + j5)) ? 1 : -1
        return ClassicalResult(res_6j.sign * phase, res_6j.sq_val * (dim_mult^2 // 1))
    end
    error("Unknown mode: $mode")
end

"""
    gsymbol(j1, j2, j3, j4, j5, j6, [k]; mode=:cyclo, T=Float64, prec=256)

Evaluates the fully symmetric G-Symbol. Restores full tetrahedral symmetry 
by completely absorbing the quantum dimensions of all 6 edges into the 6j-symbol weight.
"""
function gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::OptInt=nothing; 
                 mode=nothing, T::Type{<:AbstractFloat}=Float64, prec=256)
    mode = isnothing(mode) ? (isnothing(k) ? :cyclo : :numeric) : mode

    if mode == :cyclo
        return gsymbol_cyclo(j1, j2, j3, j4, j5, j6)
    elseif mode == :exact
        isnothing(k) && throw(ArgumentError("Mode :exact requires level k."))
        return cyclo_to_exact(gsymbol_cyclo(j1, j2, j3, j4, j5, j6), k)
    elseif mode == :numeric
        isnothing(k) && throw(ArgumentError("Mode :numeric requires level k."))
        model = NumericSU2kModel(k; T=T, prec=prec)
        return gsymbol_numeric(model, j1, j2, j3, j4, j5, j6)
    elseif mode == :classical
        res_6j = q6j_classical_exact(j1, j2, j3, j4, j5, j6)
        dim_mult = (2j1 + 1) * (2j2 + 1) * (2j3 + 1) * (2j4 + 1) * (2j5 + 1) * (2j6 + 1)
        return ClassicalResult(res_6j.sign, res_6j.sq_val * (dim_mult^2 // 1))
    end
    error("Unknown mode: $mode")
end