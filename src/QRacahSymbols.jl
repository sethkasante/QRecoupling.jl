# src/QRacahSymbols.jl
module QRacahSymbols

using Nemo
using LRUCache

# Data Structures and Types
export CycloMonomial, ExactSU2kModel, NumericSU2kModel
export GenericResult, ExactResult, ExactValue
export Spin

# Internal Handlers for advanced usage
export qracah6j_generic, qracah6j_exact, qracah6j_numeric, qracah6j_classical
export qracah3j_generic, qracah3j_exact, qracah3j_numeric, qracah3j_classical
export clear_caches!


# Primary API Functions 
export q6j, q3j, fsymbol, rmatrix, gsymbol, qint, qdim
export evaluate_exact, evaluate_classical, evaluate_generic

# Include architecture files
include("types.jl")
include("admissibility.jl")
include("symmetry.jl")
include("generic_engine.jl")
include("exact_engine.jl")
include("numerics_engine.jl")
include("tqft.jl")

# Initialize the global Evaluation & Model Caches
const Q6J_NUMERIC_CACHE = LRU{Tuple{NTuple{6, Float64}, Int, DataType}, Any}(maxsize=50000)
const Q6J_EXACT_CACHE   = LRU{Tuple{NTuple{6, Float64}, Int}, ExactResult}(maxsize=10000)
const EXACT_MODEL_CACHE = LRU{Int, ExactSU2kModel}(maxsize=20)
const PHI_EVAL_CACHE = LRU{Tuple{Int, Int}, Tuple{BigFloat, BigFloat}}(maxsize=10000)
const LOGQFACT_CACHE = LRU{Tuple{Int, DataType, Int}, Any}(maxsize = 10240)


"""
    clear_caches!()

Empties all internal LRU caches (numeric, exact models, and evaluations).
Useful for freeing up memory during massive state-sum computations.
"""
function clear_caches!()
    empty!(Q6J_NUMERIC_CACHE)
    empty!(Q6J_EXACT_CACHE)
    empty!(EXACT_MODEL_CACHE)
    empty!(PHI_EVAL_CACHE)
    empty!(LOGQFACT_CACHE) 
    return nothing
end

# ============================================================
# Quantum 6j Symbol
# ============================================================

"""
    q6j(j1, j2, j3, j4, j5, j6, [k]; mode=:generic, T=Float64, prec=256)

Evaluates the quantum 6j-symbol. 
If `k` is omitted, only `:generic` and `:classical` modes are valid.
"""
function q6j end

function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::OptInt=nothing; 
             mode=nothing, T::Type{<:AbstractFloat}=Float64, prec=256)
    
    # --- select default mode ---
    if mode === nothing
        mode = (k === nothing) ? :generic : :numeric
    end

    # --- k-independent computations---
    if mode == :generic || mode == :classical
        if !δtet(j1, j2, j3, j4, j5, j6)
            return mode == :generic ? GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[]) : 0.0
        end
        return mode == :generic ? qracah6j_generic(j1, j2, j3, j4, j5, j6) : qracah6j_classical(j1, j2, j3, j4, j5, j6)
    end
    
    # --- level k-dependent computations---
    if k === nothing
        throw(ArgumentError("Mode :$mode requires a level k. Try q6j(..., k; mode=:$mode)"))
    end

    if !qδtet(j1, j2, j3, j4, j5, j6, k) 
        if mode == :exact
            # zero-allocation return if possible
            if haskey(EXACT_MODEL_CACHE, k)
                return ExactResult(k, EXACT_MODEL_CACHE[k].K(0), EXACT_MODEL_CACHE[k].K(0))
            else
                K, _ = Nemo.cyclotomic_field(2*(k+2), "ζ")
                return ExactResult(k, K(0), K(0))
            end
        else
            return zero(T)
        end
    end

    c_spins = canonical_spins(j1, j2, j3, j4, j5, j6)
    
    if mode == :exact
        return get!(Q6J_EXACT_CACHE, (c_spins, k)) do
            model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
            qracah6j_exact(model, c_spins...)
        end
    elseif mode == :numeric
        return get!(Q6J_NUMERIC_CACHE, (c_spins, k, T)) do
            model = NumericSU2kModel(k; T=T, prec=prec)
            _qracah6j_stable(model, c_spins...)
        end
    else
        throw(ArgumentError("Unknown mode: $mode"))
    end
end


# ============================================================
# Quantum 3j Symbol
# ============================================================

"""
    q3j(j1, j2, j3, m1, m2, [m3], [k]; mode=:generic)

Evaluates the quantum 3j-symbol. 
If `m3` is omitted, it defaults to `-m1-m2`.
"""
function q3j end

# Helper to auto-complete m3 if it's omitted
# Note: If the user passes `k`, they MUST provide `m3` explicitly, otherwise 
# Julia will confuse `k` for `m3` due to positional argument rules.
q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin; kwargs...) = 
    q3j(j1, j2, j3, m1, m2, -m1-m2, nothing; kwargs...)


function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin, k::OptInt=nothing; 
             mode=nothing, T::Type{<:AbstractFloat}=Float64, prec=256)
    
    # --- select default mode ---
    if mode === nothing
        mode = (k === nothing) ? :generic : :numeric
    end

    # --- k-independent computations---
    if mode == :generic || mode == :classical
        if !δ(j1, j2, j3) || !iszero(m1 + m2 + m3)
            return mode == :generic ? GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[]) : 0.0
        end
        return mode == :generic ? qracah3j_generic(j1, j2, j3, m1, m2, m3) : qracah3j_classical(j1, j2, j3, m1, m2, m3)
    end
    
    if k === nothing
        throw(ArgumentError("Mode :$mode requires a level k."))
    end

    # --- level k-dependent computations---
    if !qδ(j1, j2, j3, k) || !iszero(m1 + m2 + m3) # quantum admissibility and zero checks
        return (mode == :exact ? ExactResult(k, 0, 0) : zero(T))
    end

    if mode == :exact
        model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
        return qracah3j_exact(model, j1, j2, j3, m1, m2)
    elseif mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return _qracah3j_stable(model, j1, j2, j3, m1, m2)
    else
        throw(ArgumentError("Unknown mode: $mode"))
    end
end



# ============================================================
# Quantum Dimensions
# ============================================================

"""
    qint(n::Int; mode=:generic)
    qint(n::Int, k::Int; mode=:exact, T=Float64)

Returns a representation of the quantum integer [n]_q. 
The `mode` determines the choice of representation:
- `:generic` : Returns a symbolic `CycloMonomial` factorization (requires only `n`).
- `:exact`   : Returns the exact algebraic representation in the cyclotomic field Q(ζ).
- `:numeric` : Returns the numerical floating-point value using the stable sine formula.
"""
function qint end



# ============================================================
# Public API for TQFT Symbols
# ============================================================

"""
    qdim(j, [k]; mode=:generic)
    fsymbol(j1, j2, j3, j4, j5, j6, [k]; mode=:numeric)
    rmatrix(j1, j2, j3, [k]; mode=:generic)

Access standard SU(2)_k category data.
"""
function qdim end

function fsymbol end

function rmatrix end

function gsymbol end



end # module