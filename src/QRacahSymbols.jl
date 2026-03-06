# src/QRacahSymbols.jl
module QRacahSymbols

using Nemo
using LRUCache

# Export Data Structures and Types
export CycloMonomial, ExactSU2kModel, NumericSU2kModel
export GenericResult, ExactResult
export Spin

# Export Primary API Functions
export q6j, q3j, fsymbol, rmatrix, gsymbol, qdim
export evaluate_exact, evaluate_classical

# Export Internal Handlers (Optional, for advanced users)
export qracah6j_generic, qracah6j_exact, qracah6j_numeric, qracah6j_classical
export qracah3j_generic, qracah3j_exact, qracah3j_numeric

# Include architecture files
include("Types.jl")
include("admissible.jl")
include("Symmetry.jl")
include("Symbolics.jl")
include("ExactAlgebra.jl")
include("Numerics.jl")
include("TQFT.jl")

# Initialize the global Evaluation & Model Caches
const Q6J_NUMERIC_CACHE = LRU{Tuple{NTuple{6, Float64}, Int, DataType}, Any}(maxsize=50000)
const Q6J_EXACT_CACHE   = LRU{Tuple{NTuple{6, Float64}, Int}, ExactResult}(maxsize=10000)
const EXACT_MODEL_CACHE = LRU{Int, ExactSU2kModel}(maxsize=20)


# ============================================================
# 1. Quantum 6j Symbol
# ============================================================

"""
    q6j(j1, j2, j3, j4, j5, j6; mode=:generic)

Evaluates the k-independent 6j-symbol. 
Valid modes: `:generic` (Symbolic CycloMonomials) or `:classical` (Ponzano-Regge limit).
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; mode=:generic)
    if mode == :classical
        return qracah6j_classical(j1, j2, j3, j4, j5, j6)
    elseif mode == :generic
        return qracah6j_generic(j1, j2, j3, j4, j5, j6)
    else
        error("Mode :$mode requires a level `k`. Call as `q6j(..., k; mode=:$mode)`.")
    end
end

"""
    q6j(j1, j2, j3, j4, j5, j6, k; mode=:numeric, T=Float64, prec=256)

Evaluates the quantum 6j-symbol at level k.
Valid modes: `:numeric` (Fast Float/BigFloat) or `:exact` (Nemo Field).
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    # Route k-independent requests safely back to the other dispatcher
    if mode == :generic || mode == :classical
        return q6j(j1, j2, j3, j4, j5, j6; mode=mode) 
    end

    # Admissibility Check
    if !qδtet(j1, j2, j3, j4, j5, j6, k) 
        if mode == :exact
            # Retrieve the exact model to generate the mathematically correct field zero
            model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
            return ExactResult(k, model.K(0), model.K(0))
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
        error("Unknown mode: $mode")
    end
end


# ============================================================
# 2. Quantum 3j Symbol
# ============================================================

# Auto-completes m3 if omitted (k-independent)
q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin; mode=:generic) = 
    q3j(j1, j2, j3, m1, m2, -m1-m2; mode=mode)

function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin; mode=:generic)
    if mode == :generic
        return qracah3j_generic(j1, j2, j3, m1, m2, m3)
    else
        error("Mode :$mode requires a level `k`. Call as `q3j(..., m3, k; mode=:$mode)`.")
    end
end

function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin, k::Int; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    if mode == :generic
        return q3j(j1, j2, j3, m1, m2, m3; mode=mode)
    end
    
    # Admissibility Check
    if !qδ(j1, j2, j3, k) || !iszero(m1 + m2 + m3)
        if mode == :exact
            model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
            return ExactResult(k, model.K(0), model.K(0))
        else
            return zero(T)
        end
    end

    if mode == :exact
        model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
        return qracah3j_exact(model, j1, j2, j3, m1, m2)
    elseif mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return _qracah3j_stable(model, j1, j2, j3, m1, m2)
    else
        error("Unknown mode: $mode")
    end
end

# ============================================================
# 3. Quantum Dimensions
# ============================================================

qdim(j::Spin; mode=:generic) = qdim_symb(j)

function qdim(j::Spin, k::Int; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    if mode == :generic return qdim(j; mode=:generic) end
    if mode == :exact
        # Use () -> syntax for clean one-liners
        model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
        return qdim_exact(model, j)
    elseif mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return qdim_numeric(j, model)
    end
    error("Unknown mode: $mode")
end

# ============================================================
# 4. Braiding (R-Matrix)
# ============================================================

rmatrix(j1::Spin, j2::Spin, j3::Spin; mode=:generic) = rmatrix_symb(j1, j2, j3)

function rmatrix(j1::Spin, j2::Spin, j3::Spin, k::Int; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    if mode == :generic return rmatrix(j1, j2, j3; mode=:generic) end
    if mode == :exact
        model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
        return rmatrix_exact(model, j1, j2, j3)
    elseif mode == :numeric
        return rmatrix_numeric(j1, j2, j3, k; T=T)
    end
    error("Unknown mode: $mode")
end

# ============================================================
# 5. Fusion & State Sum Weights (F-Symbol & G-Symbol)
# ============================================================

fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; mode=:generic) = 
    fsymbol_generic(j1, j2, j3, j4, j5, j6)

function fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    if mode == :generic return fsymbol(j1, j2, j3, j4, j5, j6; mode=:generic) end
    if mode == :exact
        model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
        return fsymbol_exact(model, j1, j2, j3, j4, j5, j6)
    elseif mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return fsymbol_numeric(model, j1, j2, j3, j4, j5, j6)
    end
    error("Unknown mode: $mode")
end

gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; mode=:generic) = 
    gsymbol_generic(j1, j2, j3, j4, j5, j6)

function gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    if mode == :generic return gsymbol(j1, j2, j3, j4, j5, j6; mode=:generic) end
    if mode == :exact
        model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
        return gsymbol_exact(model, j1, j2, j3, j4, j5, j6)
    elseif mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return gsymbol_numeric(model, j1, j2, j3, j4, j5, j6)
    end
    error("Unknown mode: $mode")
end

end # module