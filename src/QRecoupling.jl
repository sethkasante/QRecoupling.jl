module QRecoupling

using Nemo, LRUCache, Base.GMP.MPZ

const Spin = Real

# --- Internal Includes ---
include("cyclotomic_types.jl")
include("projection_types.jl")  
include("admissibility.jl")    # δ and δtet checks
include("symmetries.jl")       # canonical_spins
include("assembler_dcr.jl")   
include("builders_dcr.jl") # symbols and invariants
include("project_cl.jl")
include("project_direct.jl")
include("project_exact.jl")
include("project_analytic.jl")       # Discrete, Analytic, Classical, Exact
include("direct_numeric_lse.jl")    # q6j_direct
include("eager_ex.jl")      # q6j_exact

# --- The Unified API ---

# 4. The User-Facing API
include("apis.jl")              # Unified q6j, q3j, fsymbol, etc.


# Export the high-level API
export q6j, q3j, fsymbol, gsymbol, theta_value, ExactResult, ClassicalResult

export q6j_dcr, q3j_dcr, fsymbol_dcr, gsymbol_dcr

end



# """
#     q6j(j1, j2, j3, j4, j5, j6; k=nothing, mode=:discrete, eager=false, prec=Float64)

# The master interface for the SU(2) 6j-symbol.
# Modes:
# - :discrete  (Root of unity evaluation at level k)
# - :analytic  (Generic complex q evaluation)
# - :classical (Ponzano-Regge limit q=1)
# - :exact     (Algebraic Nemo evaluation)
# """
# function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
#              k=nothing, q=nothing, mode=:discrete, eager=false, T=Float64)
    
#     # 1. Canonicalize for Caching Stability
#     # Maps all 24 tetrahedral permutations to one unique key
#     js = canonical_spins(j1, j2, j3, j4, j5, j6)
    
#     # 2. Check Admissibility
#     if mode == :discrete || mode == :exact
#         isnothing(k) && error("Level k required for $mode mode.")
#         !qδtet(js..., k) && return zero(T)
#     else
#         !δtet(js...) && return zero(T)
#     end

#     # 3. Routing Logic
#     if eager
#         # Eager path: Direct summation (Numeric LSE or Nemo)
#         if mode == :discrete
#             return q6j_direct(js..., k; T=T)
#         elseif mode == :exact
#             return q6j_exact(js..., k)
#         elseif mode == :classical
#             return project_classical_exact(q6j_dcr(js...)) # Reuse exact logic
#         end
#     else
#         # Deferred path: Build DCR and project
#         dcr = q6j_dcr(js...)
#         if mode == :discrete
#             return project_discrete(dcr, k, T)
#         elseif mode == :analytic
#             isnothing(q) && error("Complex q required for analytic mode.")
#             return project_analytic(dcr, Complex{T}(q))
#         elseif mode == :classical
#             return project_classical(dcr, T)
#         elseif mode == :exact
#             return project_exact(dcr, k)
#         end
#     end
# end