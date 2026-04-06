module QRecoupling

using Nemo, LRUCache

const Spin = Real

# --- Internal Includes ---
include("types_cyclotomic.jl")
include("types_projection.jl")  
include("admissibility.jl")    
include("symmetries.jl")    
include("assembler_dcr.jl")   
include("builders_dcr.jl") # symbols and tqft invariants
include("projection_classical.jl")
include("projection_direct.jl")
include("projection_exact.jl")
include("projection_analytic.jl")       
include("eager_discrete.jl")  
include("eager_exact.jl")    
include("collection_apis.jl")  


# Export the high-level API
export q6j, q3j, fsymbol, gsymbol, tetrahedron, theta_value
export ExactResult, ClassicalResult, DCR, CyclotomicMonomial
export q6j_dcr, q3j_dcr, fsymbol_dcr, gsymbol_dcr
export qdim, rmatrix, rmatrix_mono
export evaluate_exact, project_analytic, project_discrete
export project_classical, project_classical_exact




end
