module QRecoupling

using Nemo, LRUCache

const Spin = Union{Int, Rational, Float64}

# --- Internal Includes ---
include("types_cyclotomic.jl")
include("types_projection.jl")  
include("admissibility.jl")    
include("symmetries.jl")    

include("generic_series_dcr.jl")   
include("recoupling_dcr.jl") # recoupling symbols and tqft invariants

include("proj_classical.jl")
include("projection_discrete.jl")
include("projection_exact.jl")
include("projection_analytic.jl")       
include("eager_discrete.jl")  
include("eager_exact.jl")    

include("recoupling_api.jl")  


# Export physics tqft and recoupling symbols api
export q6j, q3j, fsymbol, gsymbol, rmatrix, qdim

#projection 
export project_discrete, project_exact, project_analytic
export project_classical, project_classical_exact


# Export api for generic series
export CyclotomicMonomial, DCR
export add_qint!, add_qfact!, build_dcr, build_dcr!, project_dcr

end
