module FactoredQuantum

using Primes
using LRUCache

export CyclotomicBasis,
    CycloMonomial, CycloSum, 
    qinteger, qfactorial,
    qΔ2, qδ, δ, 
    DIVISOR_CACHE, QINTEGER_CACHE, QFACTORIAL_CACHE,
    q6jsummand, qRacahsum, factored_sum,
    _qracah6j, cyclotomic_eval, q_root,
    symbolic_expression

"""
    CycloMonomial

A monomial in the cyclotomic polynomial basis.
Represents ∏ Φ_d(q)^{e_d} where Φ_d is the d-th cyclotomic polynomial.
"""
struct CycloMonomial
    powers::Dict{Int, Int}   # d => exponent of Φ_d
    
    function CycloMonomial(powers::Dict{Int, Int})
     # Filter out zero exponents for type stability
     filtered = Dict(d => e for (d, e) in powers if e != 0)
     new(filtered)
    end
end

CycloMonomial() = CycloMonomial(Dict{Int, Int}())
CycloMonomial(d::Int, e::Int=1) = CycloMonomial(Dict(d => e))

const DIVISOR_CACHE = LRU{Int, Vector{Int}}(maxsize=10_000)
const QINTEGER_CACHE = LRU{Int, CycloMonomial}(maxsize=10_000)
const QFACTORIAL_CACHE = LRU{Int, CycloMonomial}(maxsize=5_000)

@inline function divisors_cached(n::Int)::Vector{Int}
    n == 0 && return Int[]
    return get!(DIVISOR_CACHE, n) do
     divisors(n)
    end
end

function Base.:*(a::CycloMonomial, b::CycloMonomial)::CycloMonomial
    expo = copy(a.powers)
    for (d, e) in b.powers
     expo[d] = get(expo, d, 0) + e
    end
    CycloMonomial(expo)
end

function Base.:/(a::CycloMonomial, b::CycloMonomial)::CycloMonomial
    a * CycloMonomial(Dict(d => -e for (d, e) in b.powers))
end

Base.://(a::CycloMonomial, b::CycloMonomial)::CycloMonomial = a / b

function Base.:^(a::CycloMonomial, n::Int)::CycloMonomial
    CycloMonomial(Dict(d => n * e for (d, e) in a.powers))
end

function Base.:(==)(a::CycloMonomial, b::CycloMonomial)::Bool
    a.powers == b.powers
end

function Base.hash(a::CycloMonomial, h::UInt)::UInt
    hash(a.powers, h)
end

"""
    qinteger(n::Int)::CycloMonomial

Returns cyclotomic representation of [n]_q.
"""
function qinteger(n::Int)::CycloMonomial
    get!(QINTEGER_CACHE, n) do
     divs = divisors_cached(n)[2:end]  # omit d = 1
     CycloMonomial(Dict(divs .=> 1))
    end
end

"""
    qfactorial(n::Int)::CycloMonomial

Returns cyclotomic representation of [n]_q!.
"""
function qfactorial(n::Int)::CycloMonomial
    n <= 1 && return CycloMonomial()
    get!(QFACTORIAL_CACHE, n) do
     qfactorial(n - 1) * qinteger(n)
    end
end

"""
    δ(j1, j2, j3)::Bool

Checks classical triangle conditions.
"""
function δ(j1, j2, j3)::Bool
    isinteger(j1 + j2 + j3) && (j3 <= j1 + j2) && (j1 <= j2 + j3) && (j2 <= j3 + j1)
end

"""
    qδ(j1, j2, j3, k)::Bool

Checks quantum triangle conditions at level k.
"""
function qδ(j1, j2, j3, k)::Bool
    δ(j1, j2, j3) && (j1 + j2 + j3) <= k
end

function δtet(j1, j2, j3, j4, j5, j6)::Bool
    δ(j1, j2, j3) && δ(j1, j5, j6) && δ(j2, j4, j6) && δ(j3, j4, j5)
end

"""
    qΔ2(j1, j2, j3)::CycloMonomial

Quantum triangle coefficient in cyclotomic representation.
"""
function qΔ2(j1, j2, j3)::CycloMonomial
    qfactorial(j1 + j2 - j3) * qfactorial(j1 - j2 + j3) * 
     qfactorial(-j1 + j2 + j3) / qfactorial(j1 + j2 + j3 + 1)
end

"""
    CycloSum

Represents a sum ∑ (-1)^{z_i} R_i where common factors are factored out.
"""
struct CycloSum
    common_factor::CycloMonomial
    terms::Vector{CycloMonomial}
    coeffs::Vector{Int}
end

"""
    q6jsummand(z, α1, α2, α3, α4, β1, β2, β3)::CycloMonomial

Individual summand in the Racah sum.
"""
function q6jsummand(z, α1, α2, α3, α4, β1, β2, β3)::CycloMonomial
    num = qfactorial(z + 1)
    den = qfactorial(z - α1) * qfactorial(z - α2) * qfactorial(z - α3) * 
       qfactorial(z - α4) * qfactorial(β1 - z) * qfactorial(β2 - z) *
       qfactorial(β3 - z)
    num / den
end

"""
    common_factors(terms::Vector{CycloMonomial})::CycloMonomial

Extract common factors from a list of monomials.
"""
function common_factors(terms::Vector{CycloMonomial})::CycloMonomial
    isempty(terms) && return CycloMonomial()
    
    common = Dict{Int, Int}()
    all_d = Set{Int}()
    
    for m in terms
     union!(all_d, keys(m.powers))
    end
    
    for d in all_d
     cmin = minimum(get(m.powers, d, 0) for m in terms)
     if cmin > 0
         common[d] = cmin
     end
    end
    
    CycloMonomial(common)
end

"""
    factored_sum(terms::Vector{CycloMonomial})::Tuple{CycloMonomial, Vector{CycloMonomial}}

Factor out common factors from terms.
"""
function factored_sum(terms::Vector{CycloMonomial})::Tuple{CycloMonomial, Vector{CycloMonomial}}
    cm = common_factors(terms)
    residue = isempty(cm.powers) ? terms : [m / cm for m in terms]
    cm, residue
end

"""
    qRacahsum(α1, α2, α3, α4, β1, β2, β3)::Tuple{CycloMonomial, Vector{CycloMonomial}, Vector{Int}}

Racah sum as a sum over cyclotomic polynomial basis.
"""
function qRacahsum(α1, α2, α3, α4, β1, β2, β3)::Tuple{CycloMonomial, Vector{CycloMonomial}, Vector{Int}}
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3)
    sgns = Int[]
    terms = CycloMonomial[]
    
    sizehint!(terms, length(zrange))
    sizehint!(sgns, length(zrange))
    
    for z in zrange
     push!(terms, q6jsummand(z, α1, α2, α3, α4, β1, β2, β3))
     push!(sgns, iseven(z) ? 1 : -1)
    end
    
    cm, res = factored_sum(terms)
    cm, res, sgns
end

"""
    _qracah6j(j1, j2, j3, j4, j5, j6)::Tuple{CycloMonomial, CycloSum}

Compute the 6j symbol in cyclotomic form.
"""
function _qracah6j(j1, j2, j3, j4, j5, j6)::Tuple{CycloMonomial, CycloSum}
    !δtet(j1, j2, j3, j4, j5, j6) && throw(ArgumentError("Inadmissible spin values"))
    
    α1, α2, α3, α4 = j1 + j2 + j3, j1 + j5 + j6, j2 + j4 + j6, j3 + j4 + j5
    β1, β2, β3 = j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6
    
    tri2_coeff = qΔ2(j1, j2, j3) * qΔ2(j1, j5, j6) * qΔ2(j2, j4, j6) * qΔ2(j3, j4, j5)
    cm, res, sgns = qRacahsum(α1, α2, α3, α4, β1, β2, β3)
    
    R = tri2_coeff * cm^2
    S = CycloSum(R, res, sgns)
    
    R, S
end


"""
    symbolic_expression(m::CycloMonomial; var::String="q")::String

Return a symbolic expression for a cyclotomic monomial.
Example: Φ_2(q)³ * Φ_5(q)
"""
function symbolic_expression(m::CycloMonomial; var::String="q")::String
    isempty(m.powers) && return "1"
    
    factors = String[]
    for d in sort(collect(keys(m.powers)))
     e = m.powers[d]
     if e == 1
         push!(factors, "Φ_$d($var)")
     else
         push!(factors, "Φ_$d($var)^$e")
     end
    end
    
    join(factors, " * ")
end

end  # end module
