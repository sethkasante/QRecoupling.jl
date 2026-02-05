module FactoredQuantum


Transform this into an efficient code even without using logtables. We can include caching. Keep this very professional 

export QFact,QFactSum, add_vectors,
       qfactorial, qinteger,
       _qΔ2, _q6jsummand, 
       _qRacahsum, _qracah6j,
       qδ, δ, δtet, qδtet




# classical and quantum triangle conditions at level k: 
# checks admissible triple 
"""
    δ(j1, j2, j3,k) -> ::Bool
    qδ(j1, j2, j3,k) -> ::Bool

Checks the triangle conditions `j3 ≤ j1 + j2`, `j1 ≤ j2 + j3`, `j2 ≤ j3 + j1` and extra condition `j1 + j2 + j3 ≤ k` for qδ.
"""
δ(j1, j2, j3) = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)

# (j3<= j1+j2) && (j1 <= j2 + j3) && (j2 <= j3 + j1) 

qδ(j1, j2, j3, k) = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 


δtet(j1, j2, j3, j4, j5, j6) = δ(j1, j2, j3) && δ(j1, j5, j6) && δ(j2, j4, j6) && δ(j3, j4, j5)

qδtet(j1, j2, j3, j4, j5, j6) = qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)


#this represents a vector of qintegers [1],...,[n]

mutable struct  QFact
    data::Vector{Int}      # sorted multiset (qintegers)
end

QFact() = QFact(Int[])

#start from 2:n
function qfactorial(n::Int)
    QFact(ones(Int,n))
end

function add_vectors(a::Vector{Int}, b::Vector{Int})
    la, lb = length(a), length(b)

    if la ≥ lb
        c = copy(a)
        @inbounds for i in 1:lb
            c[i] += b[i]
        end
        return c
    else
        c = copy(b)
        @inbounds for i in 1:la
            c[i] += a[i]
        end
        return c
    end
end


function mul!(a::QFact, b::QFact)
    a.data = add_vectors(a.data, b.data)
    return a
end

function div!(a::QFact, b::QFact)
    a.data = add_vectors(a.data, -b.data)
    return a
end

# function Base.:*(a::QFact, b::QFact)
#     QFact(add_vectors(a.data,b.data) )
# end

# function Base.:/(a::QFact, b::QFact)
#     QFact(add_vectors(a.data,-b.data) )
# end

# Base.://(a::QFact, b::QFact) = a / b


#set q = q0^2: so that units are in q0 and monomials are in q

struct QFactSum               # global common factors
    coeffs::Vector{Int}
    terms::Vector{QFact}  # reduced monomials → coeff
end




"""
quantum triangle coefficient squared in cyclotomic form
"""
function _qΔ2(j1, j2, j3)
    δ(j1, j2, j3) || return qFactors()
    num = mul!(mul!(qfactorial(j1 + j2 - j3) , qfactorial(j1 - j2 + j3)),  
        qfactorial(-j1 + j2 + j3)) 
    den = qfactorial(j1 + j2 + j3 + 1)
    div!(num,den)
end

"""
individual elements in the Racah sum 
"""
function _q6jsummand(z,α1,α2,α3,α4,β1,β2,β3)
    num = qfactorial(z+1)
    den = mul!(qfactorial(z-α1), mul!(qfactorial(z-α2), 
            mul!( qfactorial(z-α3), mul!(qfactorial(z-α4), 
            mul!(qfactorial(β1-z), mul!(qfactorial(β2-z), 
            qfactorial(β3-z) ))))))
    return div!(num , den)
end


function _qRacahsum(α1,α2,α3,α4,β1,β2,β3)
    #range of z values
    zrange = max(α1,α2,α3,α4):min(β1,β2,β3)
    N = length(zrange)
    # n = length(zrange)
    sgns = Int[] 
    sizehint!(sgns,N)
    terms = QFact[] # add sizehint = n 
    sizehint!(terms,N)
    # unit_expr = Int[] # add sizehint = n 
    # sizehint!(unit_expr,N)
    for z in zrange
        push!(terms, _q6jsummand(z,α1,α2,α3,α4,β1,β2,β3) )
        push!(sgns, iseven(z) ? 1 : -1 )
    end
    return QFactSum(sgns,terms)
end


function _qracah6j(j1, j2, j3, j4, j5, j6)
   
    #WignerSymbols reorders this for storage purposes 
    t1 = (j1, j2, j3)
    t2 = (j1, j5, j6)
    t3 = (j2, j4, j6)
    t4 = (j3, j4, j5)

    T2 = mul!(_qΔ2(t1...), mul!(_qΔ2(t2...), mul!(_qΔ2(t3...), _qΔ2(t4...))))
    
    #WignerSymbols reorders this for storage purposes 
    α1, α2, α3, α4 = j1 + j2 + j3, j1 + j5 + j6, j2 + j4 + j6, j3 + j4 + j5
    β1, β2, β3 = j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6

    S = _qRacahsum(α1,α2,α3,α4,β1,β2,β3)
    return T2, S
end



# ----- evaluations ------- 

function eval_logqn(n::Int,k::Int)
    log(eval_qn(n,k))
end

function eval_expqn(n::Int,k::Int)
    log(eval_qn(n,k))
end

function eval_logqfact(q::QFact,k::Int)
    vec = q.data
    ss = 0.0
    for i in 1:length(vec)
        ss += vec[i]*eval_logqn(i,k)
    end
    return ss
end

function eval_qn(n::Int,k::Int)
    sinpi(n/(k+2)) / sinpi(1/(k+2)) 
end



end #module 