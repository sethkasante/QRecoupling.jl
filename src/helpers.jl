# Small helpers shared across the package.

# Printing utilities
const SUBSCRIPTS = Dict('0'=>'₀','1'=>'₁','2'=>'₂','3'=>'₃','4'=>'₄',
                        '5'=>'₅','6'=>'₆','7'=>'₇','8'=>'₈','9'=>'₉')

const SUPERSCRIPTS = Dict('0'=>'⁰','1'=>'¹','2'=>'²','3'=>'³','4'=>'⁴',
                          '5'=>'⁵','6'=>'⁶','7'=>'⁷','8'=>'⁸','9'=>'⁹','-'=>'⁻')

to_subscript(n::Integer)   = map(c -> SUBSCRIPTS[c], string(n))
to_superscript(n::Integer) = map(c -> SUPERSCRIPTS[c], string(n))

"Split an exponent as e = 2q + r with r in {0, 1}."
@inline function split_exp(e::Int)
    r = e & 1
    q = (e - r) >> 1
    return q, r
end

"Euler totient of a single n ≥ 1."
function _totient(n::Int)
    n < 1 && throw(DomainError(n, "totient requires n ≥ 1"))
    r = n; m = n; p = 2
    while p * p <= m
        if m % p == 0
            while m % p == 0
                m ÷= p
            end
            r -= r ÷ p
        end
        p += 1
    end
    m > 1 && (r -= r ÷ m)
    return r
end

"Validate evaluation targets and supply the classical default without overriding a level."
function _evaluation_q(k,q,exact::Bool)
    isnothing(k) || isnothing(q) || throw(ArgumentError("specify k or q, not both"))
    if !isnothing(k)
        if k isa AbstractVector
            all(x -> x isa Integer && x >= 0,k) || throw(ArgumentError("levels must be nonnegative integers"))
        else
            k isa Integer || throw(ArgumentError("k must be an integer level"))
            k >= 0 || throw(DomainError(k,"level must be nonnegative"))
        end
        return nothing
    end
    qq = isnothing(q) ? 1 : q
    exact && qq != 1 && throw(ArgumentError("exact evaluation requires k or the classical target q=1"))
    return qq
end
