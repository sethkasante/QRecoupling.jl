
# Proving Topological Identities

The true power of `QRecoupling.jl` lies in its ability to effortlessly verify the foundational axioms of 3D topological quantum field theories (TQFTs). 

In this tutorial, we will use the package's multi-paradigm engines to prove the two most important topological identities in the Turaev-Viro/Ponzano-Regge state sum models: The **Orthogonality Relation** (Bubble Move) and the **Biedenharn-Elliott Pentagon Identity** (Pachner 2-3 Move).

---

## 1. Exact Orthogonality (The Bubble Move)

In a triangulated 3-manifold, integrating out an internal edge corresponds to the "Bubble Move" (or Pachner 1-4 move). Mathematically, this enforces the orthogonality of the $6j$-symbols.

We can prove this *exactly* (without floating-point errors) using the `:exact` cyclotomic engine.

```julia
using QRecoupling
using Nemo

# Define spins and topological level
j1, j2, j3, j4, j5 = 1, 1, 1, 1, 1
k = 6

# Set up the exact cyclotomic field for SU(2)_k
h = k + 2
K, z = cyclotomic_field(2 * h, "ζ")

# Accumulate the sum: Σ_x [2x+1]_q * {j1 j2 j3; j4 j5 x}^2
sum_val = K(0)

x_min = max(abs(j1 - j5), abs(j2 - j4))
x_max = min(j1 + j5, j2 + j4, k - max(j1+j5, j2+j4))

for x in x_min:1:x_max
    # Compute the exact symbols
    res_A = q6j(j1, j2, j3, j4, j5, x, k; mode=:exact)
    res_B = q6j(j1, j2, j3, j4, j5, x, k; mode=:exact)
    
    # Square the symbol and multiply by the quantum dimension
    squared_symbol = res_A * res_B 
    dim_x = qdim(x, k; mode=:exact)
    
    # Add to the running total (extracting factor from the ExactResult)
    sum_val += (squared_symbol * dim_x).factor
end

# The exact right-hand side is 1 / [2j3+1]_q
dim_j3 = qdim(j3, k; mode=:exact);
expected_rhs = inv(dim_j3.factor);

julia> println("Bubble Sum evaluates to: ", sum_val)
Bubble Sum evaluates to: Exact SU(2)₆ Symbol:
  Value: -ζ^6 + ζ^2 - 1

julia> println("Theoretical RHS is: ", expected_rhs)
Theoretical RHS is: -ζ^6 + ζ^2 - 1


julia> println("Exact Match?: ", sum_val.factor == expected_rhs)
Exact Match?: true
```
By leveraging `QRecoupling.jl`, the radical prefactors of the $6j$-symbols perfectly annihilate each other during multiplication, keeping the entire computation division-free and strictly inside the cyclotomic field $\mathbb{Q}(\zeta)$.

--- 

## 2. The Biedenharn-Elliott Pentagon Identity
The Pentagon Identity guarantees that topological invariants are independent of the chosen triangulation. It relates the product of two $6j$-symbols to a sum over the product of three $6j$-symbols.

Because this involves a large summation and many multiplications, it is the perfect candidate for our high-speed `:numeric` engine.
```julia
using QRecoupling

j1, j2, j3 = 1, 1, 1
l1, l2, l3 = 1, 1, 1
j23, j12 = 1, 1
k = 8

# Left Hand Side: Σ_x (-1)^Φ [2x+1]_q { } { } { }
lhs_sum = 0.0

for x in 0:k
    # The numeric engine automatically returns 0.0 for inadmissible bounds
    w1 = q6j(j1, j2, j12, l1, l2, x, k; mode=:numeric)
    w2 = q6j(j1, x, l2, l3, j3, j23, k; mode=:numeric)
    w3 = q6j(l1, j2, x, j3, l3, l2, k; mode=:numeric)
    
    if w1 != 0.0 && w2 != 0.0 && w3 != 0.0
        dim_x = qdim(x, k; mode=:numeric)
        
        # Topological phase shift
        phase = iseven(round(Int, j12 + j23 + x + l2)) ? 1.0 : -1.0
        
        lhs_sum += phase * dim_x * w1 * w2 * w3
    end
end

# Right Hand Side: { } * { }
rhs_w1 = q6j(j12, j3, j23, l3, l1, l2, k; mode=:numeric)
rhs_w2 = q6j(j1, j2, j12, j3, j23, l1, k; mode=:numeric)
rhs_val = rhs_w1 * rhs_w2

println("Pentagon LHS : ", lhs_sum)
println("Pentagon RHS : ", rhs_val)
println("Difference   : ", abs(lhs_sum - rhs_val))
```