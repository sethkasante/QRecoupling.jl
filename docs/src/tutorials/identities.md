
# Proving Topological Identities

Another advantage of `QRecoupling.jl` is its ability to verify axioms of 3D topological quantum field theories (TQFTs). 

In this tutorial, we will prove the two most important topological identities in the Turaev-Viro/Ponzano-Regge state sum models: The **Orthogonality Relation** (Bubble Move) and the **Biedenharn-Elliott Pentagon Identity** (Pachner 2-3 Move).

---

### Proof of Orthogonality 

In a triangulated 3-manifold, integrating out an internal edge corresponds to the "Bubble Move". Mathematically, this enforces the orthogonality of the $6j$-symbols.

We can prove this *exactly* (without floating-point errors) using the `exact=true` condition.

```julia
using QRecoupling

function test_orthogonality(j1, j2, j3, j4, j5, j6;k=k)
    # determine valid bounds for x
    x_min = max(abs(j1 - j2), abs(j3 - j4))
    x_max = min(j1 + j2, j3 + j4)
    
    LHS = 0
    
    for x in x_min:x_max
        # sum over admissible spins 
        if (j1 + j2 + x) <= k && (j3 + j4 + x) <= k
            dim_x = qdim(x, k=k, exact=true)
            sym1  = q6j(j1, j2, x, j3, j4, j5, k=k, exact=true)
            sym2  = q6j(j1, j2, x, j3, j4, j6, k=k, exact=true)
            
            LHS += dim_x * sym1 * sym2
        end
    end
    
    # evaluate exact RHS
    RHS = 0
    if j5 == j6 && abs(j1 - j2) <= j5 <= (j1 + j2) && (j1 + j2 + j5) <= k && abs(j3 - j4) <= j5 <= (j3 + j4) && (j3 + j4 + j5) <= k
        RHS = 1/qdim(j5, k=k, exact=true)
    end
    
    diff = LHS-RHS

    println("Orthogonality Check for SU(2)_$k")
    println("--------------------------------")
    println("LHS (Sum): ", LHS)
    println("RHS (Exact): ", RHS)
    println("LHS - RHS = ", diff)
    println(iszero(diff)  ? "Identity Holds!" : "Identity Failed!")
end


julia> test_orthogonality(1,1,1,1,1,1, k=5)
Orthogonality Check for SU(2)_5
--------------------------------
LHS (Sum): Exact Algebraic Result in ℚ(ζ₁₄):
  Value: (-ζ^4 + ζ^3)
RHS (Exact): -ζ^4 + ζ^3
LHS - RHS = Exact Algebraic Result in ℚ(ζ₁₄):
  Value: 0
Identity Holds!

julia> test_orthogonality(0.5, 0.5, 0.5, 0.5, 1.0, 1.0, k=10)
Orthogonality Check for SU(2)_10
--------------------------------
LHS (Sum): Exact Algebraic Result in ℚ(ζ₂₄):
  Value: (-1//2*ζ^6 + ζ^2 - 1//2)
RHS (Exact): -1//2*ζ^6 + ζ^2 - 1//2
LHS - RHS = Exact Algebraic Result in ℚ(ζ₂₄):
  Value: 0
Identity Holds!


julia> test_orthogonality(12,15,17,18,13,14, k=60)
Orthogonality Check for SU(2)_60
--------------------------------
LHS (Sum): Exact Algebraic Result in ℚ(ζ₁₂₄):
  Value: 0
RHS (Exact): 0
LHS - RHS = Exact Algebraic Result in ℚ(ζ₁₂₄):
  Value: 0
Identity Holds!
```
By leveraging `QRecoupling.jl`, the radical prefactors of the $6j$-symbols perfectly annihilate each other during multiplication, keeping the entire computation division-free and strictly inside the cyclotomic field $\mathbb{Q}(\zeta)$.

--- 

### The Biedenharn-Elliott (Pentagon) Identity
The Pentagon Identity guarantees that topological invariants are independent of the chosen triangulation. It relates the product of two $6j$-symbols to a sum over the product of three $6j$-symbols.

Because this involves a large summation and many multiplications, it is the perfect candidate for our high-speed `:numeric` engine.
```julia
using QRecoupling

#TODO: test  ∑ {6j}{6j}{6j} = {6j}{6j} 
```