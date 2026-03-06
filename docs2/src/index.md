# QRacahSymbols.jl Documentation

Welcome to the documentation for **QRacahSymbols.jl**, a state-of-the-art computational engine for $SU(2)_k$ topological quantum field theory data.

## Overview

The calculation of quantum $6j$-symbols relies on the $q$-deformed Kirillov-Reshetikhin formula, an alternating hypergeometric sum over quantum factorials. For large spins and high levels, standard computational libraries inevitably fail due to catastrophic cancellation.

`QRacahSymbols.jl` resolves this by exposing a rigorous **Three-Tier Architecture**:

1. **`mode=:numeric`**: Employs an $O(1)$ table-lookup system fortified by a `log-sum-exp` shift, allowing for nanosecond evaluation of tensor network invariants while remaining perfectly stable up to massive spin bounds.
2. **`mode=:exact`**: Leverages `Nemo.jl` to map the evaluation directly into the cyclotomic number field $\mathbb{Q}(\zeta_{2N})$, guaranteeing zero precision loss and returning exact topological invariants suitable for formal mathematical proofs.
3. **`mode=:generic`**: A highly optimized prime-factorization vector engine that returns the abstract polynomial decomposition of the symbol into cyclotomic roots ($\Phi_d(q)$).

`QRacahSymbols.jl` resolves this by exposing a rigorous **Three-Tier Architecture**:
1. **Numeric**: Employs an $O(1)$ table-lookup system fortified by a `log-sum-exp` shift, allowing for nanosecond evaluation of tensor network invariants.
2. **Exact**: Leverages `Nemo.jl` to map the evaluation directly into the cyclotomic number field $\mathbb{Q}(\zeta_{2N})$, guaranteeing zero precision loss.
3. **Generic**: A highly optimized prime-factorization vector engine that returns the abstract polynomial decomposition.

## Next Steps

* Visit the [Mathematical Framework](@ref) page to understand the symmetry canonicalization and hypergeometric ratio optimizations utilized under the hood.
* Visit the [API Reference](@ref) page for exhaustive details on function signatures, TQFT evaluators, and type stabilization strategies.