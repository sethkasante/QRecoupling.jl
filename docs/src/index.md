```@meta
CurrentModule = QRecoupling
```
# QRecoupling.jl

*Stable and Exact Evaluation of quantum recoupling coefficients and q-Hypergeometric Series via Cyclotomic Factorization.*

**QRecoupling.jl** is a high-performance Julia library for the exact, stable, and scalable evaluation of quantum recoupling coefficients and generic $q$-hypergeometric series. 

> **Notice of Migration:** `QRecoupling.jl` is the official, expanded successor to the deprecated package `QRacahSymbols.jl`. It has been completely rewritten and renamed to reflect its broader scope in topological quantum field theory (TQFT) and mathematics.

## The Core Philosophy 

The main idea is to 
> Separate Algebra from Evaluation

The numerical evaluation of highly oscillatory $q$-hypergeometric series (like the quantum $6j$-symbols) is susceptible to catastrophic cancellation and intermediate expression swell. 

**QRecoupling.jl** solves this by separating the construction of the series from its evaluation. All quantities are first evaluated symbolically into a **`Deferred Cyclotomic Representation (DCR)`**. This sparse, combinatorial data structure encodes the full algebraic skeleton of the series using exact cyclotomic factorization.

By replacing dense rational polynomial operations with highly optimized integer vector arithmetic, we eliminate intermediate precision loss. Explicit values are only obtained later by *projecting* the DCR into a chosen target field.

---
## Universal Projections

A single DCR object in `QRecoupling.jl` can be evaluated natively across multiple regimes:

* **Discrete (Root of Unity):** Fast evaluation for Turaev-Viro topological invariants ($k \in \mathbb{Z}$).
* **Exact Algebraic:** Zero-loss evaluation in cyclotomic fields $\mathbb Q(\zeta)$ using `Nemo.jl`.
* **Complex Analytic:** Efficient evaluation for any continuous deformation parameter $q \in \mathbb{C}$.
* **Classical Limit:** Exact geometric reduction to Ponzano-Regge invariants ($q \to 1$).

---

## Installation
You can install the package directly from the Julia REPL. Press `]` to enter the Pkg prompt, and run:
```julia
pkg> add QRecoupling
```

---
## Outline 
Please navigate through the documentation to learn more about the package
```@contents
Pages = [
    "series.md",
    "tqft.md",
    "tutorials/identities.md",
    "api.md"
]
Depth = 2
```
---

## Citation

If you use `QRecoupling.jl` in your research, please cite the following paper:

```bibtex
@misc{Asante2026dcr,
      title={Deferred Cyclotomic Representation for Stable and Exact Evaluation of q-Hypergeometric Series}, 
      author={Seth K. Asante},
      year={2026},
      eprint={2604.13196},
      archivePrefix={arXiv},
      primaryClass={math-ph}
}
```