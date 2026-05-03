# API Reference

```@meta
CurrentModule = QRecoupling
```
This page details the public API for `QRecoupling.jl`. Internal evaluations are abstracted away to provide a clean, unified interface.

## Core Recoupling Symbols
The master dispatchers for evaluating $3j$ and $6j$ symbols. The computational regime is dynamically routed based on the provided target: integer level $k$ (Turaev-Viro), parameter $q$ (complex analytic), or $q=1$ (classical Ponzano-Regge limit). Precision is controlled via the exact boolean flag.
```@docs
q6j
q3j
```

---

### Other TQFT Tensors
Composite tensors and network invariants used in the construction of 3D TQFTs, spin foam models, and string-net models.
```@docs
qdim
rmatrix
fsymbol
gsymbol
```
---

### $q$-Hypergeometric Primitives
Public primitives to easily build algebraic representations of $q$-deformed integers, factorials, and binomial coefficients. These functions return `CyclotomicMonomial` objects and are designed to be used inside the `qseries` sequence builder.
```@docs
qint
qfact
qbinomial
```

---

## Generic Series & Universal Evaluation
`QRecoupling.jl` separates algebraic construction from field evaluation. These functions allow you to construct custom $q$-hypergeometric series and project abstract Deferred Cyclotomic Representation (DCR) objects into concrete target fields.
```@docs
qseries
build_series
build_dcr!
qeval
CyclotomicMonomial
DCR
```

--- 

### Low-Level Buffer Operations & Projections
For developers building advanced, memory-optimized loops or requesting specific projection regimes directly without routing through `qeval`.
```@docs
add_qint!
add_qfact!
project_discrete
project_exact
project_classical
project_classical_exact
project_analytic
```

---

## Memory & Cache Management
When changing the topological level $k$ drastically, performing exact computations in cyclotomic fields, or benchmarking tight loops, it is recommended to clear these caches to free RAM.
```@docs
empty_caches!
```