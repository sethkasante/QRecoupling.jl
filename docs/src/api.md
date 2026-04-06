# API Reference

```@meta
CurrentModule = QRacahSymbols
```
This page details the public API for `QRacahSymbols.jl`. Internal evaluation engines and algorithmic caches are heavily abstracted away to provide a clean, unified interface.

## Core Recoupling Symbols
The master dispatchers for evaluating $3j$ and $6j$ symbols. By altering the `mode` keyword, users can swap the underlying computational architecture.
```@docs
q6j
q3j
```

---

## Topological Category Tensors
Tensors used directly in the construction of 3D topological quantum field theories (TQFTs) and string-net models.
```@docs
qdim
rmatrix
fsymbol
gsymbol
```
---

## Advanced Evaluation
To separate the topological construction phase from the numerical evaluation phase, these functions allow you to project deferred CycloResult objects into numerical regimes.
```@docs
evaluate_cyclo
```

---

## Memory & Cache Management
`QRacahSymbols.jl` is designed for massive state sums and aggressively caches prime factorizations, integer phases, and dense logarithmic arrays in memory. When changing the topological level `k` drastically, it is recommended to clear these caches.
```@docs
clear_caches!
clear_numeric_caches!
clear_exact_caches!
```