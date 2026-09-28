# Building the documentation

From the repository root, use the checkout rather than an installed release:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/check_readme.jl
julia --project=docs docs/make.jl
```

The local build writes `docs/build/` and does not deploy by default. Open
`docs/build/index.html` to read the generated pages. Set `QRECOUPLING_DOCS_BUILD_DIR` to
an absolute path for an isolated preview build. CI tests documentation on Julia 1.10
and 1.12 and explicitly enables deployment only for push builds on the 1.12 job.

Tutorials use Documenter `@example` blocks with assertions for their mathematical claims.
An example failure, broken internal cross-reference, or missing documented-export inclusion
fails the build. When editing a
numerical example, assert a justified tolerance rather than copying a fragile decimal
transcript. Use exact residual checks where the example claims an exact identity.

Keep the README, migration guide, and unreleased changelog consistent with source behavior.
Do not document development-only APIs that are no longer present. Distinguish numerical
agreement, exact coefficient arithmetic, and proof of vanishing.
