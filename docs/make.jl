# docs/make.jl
using Documenter
using QRecoupling

DocMeta.setdocmeta!(QRecoupling, :DocTestSetup, :(using QRecoupling); recursive=true)

makedocs(;
    sitename = "QRecoupling.jl",
    authors = "Seth K Asante <seth.kurankyi@gmail.com>",
    modules = [QRecoupling],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://sethkasante.github.io/QRecoupling.jl",
        collapselevel = 1,
        sidebar_sitename = false,
        edit_link="main",
        assets=String[],
    ),
    pages = [
        "Home" => "index.md",
        "Theory & Architecture" => "series.md",
        "Recoupling Symbols" => "tqft.md",
        "Tutorials" => [
            "Proving Identities" => "tutorials/identities.md",
        ],
        "API Reference" => "api.md"
    ],
    warnonly = [:missing_docs, :cross_references]
)

deploydocs(;
    repo = "github.com/sethkasante/QRecoupling.jl.git",
    devbranch = "main",
    push_preview = true,
)