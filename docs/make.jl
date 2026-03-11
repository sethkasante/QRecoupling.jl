# docs/make.jl
using Documenter
using QRacahSymbols

DocMeta.setdocmeta!(QRacahSymbols, :DocTestSetup, :(using QRacahSymbols); recursive=true)

makedocs(;
    sitename = "QRacahSymbols.jl",
    authors = "Seth K Asante <seth.kurankyi@gmail.com>",
    modules = [QRacahSymbols],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://sethkasante.github.io/QRacahSymbols.jl",
        collapselevel = 1,
        sidebar_sitename = false,
        edit_link="main",
        assets=String[],
    ),
    pages = [
        "Home" => "index.md",
        "Theory & Architecture" => "math.md",
        "API Reference" => "api.md"
    ],
    warnonly = [:missing_docs, :cross_references]
)

deploydocs(;
    repo = "github.com/sethkasante/QRacahSymbols.jl.git",
    devbranch = "main",
    push_preview = true,
)