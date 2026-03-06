# docs/make.jl
using Documenter
using QRacahSymbols

makedocs(;
    sitename = "QRacahSymbols.jl",
    authors = "Seth K Asante <seth.kurankyi@gmail.com>",
    modules = [QRacahSymbols],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://sethkasante.github.io/QRacahSymbols.jl",
        collapselevel = 1,
        sidebar_sitename = false,
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Mathematical Framework" => "math.md",
        "API Reference" => "api.md"
    ],
    warnonly = [:missing_docs, :cross_references]
)

deploydocs(;
    repo = "github.com/sethkasante/QRacahSymbols.jl.git",
    devbranch = "main",
    push_preview = true,
)