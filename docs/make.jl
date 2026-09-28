# docs/make.jl
using Documenter
using QRecoupling

DocMeta.setdocmeta!(QRecoupling, :DocTestSetup, :(using QRecoupling); recursive=true)

makedocs(;
    sitename = "QRecoupling.jl",
    authors = "Seth K Asante <seth.kurankyi@gmail.com>",
    modules = [QRecoupling],
    build = get(ENV, "QRECOUPLING_DOCS_BUILD_DIR", "build"),
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
        "Getting Started" => "getting_started.md",
        "Recoupling Symbols" => "tqft.md",
        "Applications" => "applications.md",
        "Tutorials" => [
            "Tensor Networks & Modular Data" => "tutorials/tensor_networks.md",
            "Exact Forms" => "tutorials/exact_forms.md",
            "Finite Series" => "tutorials/finite_series.md",
            "Checking & Proving Identities" => "tutorials/identities.md",
        ],
        "Factorial Rules & Architecture" => "series.md",
        "Accuracy & Performance" => "performance.md",
        "API Reference" => "api.md",
        "Migrating to v0.4" => "migration.md"
    ],
    checkdocs = :exports
)

# Local builds execute examples without publishing; CI enables deployment explicitly.
if get(ENV, "QRECOUPLING_DOCS_DEPLOY", "false") == "true"
    deploydocs(;
        repo = "github.com/sethkasante/QRecoupling.jl.git",
        devbranch = "main",
        push_preview = false,
    )
end
