# Execute the README's Julia examples; installing the package is intentionally excluded.
using QRecoupling

readme = joinpath(@__DIR__, "..", "README.md")
examples = Module(:ReadmeExamples)
count = 0
for block in eachmatch(r"```julia\n(.*?)```"s, read(readme, String))
    code = block.captures[1]
    startswith(strip(code), "using Pkg") && continue
    include_string(examples, code, readme)
    global count += 1
end
println("README examples passed: ", count, " blocks")
