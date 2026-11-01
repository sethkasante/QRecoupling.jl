# Run with `julia test/core_environment.jl`. Resolve a genuinely Nemo-free environment,
# then test it in a fresh process with no fallback to the user's global environment.
using Pkg
mktempdir() do env
    Pkg.activate(env)
    Pkg.develop(path = dirname(@__DIR__))
    @assert !any(p -> p.name == "Nemo", values(Pkg.dependencies()))
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$env $(joinpath(@__DIR__, "core_tests.jl"))`
    run(addenv(cmd, "JULIA_LOAD_PATH" => join(("@", "@stdlib"), Sys.iswindows() ? ';' : ':')))
end
