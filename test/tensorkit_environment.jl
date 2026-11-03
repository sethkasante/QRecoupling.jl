# Run with `julia test/tensorkit_environment.jl`. Resolve an environment with TensorKit and without Nemo, then test the
# sector extension in a fresh process with no fallback to the user's global environment.
using Pkg
mktempdir() do env
    Pkg.activate(env)
    Pkg.develop(path = dirname(@__DIR__))
    Pkg.add(["TensorKit", "TensorKitSectors", "TensorOperations", "TestExtras", "Test", "Random", "LinearAlgebra"])
    @assert !any(p -> p.name == "Nemo", values(Pkg.dependencies()))
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$env $(joinpath(@__DIR__, "tensorkit_tests.jl"))`
    run(addenv(cmd, "JULIA_LOAD_PATH" => join(("@", "@stdlib"), Sys.iswindows() ? ';' : ':')))
end
