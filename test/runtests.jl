using Test
using ConfidenceBands

using CSV
using CodecZlib
using ConfidenceBands: RandnPool, _globalrandnpool
using GLM

const tests = [
    "plugin",
    "bootstrap"
]

function testdata(name)
    path = (@__DIR__)*"/../data/$name.csv.gz"
    open(path) |> GzipDecompressorStream |> read |> CSV.File
end

printstyled("Running tests:\n", color=:blue, bold=true)

@time for test in tests
    include("$test.jl")
    println("\033[1m\033[32mPASSED\033[0m: $(test)")
end
