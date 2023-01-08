using CSV
using StableRNGs

rng = StableRNG(123)

function testboot()
    N = 1000
    v1 = 0.495 .+ 0.01 .* rand(rng, N)
    v2 = 2 .* v1
    v3 = rand(rng, N)
    v4 = 2 .* v1 .+ 0.001 .* randn(rng, N)
    tb = (v1=v1, v2=v2, v3=v3, v4=v4)
    CSV.write("data/testboot.csv.gz", tb, compress=true)
    return tb
end

function main()
    testboot()
end

main()
