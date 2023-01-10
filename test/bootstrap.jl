@testset "BootstrapConfidenceBand" begin
    data = testdata("testboot")
    # Results are slightly different from Matlab because quantile in Matlab works differently
    sample1 = [data.v1'; data.v2']
    # Perfectly correlated v1 and v2
    lb, ub, pwlevel = confint(SuptQuantileBootBand(), sample1, level=0.95)
    @test lb ≈ [0.4952267777068162, 0.9904535554136324] atol = 1e-7
    @test ub ≈ [0.5047150470461396, 1.0094300940922791] atol = 1e-7
    @test pwlevel ≈ 0.95 atol = 1e-5
    lb1, ub1, pwlevel1 = confint(SuptQuantileBootBand(1), sample1, level=0.95,
        x0=(0.95, 1-1e-3))
    @test (lb1, ub1, pwlevel1) == (lb, ub, pwlevel)
    lb, ub, cv = confint(SuptCVBootBand(), [0.5, 1], sample1, level=0.95)
    @test lb ≈ [0.4997396923218894, 0.9994793846437788] atol = 1e-7
    @test ub ≈ [0.5002603076781106, 1.0005206153562212] atol = 1e-7
    @test cv ≈ 0.08902680033261515 atol = 1e-8
    lb1, ub1, cv1 = confint(SuptCVBootBand(1), [0.5, 1], sample1, level=0.95)
    @test (lb1, ub1, cv1) == (lb, ub, cv)

    sample2 = [data.v1'; data.v3']
    # Independent v1 and v3
    lb, ub, pwlevel = confint(SuptQuantileBootBand(), sample2, level=0.95)
    @test pwlevel ≈ 0.9759759759759662 atol = 1e-5
    lb, ub, pwlevel = confint(SuptQuantileBootBand(1), sample2, level=0.95)
    @test pwlevel ≈ 0.7777777777777669 atol = 1e-5
    lb, ub, cv = confint(SuptCVBootBand(), [0.5, 1], sample2, level=0.95)
    @test cv ≈ 0.5350436800458136 atol = 1e-7
    lb, ub, cv = confint(SuptCVBootBand(1), [0.5, 1], sample2, level=0.95)
    @test cv ≈ 0.053233373674966104 atol = 1e-8

    sample3 = [data.v1'; data.v4']
    lb, ub, pwlevel = confint(SuptQuantileBootBand(), sample3, level=0.95)
    @test pwlevel ≈ 0.9714328808446455 atol = 1e-5
    lb, ub, pwlevel = confint(SuptQuantileBootBand(1), sample3, level=0.95)
    @test pwlevel ≈ 0.9142511550308007 atol = 1e-5
    lb, ub, cv = confint(SuptCVBootBand(), [0.5, 1], sample3, level=0.95)
    @test cv ≈ 0.14906734651562314 atol = 1e-7
    lb, ub, cv = confint(SuptCVBootBand(1), [0.5, 1], sample3, level=0.95)
    @test cv ≈ 0.05326666952510811 atol = 1e-7

    @test_throws ArgumentError confint(SuptQuantileBootBand(), sample3, level=1)
    @test_throws ArgumentError confint(SuptQuantileBootBand(2), sample3)
    @test_throws ArgumentError confint(SuptCVBootBand(), [0.5, 1], sample3, level=1)
    @test_throws ArgumentError confint(SuptCVBootBand(2), [0.5, 1], sample3)
end
