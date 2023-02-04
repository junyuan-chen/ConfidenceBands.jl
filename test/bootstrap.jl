@testset "BootstrapConfidenceBand" begin
    data = testdata("testboot")
    # Results are slightly different from Matlab because quantile in Matlab works differently
    sample1 = [data.v1'; data.v2']
    # Perfectly correlated v1 and v2
    pwlb, pwub = confint(PointwiseQuantileBootBand(), sample1, level=0.95)
    lb, ub, pwlevel = confint(SuptQuantileBootBand(), sample1, level=0.95)
    @test lb ≈ pwlb
    @test ub ≈ pwub
    @test lb ≈ [0.4952267777068162, 0.9904535554136324] atol = 1e-7
    @test ub ≈ [0.5047150470461396, 1.0094300940922791] atol = 1e-7
    @test pwlevel ≈ 0.95 atol = 1e-5
    lb1, ub1, pwlevel1 = confint(SuptQuantileBootBand(1), sample1, level=0.95,
        x0=(0.95, 1-1e-3))
    @test (lb1, ub1, pwlevel1) == (lb, ub, pwlevel)
    pwlb, pwub = confint(PointwiseCVBootBand(), [0.5, 1], sample1, level=0.95)
    lb, ub, cv = confint(SuptCVBootBand(), [0.5, 1], sample1, level=0.95)
    @test lb ≈ pwlb
    @test ub ≈ pwub
    @test lb ≈ [0.4952552898980708, 0.9905105797961417] atol = 1e-7
    @test ub ≈ [0.5047447101019292, 1.0094894202038585] atol = 1e-7
    @test cv ≈ 1.6227195522876259 atol = 1e-8
    lb1, ub1, cv1 = confint(SuptCVBootBand(1), [0.5, 1], sample1, level=0.95)
    @test (lb1, ub1, cv1) == (lb, ub, cv)

    sample2 = [data.v1'; data.v3']
    # Independent v1 and v3
    pwlb, pwub = confint(PointwiseQuantileBootBand(), sample2, level=0.9759759759759662)
    lb, ub, pwlevel = confint(SuptQuantileBootBand(), sample2, level=0.95)
    @test lb ≈ pwlb
    @test ub ≈ pwub
    @test pwlevel ≈ 0.9759759759759662 atol = 1e-5
    lb, ub, pwlevel = confint(SuptQuantileBootBand(1), sample2, level=0.95)
    @test pwlevel ≈ 0.7777777777777669 atol = 1e-5
    lb, ub, cv = confint(SuptCVBootBand(), [0.5, 0.5], sample2, level=0.95)
    @test cv ≈ 1.6784807901717407 atol = 1e-7
    @test pwcoverage(lb, ub, sample2) ≈ [0.982, 0.968]
    pwlb, pwub = confint(PointwiseCVBootBand(), [0.5, 0.5], sample2, level=0.982)
    @test pwlb ≈ [0.4951010010008407, 0.008005639579010737] atol = 1e-7
    @test pwub ≈ [0.5048989989991594, 0.9919943604209893] atol = 1e-7
    lb, ub, cv = confint(SuptCVBootBand(1), [0.5, 0.5], sample2, level=0.95)
    @test cv ≈ 1.3531559426266788 atol = 1e-8

    sample3 = [data.v1'; data.v4']
    lb, ub, pwlevel = confint(SuptQuantileBootBand(), sample3, level=0.95)
    @test pwlevel ≈ 0.9714328808446455 atol = 1e-5
    lb, ub, pwlevel = confint(SuptQuantileBootBand(1), sample3, level=0.95)
    @test pwlevel ≈ 0.9142511550308007 atol = 1e-5
    lb, ub, cv = confint(SuptCVBootBand(), [0.5, 1], sample3, level=0.95)
    @test cv ≈ 1.6756228602286503 atol = 1e-7
    lb, ub, cv = confint(SuptCVBootBand(1), [0.5, 1], sample3, level=0.95)
    @test cv ≈ 1.5603502850509385 atol = 1e-7

    @test_throws ArgumentError confint(PointwiseQuantileBootBand(), sample3, level=1)
    @test_throws ArgumentError confint(PointwiseCVBootBand(), [0.5, 1], sample3, level=1)
    @test_throws ArgumentError confint(SuptQuantileBootBand(), sample3, level=1)
    @test_throws ArgumentError confint(SuptQuantileBootBand(2), sample3)
    @test_throws ArgumentError confint(SuptCVBootBand(), [0.5, 1], sample3, level=1)
    @test_throws ArgumentError confint(SuptCVBootBand(2), [0.5, 1], sample3)
end
