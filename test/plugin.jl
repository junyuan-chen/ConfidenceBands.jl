@testset "PlugInConfidenceBand" begin
    # Compare critical values with those from Matlab
    # [Pwise,Sidak,Bonferroni,thetaproj,muproj] = SimInference.critvalues(0.9, 10, 5)
    @test criticalvalue(PointwiseBand(), 0.9, 5) ≈ 1.644853626951472 atol = 1e-10
    @test criticalvalue(SidakBand(), 0.9, 5) ≈ 2.310660084280671 atol = 1e-10
    @test criticalvalue(BonferroniBand(), 0.9, 5) ≈ 2.326347874040841 atol = 1e-10
    @test criticalvalue(ProjectionBand(), 0.9, 5) ≈ 3.039137525644590 atol = 1e-10
    @test criticalvalue(ProjectionBand(10), 0.9, 5) ≈ 3.998397075342225 atol = 1e-10

    rp = RandnPool(10)
    @test length(rp.v) == 10

    Σ = [0.656902 0.61502 0.125656;
        0.870256 0.358389 0.239537;
        0.52512  0.771079 0.287328]
    Σ = Σ'Σ
    @test length(_globalrandnpool.v) == 40_000_000
    # Compare with SimInference.suptcritval_plugin(0.9, S, 1e8)
    @test criticalvalue(SuptBand(), 0.9, Σ) ≈ 1.862 atol = 1e-2
    @test criticalvalue(SuptBand(ndraw=15_000_000), 0.9, Σ) ≈ 1.862 atol = 5e-3
    @test length(_globalrandnpool.v) == 45_000_000
    resize!(_globalrandnpool.v, 40_000_000)
    @test criticalvalue(SuptBand(1), 0.9, Σ) ≈ 1.252 atol = 1e-2
    @test_throws ArgumentError criticalvalue(SuptBand(), 1, Σ)
    @test_throws ArgumentError criticalvalue(SuptBand(3), 0.9, Σ)

    N = 1000
    x = randn(N)
    y = x .+ 0.01.*randn(N)
    data = (y=y, x=x)
    m = lm(@formula(y ~ x), data)
    lb, ub = confint(PointwiseBand(), m)
    @test lb ≈ [0, 1] atol = 1e-2
    @test ub ≈ [0, 1] atol = 1e-2

    lb, ub = confint(SuptBand(), m)
    @test lb ≈ [0, 1] atol = 1e-2
    @test ub ≈ [0, 1] atol = 1e-2
    lb, ub = confint(SuptBand(1), m)
    @test lb ≈ [0, 1] atol = 1e-2
    @test ub ≈ [0, 1] atol = 1e-2
end
