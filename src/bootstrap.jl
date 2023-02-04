"""
    BootstrapConfidenceBand <: AbstractConfidenceBand

Supertype for all confidence bands that require
bootstrap draws of parameters of interest.
"""
abstract type BootstrapConfidenceBand <: AbstractConfidenceBand end

"""
    simulcoverage(lb::AbstractVector, ub::AbstractVector, draws::AbstractMatrix, nuncovered=0)

Compute the share of point estimates in bootstrap `draws`
that are simultaneously covered by the confidence band
with lower bound `lb` and upper bound `ub`,
except for at most `nuncovered` parameters.
See also [`pwcoverage`](@ref).
"""
function simulcoverage(lb::AbstractVector, ub::AbstractVector, draws::AbstractMatrix,
        nuncovered::Integer=0)
    K = sum(lb .<= draws .<= ub, dims=1)
    M, N = size(draws)
    return sum(K .>= (M-nuncovered)) / N
end

"""
    pwcoverage(lb::AbstractVector, ub::AbstractVector, draws::AbstractMatrix)

Compute the share of point estimates in bootstrap `draws`
that are covered by the corresponding confidence intervals
with lower bounds `lb` and upper bounds `ub` for each parameter separately.
See also [`simulcoverage`](@ref).
"""
pwcoverage(lb::AbstractVector, ub::AbstractVector, draws::AbstractMatrix) =
    vec(sum(lb .<= draws .<= ub, dims=2)) ./ size(draws, 2)

"""
    PointwiseQuantileBootBand <: BootstrapConfidenceBand

Efron's pointwise equal-tailed percentile bootstrap confidence intervals.
"""
struct PointwiseQuantileBootBand <: BootstrapConfidenceBand end

function pwquantilebands!(lb::AbstractVector, ub::AbstractVector, draws::AbstractMatrix,
        pointwiselevel::Real; sorted::Bool=false)
    α = (1 - pointwiselevel) / 2
    # Modified from quantile!
    sorteddraws = sorted ? draws : collect(draws)
    @inbounds for i in axes(draws,1)
        v = view(sorteddraws, i, :)
        _quantilesort!(v, sorted, α, 1-α)
        lb[i] = _quantile(v, α)
        ub[i] = _quantile(v, 1-α)
    end
    return lb, ub
end

function confint(::PointwiseQuantileBootBand, draws::AbstractMatrix; level::Real=0.9)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    K = size(draws, 1)
    lb = Vector{Float64}(undef, K)
    ub = Vector{Float64}(undef, K)
    sorteddraws = sort(draws, dims=2)
    pwquantilebands!(lb, ub, sorteddraws, level, sorted=true)
    return lb, ub
end

"""
    PointwiseCVBootBand <: BootstrapConfidenceBand

Pointwise bootstrap confidence intervals based on the quantile of t-statistics.
"""
struct PointwiseCVBootBand <: BootstrapConfidenceBand end

function confint(::PointwiseCVBootBand, θ0::AbstractVector, draws::AbstractMatrix;
        level::Real=0.9)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    σs = vec(std(draws, dims=2))
    ts = abs.(draws .- θ0) ./ σs
    cvs = map(i->quantile!(view(ts, i, :), level), axes(ts, 1))
    return θ0 .- cvs .* σs, θ0 .+ cvs .* σs
end

"""
    SuptQuantileBootBand <: BootstrapConfidenceBand

Quantile-based bootstrap implementation of sup-t confidence band.
Implementation follows Montiel Olea and Plagborg-Møller (2019) Algorithm 2
and may allow generalized error rate control.

$_MOPMreference
"""
struct SuptQuantileBootBand <: BootstrapConfidenceBand
    nuncovered::Int
    function SuptQuantileBootBand(nuncovered::Real=0)
        nuncovered < 0 && throw(ArgumentError("nuncovered cannot be negative"))
        return new(Int(nuncovered))
    end
end

_default_root_solver = Brent()

"""
    confint(cb::SuptQuantileBootBand, draws::AbstractMatrix; level::Real=0.9, kwargs...)

Compute a sup-t confidence band with quantile-based bootstrap implementation
based on Montiel Olea and Plagborg-Møller (2019) Algorithm 2.
The bootstrap `draws` of point estimates need to be in a matrix
with each column being a vector of point estimates from the same draw.
In addition to the lower and upper bounds,
the pointwise confidence level (when the intervals from the confidence band
are viewed as pointwise confidence intervals)
is returned as the third object.

The procedure involves solving a root-finding problem for
seeking the band with the specified confidence level.
This is accomplished with the `find_zero` function from
[`Roots.jl`](https://github.com/JuliaMath/Roots.jl).
The default bracketing interval (or starting point) used to solve this problem
can be overriden by specifying the keyword argument `x0`.
A solver from `Roots.jl` can be specified with keyword argument `solver`.
Any additional keyword argument will be passed to `find_zero`.

$_MOPMreference
"""
function confint(cb::SuptQuantileBootBand, draws::AbstractMatrix;
        level::Real=0.9, solver=_default_root_solver,
        x0=(cb.nuncovered==0 ? level : level/size(draws,1), 1-1e-3), kwargs...)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    K = size(draws, 1)
    nexc = cb.nuncovered
    nexc < K || throw(ArgumentError(
        "value of cb.nuncovered ($nexc) must be smaller than the number of parameters ($K)"))
    lb = Vector{Float64}(undef, K)
    ub = Vector{Float64}(undef, K)
    sorteddraws = sort(draws, dims=2)
    function f(x)
        pwquantilebands!(lb, ub, sorteddraws, x, sorted=true)
        return simulcoverage(lb, ub, draws, nexc) - level
    end
    pwlevel = find_zero(f, x0, solver; kwargs...)
    # Last call by solver is not always at the solution
    pwquantilebands!(lb, ub, sorteddraws, pwlevel, sorted=true)
    return lb, ub, pwlevel
end

"""
    SuptCVBootBand <: BootstrapConfidenceBand

Critical-value-based bootstrap implementation of sup-t confidence band.
Implementation follows Montiel Olea and Plagborg-Møller (2019)
Algorithm 3 in appendix and may allow generalized error rate control.

$_MOPMreference
"""
struct SuptCVBootBand <: BootstrapConfidenceBand
    nuncovered::Int
    function SuptCVBootBand(nuncovered::Real=0)
        nuncovered < 0 && throw(ArgumentError("nuncovered cannot be negative"))
        return new(Int(nuncovered))
    end
end

"""
    confint(cb::SuptCVBootBand, θ0::AbstractVector, draws::AbstractMatrix; level::Real=0.9)

Compute a sup-t confidence band with critical-value-based bootstrap implementation
based on Montiel Olea and Plagborg-Møller (2019) Algorithm 3 in appendix.
`θ0` is a vector of point estimates to be used as the middle points of the band.
The bootstrap `draws` of point estimates need to be in a matrix
with each column being a vector of point estimates from the same draw.
In addition to the lower and upper bounds,
the critical value is returned as the third object.

$_MOPMreference
"""
function confint(cb::SuptCVBootBand, θ0::AbstractVector, draws::AbstractMatrix;
        level::Real=0.9)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    σs = view(std(draws, dims=2), :)
    nexc = cb.nuncovered
    ts = abs.(draws .- θ0) ./ σs
    if iszero(nexc)
        mts = maximum(ts, dims=1)
    else
        K = size(draws, 1)
        nexc < K || throw(ArgumentError(
            "value of cb.nuncovered ($nexc) must be smaller than the number of parameters ($K)"))
        mts = map(i->partialsort(view(ts,:,i), nexc+1, rev=true), axes(ts, 2))
    end
    cv = quantile!(view(mts, :), level)
    return θ0 .- cv .* σs, θ0 .+ cv .* σs, cv
end
