"""
    BootstrapConfidenceBand <: AbstractConfidenceBand

Supertype for all confidence bands that require
a bootstrap sample of parameters of interest.
"""
abstract type BootstrapConfidenceBand <: AbstractConfidenceBand end

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

function quantilebands!(lb::AbstractVector, ub::AbstractVector, draws::AbstractMatrix,
        pointwiselevel::Real; sorted::Bool=false)
    α = (1 - pointwiselevel) / 2
    map!(i->quantile(view(draws,i,:), α, sorted=sorted), lb, axes(draws,1))
    map!(i->quantile(view(draws,i,:), 1-α, sorted=sorted), ub, axes(draws,1))
    return lb, ub
end

function bandcoverage(lb::AbstractVector, ub::AbstractVector, draws::AbstractMatrix,
        nuncovered::Integer)
    K = sum(lb .<= draws .<= ub, dims=1)
    M, N = size(draws)
    return sum(K .>= (M-nuncovered)) / N
end

_default_root_solver = Brent()

"""
    confint(cb::SuptQuantileBootBand, draws::AbstractMatrix; level::Real=0.9, kwargs...)

Compute a sup-t confidence band with quantile-based bootstrap implementation
based on Montiel Olea and Plagborg-Møller (2019) Algorithm 2.
The bootstrap sample of point estimates `draws` is required to be a matrix
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
    sorted = sort(draws, dims=2)
    function f(x)
        quantilebands!(lb, ub, sorted, x, sorted=true)
        return bandcoverage(lb, ub, draws, nexc) - level
    end
    pwlevel = find_zero(f, x0, solver; kwargs...)
    # Last call by solver is not always at the solution
    quantilebands!(lb, ub, sorted, pwlevel, sorted=true)
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
The bootstrap sample of point estimates `draws` is required to be a matrix
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
    cv = quantile!(view(mts, :), 1-level)
    return θ0 .- cv .* σs, θ0 .+ cv .* σs, cv
end
