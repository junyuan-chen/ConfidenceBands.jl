abstract type BootstrapConfidenceBand <: AbstractConfidenceBand end

struct SuptQuantileBootBand <: BootstrapConfidenceBand
    nuncovered::Int
    function SuptQuantileBootBand(nuncovered::Real=0)
        nuncovered < 0 && throw(ArgumentError("nuncovered cannot be negative"))
        return new(Int(nuncovered))
    end
end

function quantilebands!(lb::AbstractVector, ub::AbstractVector, sample::AbstractMatrix,
        pointwiselevel::Real; sorted::Bool=false)
    α = (1 - pointwiselevel) / 2
    map!(i->quantile(view(sample,i,:), α, sorted=sorted), lb, axes(sample,1))
    map!(i->quantile(view(sample,i,:), 1-α, sorted=sorted), ub, axes(sample,1))
    return lb, ub
end

function bandcoverage(lb::AbstractVector, ub::AbstractVector, sample::AbstractMatrix,
        nuncovered::Integer)
    K = sum(lb .<= sample .<= ub, dims=1)
    M, N = size(sample)
    return sum(K .>= (M-nuncovered)) / N
end

_default_root_solver = Brent()

function confint(cb::SuptQuantileBootBand, sample::AbstractMatrix;
        level::Real=0.9, solver=_default_root_solver,
        x0=(level, 1-1e-3), kwargs...)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    K = size(sample, 1)
    lb = Vector{Float64}(undef, K)
    ub = Vector{Float64}(undef, K)
    sorted = sort(sample, dims=2)
    function f(x)
        quantilebands!(lb, ub, sorted, x, sorted=true)
        return bandcoverage(lb, ub, sample, cb.nuncovered) - level
    end
    pwlevel = find_zero(f, x0, solver; kwargs...)
    # Last call by solver is not always at the solution
    quantilebands!(lb, ub, sorted, pwlevel, sorted=true)
    return lb, ub, pwlevel
end

struct SuptCVBootBand <: BootstrapConfidenceBand
    nuncovered::Int
    function SuptCVBootBand(nuncovered::Real=0)
        nuncovered < 0 && throw(ArgumentError("nuncovered cannot be negative"))
        return new(Int(nuncovered))
    end
end

function confint(cb::SuptCVBootBand, θ0::AbstractVector, sample::AbstractMatrix;
        level::Real=0.9)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    σs = view(std(sample, dims=2), :)
    nexc = cb.nuncovered
    ts = abs.(sample .- θ0) ./ σs
    if iszero(nexc)
        mts = maximum(ts, dims=1)
    else
        mts = map(i->partialsort(view(sample,:,i), nexc+1, rev=true), axes(sample, 2))
    end
    cv = quantile!(view(mts, :), 1-level)
    return θ0 .- cv .* σs, θ0 .+ cv .* σs, cv
end
