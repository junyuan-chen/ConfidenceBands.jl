abstract type AbstractConfidenceBand end

abstract type PlugInConfidenceBand <: AbstractConfidenceBand end

criticalvalue(cb::PlugInConfidenceBand, level::Real, Σ::AbstractMatrix) =
    criticalvalue(cb, level, size(Σ, 1))

function confint(cb::PlugInConfidenceBand,
        θ::AbstractVector, Σ::AbstractMatrix; level::Real=0.9)
    se = sqrt.(diag(Σ))
    cv = criticalvalue(cb, level, Σ)
    return θ .- cv .* se, θ .+ cv .* se
end

confint(cb::PlugInConfidenceBand, m::StatisticalModel; level::Real=0.9) =
    confint(cb, coef(m), vcov(m), level=level)

struct PointwiseBand <: PlugInConfidenceBand end

criticalvalue(::PointwiseBand, level::Real, npara::Integer) =
    norminvccdf((1-level)/2)

struct SuptBand <: PlugInConfidenceBand
    nuncovered::Int
    ndraw::Int
    function SuptBand(nuncovered::Real=0; ndraw::Real=1_000_000)
        nuncovered < 0 && throw(ArgumentError("nuncovered cannot be negative"))
        ndraw > 0 || throw(ArgumentError("ndraw must be positive"))
        return new(Int(nuncovered), Int(ndraw))
    end
end

struct RandnPool
    v::Vector{Float64}
    RandnPool(N::Int) = new(randn(N))
end

function getmat!(p::RandnPool, N1::Int, N2::Int)
    N = N1 * N2
    N0 = length(p.v)
    if N > N0
        resize!(p.v, N)
        randn!(view(p.v, N0+1:N))
    end
    return reshape(view(p.v, 1:N), N1, N2)
end

const _globalrandnpool = RandnPool(40_000_000)

function criticalvalue(cb::SuptBand, level::Real, Σ::AbstractMatrix)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    nexc = cb.nuncovered
    K = size(Σ, 1)
    nexc < K || throw(ArgumentError(
        "value of cb.nuncovered ($(cb.nuncovered)) must be smaller than the number of parameters ($K)"))
    se = sqrt.(diag(Σ))
    ikept = se .> eps()
    corr = view(Σ, ikept, ikept) ./ view(se, ikept) ./ view(se, ikept)'
    F = eigen!(corr)
    A = F.vectors
    A .*= sqrt.(F.values')
    draws = getmat!(_globalrandnpool, size(corr, 1), cb.ndraw)
    ts = A * draws
    if iszero(nexc)
        return quantile!(view(maximum(abs, ts, dims=1), :), level)
    else
        f!(i) = partialsort!(view(ts,:,i), nexc+1, by=abs, rev=true)
        return quantile!(map(f!, axes(ts,2)), level)
    end
end

struct SidakBand <: PlugInConfidenceBand end 

criticalvalue(::SidakBand, level::Real, npara::Integer) =
    chisqinvcdf(1, level^(1/npara))^0.5


struct BonferroniBand <: PlugInConfidenceBand end

criticalvalue(::BonferroniBand, level::Real, npara::Integer) =
    norminvccdf((1-level)/(2*npara))

struct ProjectionBand{T} <: PlugInConfidenceBand
    npara::T
    function ProjectionBand(npara::Integer)
        N = Int(npara)
        N > 0 || throw(ArgumentError("npara must be positive"))
        return new{Int}(N)
    end
    ProjectionBand() = new{Nothing}(nothing)
end

criticalvalue(::ProjectionBand{Nothing}, level::Real, npara::Integer) =
    chisqinvcdf(npara, level)^0.5

# The npara argument is redundant here
criticalvalue(cb::ProjectionBand{Int}, level::Real, npara::Integer=cb.npara) =
    chisqinvcdf(cb.npara, level)^0.5
