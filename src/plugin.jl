"""
    AbstractConfidenceBand

Supertype for all confidence bands.
"""
abstract type AbstractConfidenceBand end

"""
    PlugInConfidenceBand <: AbstractConfidenceBand

Supertype for all plug-in confidence bands.
"""
abstract type PlugInConfidenceBand <: AbstractConfidenceBand end

"""
    criticalvalue(cb::PlugInConfidenceBand, level::Real, Σ::AbstractMatrix)

Return the critical value for `cb` with confidence level `level`
when the estimates have an estimated variance-covariance matrix `Σ`.
For some types of plug-in confidence bands,
providing the number of point estimates in place of `Σ` is sufficient.
"""
criticalvalue(cb::PlugInConfidenceBand, level::Real, Σ::AbstractMatrix) =
    criticalvalue(cb, level, size(Σ, 1))

"""
    confint(cb::PlugInConfidenceBand, θ::AbstractVector, Σ::AbstractMatrix; level::Real=0.9)

Compute the specified plug-in confidence band with confidence level `level`
using point estiamtes `θ` and variance-covariance matrix `Σ`.
"""
function confint(cb::PlugInConfidenceBand,
        θ::AbstractVector, Σ::AbstractMatrix; level::Real=0.9)
    se = sqrt.(diag(Σ))
    cv = criticalvalue(cb, level, Σ)
    return θ .- cv .* se, θ .+ cv .* se
end

"""
    confint(cb::PlugInConfidenceBand, m::StatisticalModel; level::Real=0.9)

Compute the specified plug-in confidence band with confidence level `level`
for the coefficients of model `m`.
"""
confint(cb::PlugInConfidenceBand, m::StatisticalModel; level::Real=0.9) =
    confint(cb, coef(m), vcov(m), level=level)

"""
    PointwiseBand <: PlugInConfidenceBand

Pointwise confidence intervals with critical values based on a normal distribution.
"""
struct PointwiseBand <: PlugInConfidenceBand end

criticalvalue(::PointwiseBand, level::Real, npara::Integer) =
    norminvccdf((1-level)/2)

const _MOPMreference = """
    # References
    - Montiel Olea, José Luis and Mikkel Plagborg-Møller. 2019.
      "Simultaneous Confidence Bands: Theory, Implementation, and an Application to SVARs."
      Journal of Applied Econometrics 34 (1): 1-17."""

"""
    SuptBand <: PlugInConfidenceBand

Plug-in sup-t confidence band.
Implementation follows Montiel Olea and Plagborg-Møller (2019) Algorithm 1
and may allow generalized error rate control.

Critical values computed for `SuptBand` are based on
random draws from a normal distribution.
Since the random numbers are drawn only once
and stored in an unexported global object `_globalrandnpool`,
results from the same Julia session remain unchanged if executed multiple times.
However, results obtained across different sessions are not identical
because the random numbers generated vary.
See Julia manual section on [`Random`](https://docs.julialang.org/en/v1/stdlib/Random/)
for reproducibility of random numbers.

$_MOPMreference
"""
struct SuptBand <: PlugInConfidenceBand
    nuncovered::Int
    ndraw::Int
    function SuptBand(nuncovered::Real, ndraw::Real)
        nuncovered < 0 && throw(ArgumentError("nuncovered cannot be negative"))
        ndraw > 0 || throw(ArgumentError("ndraw must be positive"))
        return new(Int(nuncovered), Int(ndraw))
    end
end

"""
    SuptBand(nuncovered::Real=0; ndraw::Real=1_000_000)

Return a `SuptBand` instance that requires `ndraw` random numbers for computation.
Results tend to be more accurate with a larger value of `ndraw`.

A positive value of `nuncovered` allows generalized error rate control
described by Montiel Olea and Plagborg-Møller (2019).
Specifically, at most `nuncovered` number of point estimates
are allowed to be not covered by the confidence band
when considering the coverage.

$_MOPMreference
"""
SuptBand(nuncovered::Real=0; ndraw::Real=1_000_000) = SuptBand(nuncovered, ndraw)

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
        "value of cb.nuncovered ($nexc) must be smaller than the number of parameters ($K)"))
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

"""
    SidakBand <: PlugInConfidenceBand

Šidák band with exact asymptotic simultaneous coverage
only for point estimators that are uncorrelated elementwise.
"""
struct SidakBand <: PlugInConfidenceBand end

criticalvalue(::SidakBand, level::Real, npara::Integer) =
    chisqinvcdf(1, level^(1/npara))^0.5

"""
    BonferroniBand <: PlugInConfidenceBand

Confidence band with pointwise significance level
adjusted by a Bonferroni correction for multiple hypotheses.
"""
struct BonferroniBand <: PlugInConfidenceBand end

criticalvalue(::BonferroniBand, level::Real, npara::Integer) =
    norminvccdf((1-level)/(2*npara))

"""
    ProjectionBand{T} <: PlugInConfidenceBand

The smallest (rectangular) confidence band that contains
the Wald confidence ellipsoid for parameters with a given dimension.
"""
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
