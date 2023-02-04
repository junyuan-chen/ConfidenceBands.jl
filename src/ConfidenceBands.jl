module ConfidenceBands

using LinearAlgebra: diag, eigen!
using Random: randn!, randn
using Roots: Brent, find_zero
using Statistics: quantile!, quantile, _quantilesort!, _quantile, std
using StatsAPI: StatisticalModel, coef, vcov, stderror
using StatsFuns: norminvccdf, chisqinvcdf

import StatsAPI: confint

export confint

export AbstractConfidenceBand,
       PlugInConfidenceBand,
       PointwiseBand,
       SuptBand,
       SidakBand,
       BonferroniBand,
       ProjectionBand,
       criticalvalue,

       BootstrapConfidenceBand,
       PointwiseQuantileBootBand,
       PointwiseCVBootBand,
       SuptQuantileBootBand,
       SuptCVBootBand,
       simulcoverage,
       pwcoverage

include("plugin.jl")
include("bootstrap.jl")

end # module ConfidenceBands
