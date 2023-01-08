module ConfidenceBands

using LinearAlgebra: diag, eigen!
using Random: randn!, randn
using Roots: Brent, find_zero
using Statistics: quantile!, quantile, std
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
       SuptQuantileBootBand,
       SuptCVBootBand

include("plugin.jl")
include("bootstrap.jl")

end # module ConfidenceBands
