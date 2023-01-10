# ConfidenceBands.jl

*Confidence bands for simultaneous statistical inference*

[![CI-stable][CI-stable-img]][CI-stable-url]
[![codecov][codecov-img]][codecov-url]
[![PkgEval][pkgeval-img]][pkgeval-url]

[CI-stable-img]: https://github.com/junyuan-chen/ConfidenceBands.jl/workflows/CI-stable/badge.svg
[CI-stable-url]: https://github.com/junyuan-chen/ConfidenceBands.jl/actions?query=workflow%3ACI-stable

[codecov-img]: https://codecov.io/gh/junyuan-chen/ConfidenceBands.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/junyuan-chen/ConfidenceBands.jl

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/C/ConfidenceBands.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/C/ConfidenceBands.html

[ConfidenceBands.jl](https://github.com/junyuan-chen/ConfidenceBands.jl)
is a lightweight Julia package for computing confidence bands
that are useful for simultaneous statistical inference.
In contrast to pointwise confidence intervals computed for each parameter separately,
a confidence band treats the entire vector of parameters as a single object
and is more suitable for comparisons involving multiple parameters.

## Example Usage

ConfidenceBands.jl extends the `confint` function for computing confidence bands.
Accepted arguments may vary depending on the type of confidence band.
Details may be found from docstrings in the help mode of Julia REPL.

#### Plug-In Confidence Bands

Computation of a plug-in confidence band is based on a critical value:

```julia
using ConfidenceBands
# Compute the critical value for Bonferroni bands with 90% confidence level
# when there are five parameters
criticalvalue(BonferroniBand(), 0.9, 5)
# A variance-covariance matrix Σ is required for sup-t bands
criticalvalue(SuptBand(), 0.9, Σ)
```

To obtain confidence bands:

```julia
# First obtain point estimates θ as a vector and variance-covariance matrix Σ
lb, ub = confint(SuptBand(), θ, Σ, level=0.95)
```

#### Bootstrap Confidence Bands

Some types of confidence bands are designed for
a valid bootstrap sample provided by users.
A bootstrap sample of point estimates needs to be collected in a matrix
with each column being a vector of point estimates from the same draw.
Currently, quantile-based and critical-value-based bootstrap implementation of
sup-t bands (`SuptQuantileBootBand` and `SuptCVBootBand`)
are implemented following Montiel Olea and Plagborg-Møller (2019):

```julia
lb, ub, pwlevel = confint(SuptQuantileBootBand(), draws)
lb, ub, cv = confint(SuptCVBootBand(), θ, draws)
```

The former additionally returns the confidence level
when the intervals from the confidence band
are viewed as pointwise confidence intervals.
The latter additionally returns the critical value.

## References

**Montiel Olea, José Luis and Mikkel Plagborg-Møller.** 2019.
"Simultaneous Confidence Bands: Theory, Implementation, and an Application to SVARs."
*Journal of Applied Econometrics* 34 (1): 1-17.
