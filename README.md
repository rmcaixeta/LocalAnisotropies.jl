# LocalAnisotropies

[![Build Status][build-img]][build-url] [![Coverage][codecov-img]][codecov-url]

`LocalAnisotropies.jl` is a implementation of geostatistics methods to deal with the non-stationarity of second order moments (aka locally varying anisotropies). It is developed to be used as an extension of [`GeoStats.jl`](https://github.com/JuliaEarth/GeoStats.jl).

**Warning**: This package is still under (slow) development. Some of the implementations were not fully validated and may give inconsistent results.

## Introduction

This package offer some solutions to extract local parameters from a reference input, model non-stationary covariance models, and adapt estimation methods to be used with them. There are some extra tools, like local parameters interpolation and conversion to/between different rotation conventions. A list of current implementations:

- <u>Local parameters extraction methods</u>:
  - Gradients
  - SVD / PCA <p>
- <u>Nonstationary spatial methods</u>:
  - Moving windows
  - Kernel convolution
  - Spatial deformation <p>
- <u>Estimation methods adapted to them</u>:
  - Kriging
  - IDW

## Installation

First, it is necessary to install Julia. Installation instructions for Windows, Linux and macOS are available [here](https://julialang.org/downloads/platform/).

`LocalAnisotropies.jl` is not released yet. To use it: open the Julia REPL and then run the following command.

```julia
using Pkg; Pkg.develop(url="https://github.com/rmcaixeta/LocalAnisotropies.jl"); Pkg.add("GeoStats")
```

## Usage example

```julia
using LocalAnisotropies
using GeoStats
using Plots

# reference scenario for tests
D = georef((P=[sin(i)+j for i in 1:20, j in 1:20],))
S = sample(D, 80, replace=false)
G = CartesianGrid(20,20)

# Estimation problem
P = EstimationProblem(S, G, :P)
γ = NuggetEffect(0.1) + 0.9*ExponentialVariogram(range=60.0)

searcher = KNearestSearch(G, 10)

plot(D)
```

```julia
# get local parameters
lpars = localparameters(Gradients(), D, :P, 3)
plot(lpars)
```

```julia
# rescale magnitude and average 10 nearest local parameters
lpars = rescale_magnitude(lpars, (0.2,1.0))
lpars = smoothpars(lpars, searcher)
plot(lpars)
```

```julia
# LocalKriging (MW)
MW = LocalKriging(:P => (variogram=(:X=>γ), localpars=lpars, method=:MovingWindows))
s1 = solve(P, MW)
plot(s1)
```

```julia
# LocalKriging (KC)
KC = LocalKriging(:P => (variogram=(:X=>γ), localpars=lpars, method=:KernelConvolution))
s2 = solve(P, KC)
plot(s2)
```

```julia
# Spatial deformation: anisotropic distances
Sd1, Dd1 = deformspace(S, G, lpars, LocalAnisotropy(), anchors=1500)
Pd1 = EstimationProblem(Sd1, Dd1, :P)
s3 = solve(Pd1, Kriging(:P => (variogram=γ,)))
plot(s3)
```

```julia
# Spatial deformation: anisotropic variogram distances
Sd2, Dd2 = deformspace(S, G, lpars, LocalVariogram(), γ, anchors=1500)
Pd2 = EstimationProblem(Sd2, Dd2, :P)
s4 = solve(Pd2, Kriging(:P => (variogram=γ,)))
plot(s4)
```

```julia
# Spatial deformation: geodesic anisotropic distances
LDa = addgraph(S, G, lpars, LocalAnisotropy(), searcher)
Sd3, Dd3 = deformspace(LDa, GraphDistance(), anchors=1500)
Pd3 = EstimationProblem(Sd3, Dd3, :P)
s5 = solve(Pd3, Kriging(:P => (variogram=γ,)))
plot(s5)
```

```julia
# Spatial deformation: geodesic anisotropic variogram distances
LDv = addgraph(S, G, lpars, LocalVariogram(), γ, searcher)
Sd4, Dd4 = deformspace(LDv, GraphDistance(), anchors=1500)
Pd4 = EstimationProblem(Sd4, Dd4, :P)
s6 = solve(Pd4, Kriging(:P => (variogram=γ,)))
plot(s6)
```

Some extra tools to work with local parameters:

```julia
# import external local parameters in GSLIB convention
dummy = georef((az=1:10, r1=1:10, r2=1:10), PointSet(rand(2,10)))
pars  = LocalParameters(dummy, [:az], [:r1,:r2], :GSLIB)

# interpolate local parameters into a coarser grid
G_ = CartesianGrid((10,10),(0.5,0.5),(2.0,2.0))
lpars_ = IDWpars(lpars, searcher, G_, power=2)

# convert between different rotation conventions
angs1 = convertangles([30,30,30], :GSLIB, :Datamine)
angs2 = convertangles.(pars.rotation, :GSLIB)
```

## Documentation

Not available yet

## References

##### Introduction to local parameters extraction methods:
[Lillah & Boisvert (2015)](https://doi.org/10.1016/j.cageo.2015.05.015) Inference of locally varying anisotropy fields from diverse data sources

##### Introduction to nonstationary spatial methods:
[Sampson (2010)](https://doi.org/10.1201/9781420072884-13) Constructions for Nonstationary Spatial Processes

##### Moving windows:
[Haas (1990)](https://doi.org/10.1016/0960-1686(90)90508-K) Kriging and automated variogram modeling within a moving window <br>
[Stroet & Snepvangers (2005)](https://doi.org/10.1007/s11004-005-7310-y) Mapping curvilinear structures with local anisotropy kriging

##### Kernel convolution:
[Higdon (1998)](https://doi.org/10.1023/A:1009666805688) A process-convolution approach to modelling temperatures in the North Atlantic Ocean <br>
[Fouedjio et al. (2016)](https://doi.org/10.1016/j.spasta.2016.01.002) A generalized convolution model and estimation for non-stationary random functions

##### Spatial deformation:
[Sampson & Guttorp (1992)](https://doi.org/10.1080/01621459.1992.10475181) Nonparametric estimation of nonstationary spatial covariance structure <br>
[Boisvert (2010)](https://era.library.ualberta.ca/items/5acca59f-6e97-414d-ad13-34c8f97ce223) Geostatistics with locally varying anisotropy



[build-img]: https://img.shields.io/github/workflow/status/rmcaixeta/LocalAnisotropies.jl/CI?style=flat-square
[build-url]: https://github.com/rmcaixeta/LocalAnisotropies.jl/actions

[codecov-img]: https://codecov.io/gh/rmcaixeta/LocalAnisotropies.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/rmcaixeta/LocalAnisotropies.jl
