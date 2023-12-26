# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LocalKriging(params ...)

LocalKriging estimation solver, where `var::Symbol` is the variable name
and `param` is a `NamedTuple` containing the parameters below:

## Parameters

* `variogram` - Reference variogram model
* `mean`      - Simple Kriging mean
* `method`    - LocalKriging method. :MovingWindows or :KernelConvolution
  (default to :MovingWindows)
* `localaniso`    - Local parameters of the domain
* `localanisohd`  - Local parameters of the samples. Only necessary for
  :KernelConvolution method. They are automatically passed via NN from
  `localaniso` if not informed.
"""

struct LocalKriging <: KrigingModel
  method::Symbol
	localaniso::LocalAnisotropy
  variogram::Variogram
  skmean::Union{Number,Nothing}  # for later: local mean array
end

LocalKriging(method::Symbol, localaniso::LocalAnisotropy, γ::Variogram; skmean=nothing) =
  LocalKriging(method, localaniso, γ, skmean)

abstract type LocalKriging <: GeoStatsModel end
