# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LocalKriging(var₁=>param₁, var₂=>param₂, ...)

LocalKriging estimation solver, where `var::Symbol` is the variable name
and `param` is a `NamedTuple` containing the parameters below:

## Parameters

* `variogram` - Reference direction and variogram model informed as `axis => γ`.
  The model is fixed in that direction and rescaled in the others according to
  `LocalParameters` information. e.g. :X => ExponentialVariogram()
* `mean`      - Simple Kriging mean
* `method`    - LocalKriging method. :MovingWindows or :KernelConvolution
  (default to :MovingWindows)
* `localpars`    - Local parameters of the domain
* `localparshd`  - Local parameters of the samples. Only necessary for
  :KernelConvolution method. They are automatically passed via NN from
  `localpars` if not informed.
* `minneighbors` - Minimum number of neighbors (default to 1)
* `maxneighbors` - Maximum number of neighbors (default to 20)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)
"""
@estimsolver LocalKriging begin
  @param variogram = (:X => ExponentialVariogram())
  @param mean = nothing #0.0 or [1.0,0.8,....]
  @param method = :MovingWindows
  @param localpars = nothing
  @param localparshd = nothing
  @param minneighbors = 1
  @param maxneighbors = 20
  @param neighborhood = nothing
  @param distance = Euclidean()
end
