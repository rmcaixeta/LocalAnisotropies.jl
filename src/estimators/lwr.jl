# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsModels.jl
# ------------------------------------------------------------------

# Local LWR estimator using local anisotropies
struct LocalLWRModel{F,L} <: GeoStatsModel
  weightfun::F
  localaniso::L
end

"""
    LocalLWR(weightfun, localaniso)

The locally weighted regression model using local anisotropies.
The local parameters of the point to be estimated is used to
calculate the anisotropic distances.
"""
LocalLWR(weightfun, localaniso::LocalAnisotropy) = LocalLWRModel(weightfun, localaniso)

LocalLWR(localaniso::LocalAnisotropy) = LocalLWRModel(h -> exp(-3 * h ^ 2), localaniso)

function local_fit(model::LocalLWRModel, data; i, m)
  LA = model.localaniso

  # get local anisotropy at estimation point
  Qi = qmat(LA, i)
  d = Mahalanobis(Symmetric(Qi))

  # default LWR fit procedure
  model_ = GeoStatsModels.LWR(model.weightfun, d)
  X = GeoStatsModels.prealloc(model_, data)
  GeoStatsModels.setx!(model_, X, data)
  state = GeoStatsModels.LWRState(data, X)

  GeoStatsModels.FittedLWR(model_, state)
end

predictprob_(fitted::GeoStatsModels.FittedLWR, args...) = predictprob(fitted, args...)
