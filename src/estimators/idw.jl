# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsModels.jl
# ------------------------------------------------------------------

# Local IDW estimator using local anisotropies
struct LocalIDWModel{E,L} <: GeoStatsModel
    exponent::E
    localaniso::L
end

"""
    LocalIDW(exponent, localaniso)

The inverse distance weighting model using local anisotropies.
The local parameters of the point to be estimated is used to
calculate the anisotropic distances. If exponent is not passed, 1.0 is
used by default
"""
LocalIDW(exponent, localaniso::LocalAnisotropy) = LocalIDWModel(exponent, localaniso)

LocalIDW(localaniso::LocalAnisotropy) = LocalIDWModel(1, localaniso)

function local_fit(model::LocalIDWModel, data; i, m)
    LA = model.localaniso
    state = IDWState(data)

    # get local anisotropy at estimation point
    Qi = qmat(rotation(LA, i), magnitude(LA, i))
    d = Mahalanobis(Symmetric(Qi))

    FittedIDW(IDW(model.exponent, d), state)
end
