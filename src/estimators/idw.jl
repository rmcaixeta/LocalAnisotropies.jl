# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsModels.jl
# ------------------------------------------------------------------

# Local IDW estimator using local anisotropies

struct LocalIDWModel <: GeoStatsModel
    exponent::Float64
    localaniso::LocalAnisotropy
    hdlocalaniso::Union{AbstractVector,Nothing}
end

"""
    LocalIDW(params ...)
"""
LocalIDW(exponent::Float64, localaniso::LocalAnisotropy, ; hdlocalaniso = nothing) =
    LocalIDWModel(exponent, localaniso, hdlocalaniso)

struct LocalFittedIDW{M<:LocalIDWModel,S<:IDWState}
    model::M
    state::S
end

function local_fit(model::LocalIDWModel, data)
    # record state
    state = IDWState(data)

    # return fitted model
    LocalFittedIDW(model, state)
end


predict(fitted::LocalFittedIDW, var, uₒ) = idw(fitted, weights(fitted, uₒ), var)

predictprob(fitted::LocalFittedIDW, var, uₒ) = Dirac(predict(fitted, var, uₒ))

function idw(fitted::LocalFittedIDW, weights, var)
    d = fitted.state.data
    c = Tables.columns(values(d))
    z = Tables.getcolumn(c, var)
    w = weights
    Σw = sum(w)

    λ(i) = w[i] / Σw

    if isinf(Σw) # some distance is zero?
        z[findfirst(isinf, w)]
    else
        sum(i -> λ(i) * z[i], eachindex(z))
    end
end

function weights(fitted::LocalFittedIDW, uₒ)
    e = fitted.model.exponent
    δ = fitted.model.distance
    d = fitted.state.data
    Ω = domain(d)

    # adjust CRS of uₒ
    uₒ′ = uₒ |> Proj(crs(Ω))

    pₒ = centroid(uₒ′)
    p(i) = centroid(Ω, i)

    λ(i) = 1 / evaluate(δ, pₒ, p(i))^e

    map(λ, 1:nelements(Ω))
end
