# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# Variogram functions using local anisotropies
function qmat(q::Quaternion, m::AbstractVector)
    # Anisotropy matrix for Mahalanobis distance
    N = length(m)
    P = quat_to_dcm(q)[SOneTo(N), SOneTo(N)]
    Λ = Diagonal(SVector{N}(one(eltype(m)) ./ m .^ 2))
    Q = P' * Λ * P
    # need to explicitly set Symmetric(Q), but then det(Q) will fail with StaticArrays
    # while not fixed, setting Symmetric only when calling Mahalanobis
    # https://github.com/JuliaArrays/StaticArrays.jl/issues/955
end

qmat(L::LocalAnisotropy, i::Int) = qmat(localpair(L, i)...)
qmat(L::LocalGeoData, i::Int) = qmat(localpair(L, i)...)
qmat(L::Union{LocalAnisotropy,Nothing}) =
    isnothing(L) ? nothing : [qmat(L, i) for i = 1:nvals(L)]

function mw_estimator(μ, γ, localpar)
    # get local Mahalanobis matrix
    Q = qmat(localpar...)
    # get reference pars. and apply local anisotropy to given structures
    p = structures(γ)
    γs = map(p[3]) do γ
        Qs = ustrip.(Q ./ collect(radii(γ.ball))' .^ 2)
        γ = @set γ.ball = MetricBall(1.0, Mahalanobis(Symmetric(Qs)))
    end
    γl = NuggetEffect(p[1]) + sum(c * γ for (c, γ) in zip(p[2], γs))

    # return local estimator
    isnothing(μ) ? OrdinaryKriging(γl) : SimpleKriging(γl, μ)
end

mw_estimator(model::MW_SKModel, localpar) = mw_estimator(model.μ, model.γ, localpar)
mw_estimator(model::MW_OKModel, localpar) = mw_estimator(nothing, model.γ, localpar)
mw_estimator(model::SimpleKriging, localpar) = mw_estimator(model.μ, model.γ, localpar)
mw_estimator(model::OrdinaryKriging, localpar) = mw_estimator(nothing, model.γ, localpar)

function kcfill!(Γ, γ::Variogram, X, localaniso)
    n = length(X)
    @inbounds for j = 1:n
        xj = X[j]
        Qj = localaniso[j]
        for i = j+1:n
            xi = X[i]
            Qi = localaniso[i]
            Γ[i, j] = kccov(γ, xi, xj, Qi, Qj)
        end
        Γ[j, j] = sill(γ)
        for i = 1:j-1
            Γ[i, j] = Γ[j, i]
        end
    end
end

function kccov(γ::Variogram, xi, xj, Qi::AbstractMatrix, Qj::AbstractMatrix)
    Qij = (Qi + Qj) / 2

    ## modify variogram with average anisotropy matrix
    p = structures(γ)
    γs = map(p[3]) do γx
        Qs = ustrip.(Qij ./ collect(radii(γx.ball))' .^ 2)
        γx = @set γx.ball = MetricBall(1.0, Mahalanobis(Symmetric(Qs)))
    end
    γl = NuggetEffect(p[1]) + sum(c * γx for (c, γx) in zip(p[2], γs))

    Cij = (sill(γl) - γl(xi, xj))
    (det(Qi)^0.25) * (det(Qj)^0.25) * (det(Qij)^-0.5) * Cij
end
