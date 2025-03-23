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
qmat(L::Union{LocalAnisotropy,Nothing}) = isnothing(L) ? nothing : [qmat(L, i) for i in 1:nvals(L)]

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

mw_estimator(model::MW_SKModel, localpar) = mw_estimator(model.mean, model.fun, localpar)
mw_estimator(model::MW_OKModel, localpar) = mw_estimator(nothing, model.fun, localpar)
mw_estimator(model::SimpleKriging, localpar) = mw_estimator(model.mean, model.fun, localpar)
mw_estimator(model::OrdinaryKriging, localpar) = mw_estimator(nothing, model.fun, localpar)

function kcfill!(Γ, γ::GeoStatsFunction, domain, localaniso)
  _, (_, n, k) = GeoStatsFunctions.matrixparams(γ, domain, domain)
  @inbounds for j in 1:n
    gⱼ = domain[j]
    Qj = localaniso[j]
    sⱼ = GeoStatsFunctions._sample(γ, gⱼ)
    # lower triangular entries
    for i in (j + 1):n
      gᵢ = domain[i]
      Qi = localaniso[i]
      sᵢ = GeoStatsFunctions._sample(γ, gᵢ)
      Γᵢⱼ = ustrip.(mean(kcvario(γ, pᵢ, pⱼ, Qi, Qj) for pᵢ in sᵢ, pⱼ in sⱼ))
      Γ[((i - 1) * k + 1):(i * k), ((j - 1) * k + 1):(j * k)] .= Γᵢⱼ
    end
    # diagonal entries
    Γᵢⱼ = ustrip.(mean(kcvario(γ, pⱼ, pⱼ, Qj, Qj) for pⱼ in sⱼ, pⱼ in sⱼ))
    Γ[((j - 1) * k + 1):(j * k), ((j - 1) * k + 1):(j * k)] .= Γᵢⱼ
  end

  # upper triangular entries
  @inbounds for j in 1:(n * k)
    for i in 1:(j - 1)
      Γ[i, j] = Γ[j, i]
    end
  end

  Γ
end

function kcfill!(F, f::GeoStatsFunction, domain₁, domain₂, Qx₀, hdla)
  _, (m, n, k) = GeoStatsFunctions.matrixparams(f, domain₁, domain₂)
  @inbounds for j in 1:n
    gⱼ = domain₂[j]
    sⱼ = GeoStatsFunctions._sample(f, gⱼ)
    for i in 1:m
      gᵢ = domain₁[i]
      sᵢ = GeoStatsFunctions._sample(f, gᵢ)
      Fᵢⱼ = ustrip.(mean(kcvario(f, pᵢ, pⱼ, hdla[i], Qx₀) for pᵢ in sᵢ, pⱼ in sⱼ))
      F[((i - 1) * k + 1):(i * k), ((j - 1) * k + 1):(j * k)] .= Fᵢⱼ
    end
  end
  F
end

function kccov(γ::GeoStatsFunction, xi, xj, Qi::AbstractMatrix, Qj::AbstractMatrix)
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

function kcvario(γ::GeoStatsFunction, xi, xj, Qi::AbstractMatrix, Qj::AbstractMatrix)
  C = kccov(γ, xi, xj, Qi, Qj)
  sill(γ) - C
end