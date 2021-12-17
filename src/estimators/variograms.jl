# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# Variogram functions using local anisotropies

function qmat(q,m)
  # Anisotropy matrix for Mahalanobis distance
  N = length(m)
  P = quat_to_dcm(q)[SOneTo(N),SOneTo(N)]
  Λ = Diagonal(SVector{N}(one(eltype(m))./m.^2))
  Q = P'*Λ*P
  # need to explicitly set Symmetric(Q), but then det(Q) will fail with StaticArrays
  # while not fixed, setting Symmetric only when calling Mahalanobis
  # https://github.com/JuliaArrays/StaticArrays.jl/issues/955
end

function mwvario(estimator, localpar)
  # get local Mahalanobis matrix
  Q = qmat(localpar[1],localpar[2])
  # get reference pars. and apply local anisotropy to given structures
  p = structures(estimator.γ)
  γs = map(p[3]) do γ
    Qs = Q ./ radii(γ.ball)' .^ 2
    γ = @set γ.ball = MetricBall(1.0, Mahalanobis(Symmetric(Qs)))
  end
  γl = NuggetEffect(p[1]) + sum(c*γ for (c, γ) in zip(p[2], γs))

  # return local estimator
  if typeof(estimator) <: SimpleKriging
    return SimpleKriging(γl, estimator.mean)
  elseif typeof(estimator) <: OrdinaryKriging
    return OrdinaryKriging(γl)
  end

end

function kcfill!(Γ, γ::Variogram, X, localaniso)
  n = length(X)
  @inbounds for j in 1:n
    xj = X[j]
    Qj = localaniso[j]
    for i=j+1:n
      xi = X[i]
      Qi = localaniso[i]
      Γ[i,j] = kccov(γ, xi, xj, Qi, Qj)
    end
    Γ[j,j] = sill(γ)
    for i=1:j-1
      Γ[i,j] = Γ[j,i]
    end
  end
end

function kccov(γ::Variogram, xi, xj, Qi::AbstractMatrix, Qj::AbstractMatrix)
  Qij = (Qi+Qj)/2

  ## modify variogram with average anisotropy matrix
  p = structures(γ)
  γs = map(p[3]) do γx
    Qs = Qij ./ radii(γx.ball)' .^ 2
    γx = @set γx.ball = MetricBall(1.0, Mahalanobis(Symmetric(Qs)))
  end
  γl = NuggetEffect(p[1]) + sum(c*γx for (c, γx) in zip(p[2], γs))

  Cij = (sill(γl)-γl(xi,xj))
  (det(Qi)^0.25)*(det(Qj)^0.25)*(det(Qij)^-0.5)*Cij
end
