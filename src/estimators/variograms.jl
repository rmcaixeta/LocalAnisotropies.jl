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
  Q
end

function mwvario(estimator, localpar)
  # get local Mahalanobis
  Q = qmat(localpar[1],localpar[2])
  local_d = Mahalanobis(Q)

  # get reference pars. and apply local anisotropy to given structures
  p = structures(estimator.γ)
  γs = map(p[3]) do γ
    γ = @set γ.distance = local_d
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
  # X = neighbors coords
  # need to consider LHS[i,j] = sill(γ) - LHS[i,j] to convert vario to covario
  n = nelements(X)
  @inbounds for j=1:n
    xj = centroid(X, j)
    Qj = localaniso[j]
    for i=j+1:n
      xi = centroid(X, i)
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
  local_d = Mahalanobis(Qij)

  # get reference pars. and apply local anisotropy to given structures
  p = structures(γ)
  γs = map(p[3]) do γ
    γ = @set γ.distance = local_d
  end
  γl = NuggetEffect(p[1]) + sum(c*γ for (c, γ) in zip(p[2], γs))

  (det(Qi)^0.25)*(det(Qj)^0.25)*(det(Qij)^-0.5)*(sill(γ)-γ(xi,xj))
end
