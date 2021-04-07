# Variogram functions using local parameters

function qmat(q,m)
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

function kcfill!(Γ, γ::Variogram, X, localpars)
  # X = neighbors coords
  # need to consider LHS[i,j] = sill(γ) - LHS[i,j] to convert vario to covario
  n = nelements(X)
  @inbounds for j=1:n
    xj = centroid(X, j)
    Qj = localpars[j]
    for i=j+1:n
      xi = centroid(X, i)
      Qi = localpars[i]
      Γ[i,j] = kccov(γ, xi, xj, Qi, Qj)
    end
    Γ[j,j] = sill(γ) # kccov(γ, xj, xj, Qj, Qj)
    for i=1:j-1
      Γ[i,j] = Γ[j,i] # leverage the symmetry
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

function setref_axis(localpars::LocalParameters, ax::Symbol)
  ix = Dict(:X=>1,:Y=>2,:Z=>3)
  m = localpars.magnitude
  ref = m[ix[ax],:]
  for i in size(m, 1)
    m[i,:] ./= ref
  end
  LocalParameters(localpars.rotation,m)
end
