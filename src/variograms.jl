
"""
 Variogram function using local pars
 Local pars can be different between structures
 e.g. γns = γ1(pars,local_pars) + γ2(pars) # only vary in the first structure
"""

function get_pars(γ::Variogram)
  #### do for multistructure later
  kw = Dict(
    :range => γ.range,
    :distance => γ.distance,
    :nugget => γ.nugget,
    :sill => γ.sill
  )

  st1 = Dict(
    :type => getfield(Main, nameof(typeof(γ))),
    :kw => kw
  )

  [st1]
end

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
  ref_pars = get_pars(estimator.γ)
  ### loop multi structure after updated struct
  kwargs = ref_pars[1][:kw]
  kwargs[:distance] = local_d
  γl = ref_pars[1][:type](;kwargs...)

  # return local estimator
  if typeof(estimator) <: SimpleKriging
    return SimpleKriging(γl, estimator.mean)
  elseif typeof(estimator) <: OrdinaryKriging
    return OrdinaryKriging(γl)
  end

end

function kcfill!(Γ, γ::Variogram, X::AbstractMatrix, localpars)
  # X = neighbors coords
  # need to consider LHS[i,j] = sill(γ) - LHS[i,j] to convert vario to covario
  m, n = size(X)
  @inbounds for j=1:n
    xj = view(X, :, j)
    Qj = localpars[j]
    for i=j+1:n
      xi = view(X, :, i)
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
  ref_pars = get_pars(γ)
  ### loop multi structure after updated struct
  kwargs = ref_pars[1][:kw]
  kwargs[:distance] = local_d
  γl = ref_pars[1][:type](;kwargs...)

  (det(Qi)^0.25)*(det(Qj)^0.25)*(det(Qij)^-0.5)*(sill(γ)-γ(xi,xj))
end
