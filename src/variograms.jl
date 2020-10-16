
"""
 Variogram function using local pars
 Local pars can be different between structures
 e.g. γns = γ1(pars,local_pars) + γ2(pars) # only vary in the first structure
"""

function get_pars(γ::Variogram)

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


function mwvario(estimator, localpar)

  #if local_par <: LocalParameters # to check scaling multistructures
  m = localpar.magnitude # maybe need to check if range = 1 or not
  N = length(m)
  P = quat_to_dcm(localpar.rotation)[SOneTo(N),SOneTo(N)]
  Λ = Diagonal(SVector{N}(one(eltype(m))./m.^2))
  Q = P*Λ*P'
  local_d = Mahalanobis(Q)

  ref_pars = get_pars(estimator.γ)

  # localvario = nothing
  # for i in 1:length(ref_pars)
  #   kwargs = ref_pars[i][:kw]
  #   kwargs[:distance] = local_d
  #   localvario += ref_pars[i][:type](kwargs...)
  # end
  kwargs = ref_pars[1][:kw]
  kwargs[:distance] = local_d
  localvario = ref_pars[1][:type](;kwargs...)

  if typeof(estimator) <: SimpleKriging
    return SimpleKriging(localvario, estimator.mean)
  elseif typeof(estimator) <: OrdinaryKriging
    return OrdinaryKriging(localvario)
  end

  #SimpleKriging(varparams.variogram, varparams.mean)
  #OrdinaryKriging(varparams.variogram)
  #estimator
end

function kcfill!(LHS, γ, X, localpars)
  nothing
end
