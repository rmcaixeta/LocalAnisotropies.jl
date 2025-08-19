# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsModels.jl
# ------------------------------------------------------------------

abstract type LocalKrigingModel <: GeoStatsModel end

struct MW_OKModel{G<:GeoStatsFunction} <: LocalKrigingModel
  localaniso::LocalAnisotropy
  fun::G
end

struct MW_SKModel{G<:GeoStatsFunction,V} <: LocalKrigingModel
  localaniso::LocalAnisotropy
  fun::G
  mean::V
end

struct KC_OKModel{G<:GeoStatsFunction} <: LocalKrigingModel
  localaniso::LocalAnisotropy
  fun::G
  hdlocalaniso::Union{AbstractVector,Nothing}
end

struct KC_SKModel{G<:GeoStatsFunction,V} <: LocalKrigingModel
  localaniso::LocalAnisotropy
  fun::G
  mean::V
  hdlocalaniso::Union{AbstractVector,Nothing}
end

struct OKModel{G<:GeoStatsFunction} <: LocalKrigingModel
  localaniso::LocalAnisotropy
  fun::G
end

struct SKModel{G<:GeoStatsFunction,V} <: LocalKrigingModel
  localaniso::LocalAnisotropy
  fun::G
  mean::V
end

MWModels = Union{MW_OKModel,MW_SKModel}
KCModels = Union{KC_OKModel,KC_SKModel}
KModels = Union{OKModel,SKModel}

krig_estimator(model::MWModels, localpar) = mw_estimator(model, localpar)
krig_estimator(model::KC_OKModel, localpar=nothing) = OrdinaryKriging(model.fun)
krig_estimator(model::KC_SKModel, localpar=nothing) = SimpleKriging(model.fun, model.mean)
krig_estimator(model::OKModel, localpar=nothing) = OrdinaryKriging(model.fun)
krig_estimator(model::SKModel, localpar=nothing) = SimpleKriging(model.fun, model.mean)

neighs_localaniso(model::MWModels, m) = nothing
neighs_localaniso(model::KModels, m) = nothing
neighs_localaniso(model::KCModels, m) = view(model.hdlocalaniso, m)

function neighs_localaniso(localaniso::LocalAnisotropy, dom, neigh; method::Symbol=:KernelConvolution)
  method != :KernelConvolution ? nothing : qmat(nnpars(localaniso, dom, neigh))
end

function initmodel(model::KCModels, geotable, pdomain)
  hd = grid2hd_qmat(geotable, pdomain, model.localaniso)
  la = model.localaniso
  fun = model.fun
  model isa KC_OKModel ? KC_OKModel(la, fun, hd) : KC_SKModel(la, fun, model.mean, hd)
end

"""
    LocalKriging(method, localaniso, γ, μ=nothing, localanisohd=nothing)

LocalKriging estimation solver where `method` can be :MovingWindows
or :KernelConvolution; `γ` is the variogram model and `μ` is the mean in case
simple kriging is used. `localanisohd`is only necessary for :KernelConvolution
and it is automatically passed via NN from `localaniso` if not informed.
"""
function LocalKriging(
  method::Symbol,
  localaniso::LocalAnisotropy,
  fun::GeoStatsFunction;
  mean=nothing,
  hdlocalaniso::Union{LocalAnisotropy,Nothing}=nothing
)
  if method == :MovingWindows
    isnothing(mean) ? MW_OKModel(localaniso, fun) : MW_SKModel(localaniso, fun, mean)
  elseif method == :KernelConvolution
    hd = qmat(hdlocalaniso)
    isnothing(mean) ? KC_OKModel(localaniso, fun, hd) : KC_SKModel(localaniso, fun, mean, hd)
  elseif method == :Global
    isnothing(mean) ? OKModel(localaniso, fun) : SKModel(localaniso, fun, mean)
  else
    @assert false "method must be :MovingWindows or :KernelConvolution"
  end
end

nconstraints(::LocalKrigingModel) = isnothing(mean) ? 1 : 0 # OK otherwise SK

function local_fit(model_::LocalKrigingModel, data; i, m)
  localaniso = model_.localaniso
  localpar = localpair(localaniso, i)
  model = krig_estimator(model_, localpar)

  hdlocalaniso = neighs_localaniso(model_, m)
  Qx₀ = model_ isa KCModels ? qmat(localpar...) : nothing

  # initialize Kriging system
  LHS, RHS, nfun, miss = localinitkrig(model, data, hdlocalaniso)
  FLHS = GeoStatsModels.lhsfactorize(model, LHS)
  state = KrigingState(data, LHS, RHS, FLHS, nfun, miss)

  localfitting(model, state, Qx₀, hdlocalaniso)
end

# initialize Kriging system
function localinitkrig(model::KrigingModel, data, hdlocalaniso)
  LHS, RHS = GeoStatsModels.prealloc(model, data)

  fun = model.fun
  dom = domain(data)
  tab = values(data)

  # number of function evaluations
  nobs = nelements(dom)
  nvar = nvariates(fun)
  nfun = nobs * nvar

  # find locations with missing values
  miss = GeoStatsModels.missingindices(tab)

  # set main block with pairwise evaluation
  if isnothing(hdlocalaniso)
    GeoStatsFunctions.pairwise!(LHS, fun, dom)
  else
    kcfill!(LHS, fun, dom, hdlocalaniso)
  end

  # adjustments for numerical stability
  if isstationary(fun) && !GeoStatsModels.isbanded(fun)
    GeoStatsModels.lhsbanded!(LHS, fun, dom)
  end

  # set blocks of constraints
  GeoStatsModels.lhsconstraints!(model, LHS, dom)

  # find locations with missing values
  miss = GeoStatsModels.missingindices(tab)

  # knock out entries with missing values
  GeoStatsModels.lhsmissings!(LHS, nfun, miss)

  LHS, RHS, nfun, miss
end

localfitting(model, state, q, hd) = isnothing(q) ? FittedKriging(model, state) : LocalFittedKriging(model, state, q, hd)

struct LocalFittedKriging{Q,H}
  model::KrigingModel
  state::KrigingState
  Qx₀::Q
  hdlocalaniso::H
end

FKC(m::LocalFittedKriging) = FittedKriging(m.model, m.state)

local_status(fitted::LocalFittedKriging) = issuccess(fitted.state.FHS)
local_status(fitted::FittedKriging) = GeoStatsModels.status(fitted)

predict(fitted::LocalFittedKriging, var::AbstractString, gₒ) = predict(fitted, Symbol(var), gₒ)
predict(fitted::LocalFittedKriging, var::Symbol, gₒ) = predictmean(FKC(fitted), weights(fitted, gₒ), (var,)) |> first
predict(fitted::LocalFittedKriging, vars, gₒ) = predictmean(FKC(fitted), weights(fitted, gₒ), vars)

predictprob(fitted::LocalFittedKriging, var::AbstractString, gₒ) = predictprob(fitted, Symbol(var), gₒ)

function predictprob(fitted::LocalFittedKriging, var::Symbol, gₒ)
  w = weights(fitted, gₒ)
  μ = predictmean(fitted, w, (var,)) |> first
  σ² = predictvar(fitted, w, gₒ) |> first
  ## Ordinary kriging: inflate variance in comparison to simple kriging variance; correction for multigaussian case
  #σ² = fitted.model isa OrdinaryKriging  ? σ² + 2 * w.ν[1] : σ²
  Normal(ustrip.(μ), √σ²)
end

function predictprob(fitted::LocalFittedKriging, vars, gₒ)
  w = weights(fitted, gₒ)
  μ = predictmean(fitted, w, vars)
  Σ = predictvar(fitted, w, gₒ)
  # Correct somehow for the ordinary multivariate case?
  # Σ = fitted.model isa OrdinaryKriging  ? Σ + 2 * w.ν[1] : Σ
  MvNormal(ustrip.(μ), Σ)
end

predictvar(fitted::LocalFittedKriging, weights::KrigingWeights, gₒ) =
  GeoStatsModels.predictvar(FKC(fitted), weights::KrigingWeights, gₒ)

predictmean(fitted::LocalFittedKriging, weights::KrigingWeights, var) =
  GeoStatsModels.predictmean(FKC(fitted), weights::KrigingWeights, var)

rhsconstraints!(fitted::LocalFittedKriging, pₒ) = GeoStatsModels.rhsconstraints!(FKC(fitted), pₒ)

function weights(fitted::LocalFittedKriging, gₒ)
  LHS = fitted.state.FHS
  RHS = fitted.state.RHS
  nfun = fitted.state.nfun
  miss = fitted.state.miss
  dom = domain(fitted.state.data)
  fun = fitted.model.fun

  # adjust CRS of gₒ
  gₒ′ = gₒ |> Proj(crs(dom))

  # set main blocks with pairwise evaluation
  (; Qx₀, hdlocalaniso) = fitted
  kcfill!(RHS, fun, dom, [gₒ′], Qx₀, hdlocalaniso)

  # adjustments for numerical stability
  if isstationary(fun) && !GeoStatsModels.isbanded(fun)
    GeoStatsModels.rhsbanded!(RHS, fun, dom)
  end

  # set blocks of constraints
  rhsconstraints!(fitted, gₒ′)

  # knock out entries with missing values
  GeoStatsModels.rhsmissings!(RHS, miss)

  # solve Kriging system
  W = LHS \ RHS

  # index of first constraint
  λ = @view W[begin:nfun, :]
  ν = @view W[(nfun+1):end, :]

  KrigingWeights(λ, ν)
end
