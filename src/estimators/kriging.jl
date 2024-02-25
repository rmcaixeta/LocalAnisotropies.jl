# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsModels.jl and GeoStatsTransforms.jl
# ------------------------------------------------------------------

struct LocalInterpolate{D<:Domain,N,M} <: TableTransform
  domain::D
  selectors::Vector{ColumnSelector}
  models::Vector{LocalKrigingModel}
  minneighbors::Int
  maxneighbors::Int
  neighborhood::N
  distance::M
  point::Bool
  prob::Bool
end

LocalInterpolate(
  domain::Domain,
  selectors,
  models;
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing,
  distance=Euclidean(),
  point=true,
  prob=false
) = LocalInterpolate(
  domain,
  collect(ColumnSelector, selectors),
  collect(LocalKrigingModel, models),
  minneighbors,
  maxneighbors,
  neighborhood,
  distance,
  point,
  prob
)

LocalInterpolate(domain, pairs::Pair{<:Any,<:LocalKrigingModel}...; kwargs...) =
  LocalInterpolate(domain, selector.(first.(pairs)), last.(pairs); kwargs...)

isrevertible(::Type{<:LocalInterpolate}) = false

function apply(transform::LocalInterpolate, geotable::AbstractGeoTable)
  tab = values(geotable)
  cols = Tables.columns(tab)
  vars = Tables.columnnames(cols)

  domain = transform.domain
  selectors = transform.selectors
  models = transform.models
  minneighbors = transform.minneighbors
  maxneighbors = transform.maxneighbors
  neighborhood = transform.neighborhood
  distance = transform.distance
  point = transform.point
  prob = transform.prob
  path = LinearPath()

  interps = map(selectors, models) do selector, model
    svars = selector(vars)
    data = geotable[:, svars]
    localfitpredict(model, data, domain, point, prob, minneighbors, maxneighbors, neighborhood, distance, path)
  end

  newgeotable = reduce(hcat, interps)

  newgeotable, nothing
end

# only differences to original: local fit with point to estimate; check pars; add KC info
function localfitpredict(model::LocalKrigingModel, geotable::AbstractGeoTable, pdomain::Domain,
  point=true,  prob=false,  minneighbors=1,  maxneighbors=10,  neighborhood=nothing, distance=Euclidean(), path=LinearPath())

  table = values(geotable)
  ddomain = domain(geotable)
  vars = Tables.schema(table).names

  # adjust data
  data = if point
    pset = PointSet(centroid(ddomain, i) for i in 1:nelements(ddomain))
    GeoStatsModels._adjustunits(georef(values(geotable), pset))
  else
    GeoStatsModels._adjustunits(geotable)
  end

  # fix neighbors limits
  nobs = nrow(data)
  if maxneighbors > nobs || maxneighbors < 1
    maxneighbors = nobs
  end
  if minneighbors > maxneighbors || minneighbors < 1
    minneighbors = 1
  end

  # determine bounded search method
  searcher = if isnothing(neighborhood)
    # nearest neighbor search with a metric
    KNearestSearch(ddomain, maxneighbors; metric=distance)
  else
    # neighbor search with ball neighborhood
    KBallSearch(ddomain, maxneighbors, neighborhood)
  end

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, maxneighbors)

  # prediction order
  inds = traverse(pdomain, path)

  # predict function
  predfun = prob ? predictprob : predict

  # check pars
  okmeth = model.method in [:MovingWindows, :KernelConvolution]
  @assert okmeth "method must be :MovingWindows or :KernelConvolution"
  localaniso = model.localaniso
  oklocal1 = length(localaniso.rotation) == nvals(pdomain)
  oklocal2 = typeof(localaniso) <: LocalAnisotropy
  @assert oklocal1 "number of local anisotropies must match domain points"
  @assert oklocal2 "wrong format of local anisotropies"

  # add KC info
  if model.method == :KernelConvolution
    if model.hdlocalaniso != nothing
      hdlocalaniso = toqmat(model.hdlocalaniso)
    else
      hdlocalaniso = grid2hd_qmat(geotable,pdomain,model.localaniso)
    end
    model = LocalKrigingModel(model.method, model.localaniso, model.γ, model.skmean, hdlocalaniso)
  end

  # predict variable values
  function pred(var)
    map(inds) do ind
      # centroid of estimation
      center = centroid(pdomain, ind)

      # find neighbors with data
      nneigh = search!(neighbors, center, searcher)

      # predict if enough neighbors
      if nneigh ≥ minneighbors
        # final set of neighbors
        ninds = view(neighbors, 1:nneigh)

        # view neighborhood with data
        samples = view(data, ninds)

        # fit model to samples
        fmodel = local_fit(model, samples, i=ind, m=ninds)

        # save prediction
        geom = point ? center : pdomain[ind]
        predfun(fmodel, var, geom)
      else
        # missing prediction
        missing
      end
    end
  end

  pairs = (var => pred(var) for var in vars)
  newtab = (; pairs...) |> Tables.materializer(table)
  georef(newtab, pdomain)
end

function local_fit(model_::LocalKrigingModel, data; i, m)
  MW = (model_.method == :MovingWindows)
  localaniso = model_.localaniso
  localpar = (rotation(localaniso,i),magnitude(localaniso,i))
  if MW
    model = mwvario(model_, localpar)
  else
    model = model_.skmean == nothing ? GeoStatsModels.OrdinaryKriging(model_.γ) : GeoStatsModels.SimpleKriging(model_.γ, model_.skmean)
    hdlocalaniso = view(model_.hdlocalaniso, m)
    Qx₀ = qmat(localpar...)
  end

  γ = model.γ
  D = domain(data)

  # build Kriging system
  LHS = MW ? lhs(model, D) : local_lhs(model, D, hdlocalaniso)
  RHS = Vector{eltype(LHS)}(undef, size(LHS, 1))

  # factorize LHS
  FLHS = GeoStatsModels.factorize(model, LHS)

  # variance type
  VARTYPE = GeoStatsFunctions.returntype(γ, first(D), first(D))

  # record Kriging state
  state = KrigingState(data, FLHS, RHS, VARTYPE)

  # return fitted model
  FK = MW ? FittedKriging(model, state) : LocalFittedKriging(model, state, Qx₀, hdlocalaniso)
  FK
end

struct LocalFittedKriging#{M<:LocalKrigingModel,S<:KrigingState}
  model::KrigingModel
  state::KrigingState
  Qx₀
  hdlocalaniso
end

FKC(m::LocalFittedKriging) = FittedKriging(m.model, m.state)

status(fitted::LocalFittedKriging) = issuccess(fitted.state.LHS)

predict(fitted::LocalFittedKriging, var, uₒ) = predictmean(FKC(fitted), weights(fitted, uₒ), var)

function predictprob(fitted::LocalFittedKriging, var, uₒ)
  w = local_weights(fitted, uₒ)
  μ = predictmean(fitted, w, var)
  σ² = predictvar(fitted, w)
  Normal(μ, √σ²)
end

predictvar(fitted::LocalFittedKriging, weights::KrigingWeights) =
  GeoStatsModels.predictvar(FKC(fitted), weights::KrigingWeights)

predictmean(fitted::LocalFittedKriging, weights::KrigingWeights, var) =
  GeoStatsModels.predictmean(FKC(fitted), weights::KrigingWeights, var)

set_constraints_rhs!(fitted::LocalFittedKriging, pₒ) =
  GeoStatsModels.set_constraints_rhs!(FKC(fitted), pₒ)

function local_lhs(model::KrigingModel, domain, localaniso)
  γ = model.γ
  nobs = nvals(domain)
  ncons = nconstraints(model)

  # pre-allocate memory for LHS
  u = first(domain)
  T = GeoStatsFunctions.returntype(γ, u, u)
  m = nobs + ncons
  LHS = Matrix{T}(undef, m, m)

  # set variogram/covariance block
  kcfill!(LHS, γ, domain, localaniso)

  # set blocks of constraints
  GeoStatsModels.set_constraints_lhs!(model, LHS, domain)

  LHS
end

function set_local_rhs!(fitted::LocalFittedKriging, pₒ)
  localaniso = (fitted.Qx₀, fitted.hdlocalaniso)
  γ = fitted.model.γ
  X = domain(fitted.state.data)
  RHS = fitted.state.RHS

  # RHS variogram/covariance
  @inbounds for j in 1:nvals(X)
    xj = centroid(X, j)
    RHS[j] = kccov(γ, pₒ, xj, localaniso[1], localaniso[2][j])
  end

  set_constraints_rhs!(fitted, pₒ)
end


function weights(fitted::LocalFittedKriging, pₒ)
  localaniso = (fitted.Qx₀, fitted.hdlocalaniso)
  nobs = nvals(fitted.state.data)

  set_local_rhs!(fitted, pₒ)

  # solve Kriging system
  x = fitted.state.LHS \ fitted.state.RHS

  λ = view(x,1:nobs)
  ν = view(x,nobs+1:length(x))

  KrigingWeights(λ, ν)
end
