# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsTransforms.jl
# ------------------------------------------------------------------

"""
    LocalInterpolate(params ...)

Example:
MW = LocalKriging(:MovingWindows, lpars, γ)
data |> Select(:var) |> LocalInterpolate(domain, model=MW, kwargs...)

## Keyword Parameters
* `point`            - Perform interpolation on point support (default to `true`)
* `prob`             - Perform probabilistic interpolation (default to `false`)
* `get_only_weights` - Get only interpolator declustering weights
* `minneighbors`     - Minimum number of neighbors (default to `1`)
* `maxneighbors`     - Maximum number of neighbors (default to `10`)
* `neighborhood`     - Search neighborhood (default to `nothing`)
* `local_search`     - If rotate/rescale searches using local anisotropies
"""
struct LocalInterpolate{D<:Domain,GM<:GeoStatsModel,P,N} <: TableTransform
  domain::D
  model::GM
  path::P
  point::Bool
  prob::Bool
  get_only_weights::Bool
  local_search::Bool
  minneighbors::Int
  maxneighbors::Int
  neighborhood::N
end

LocalInterpolate(
  domain::Domain;
  model,
  path=LinearPath(),
  point=true,
  prob=false,
  get_only_weights=false,
  local_search=true,
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing
) = LocalInterpolate(
  domain,
  model,
  path,
  point,
  prob,
  get_only_weights,
  local_search,
  minneighbors,
  maxneighbors,
  neighborhood
)

LocalInterpolate(geoms::AbstractVector{<:Geometry}; kwargs...) = LocalInterpolate(GeometrySet(geoms); kwargs...)

isrevertible(::Type{<:LocalInterpolate}) = false

function apply(transform::LocalInterpolate, geotable::AbstractGeoTable)
  interp = localfitpredict(
    # forward arguments
    transform.model,
    geotable |> AbsoluteUnits(), # handle affine units
    transform.domain;
    path=transform.path,
    point=transform.point,
    prob=transform.prob,
    get_only_weights=transform.get_only_weights,
    local_search=transform.local_search,
    minneighbors=transform.minneighbors,
    maxneighbors=transform.maxneighbors,
    neighborhood=transform.neighborhood
  )

  interp, nothing
end

function localfitpredict(
  model::GeoStatsModel,
  dat::AbstractGeoTable,
  dom::Domain;
  path=LinearPath(),
  point=true,
  prob=false,
  get_only_weights=false,
  local_search=true,
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing
)
  # point or volume support
  pdat = point ? GeoStatsModels._pointsupport(dat) : dat

  if model isa KCModels && isnothing(model.hdlocalaniso)
    model = initmodel(model, dat, dom)
  end

  # scale objects for numerical stability
  #smodel, sdat, sdom, sneigh = GeoStatsModels._scale(model, pdat, dom, neighborhood)
  # or do not rescale by default
  smodel, sdat, sdom, sneigh = (model, pdat, dom, neighborhood)

  nobs = nrow(sdat)
  if maxneighbors > nobs || maxneighbors < 1
    maxneighbors = nobs
  end
  if minneighbors > maxneighbors || minneighbors < 1
    minneighbors = 1
  end

  neighbors = Vector{Int}(undef, maxneighbors)
  inds = traverse(sdom, path)
  predfun = prob ? GeoStatsModels._marginals ∘ predictprob : predict

  # check pars
  localaniso = smodel.localaniso
  oklocal1 = nvals(localaniso) == nvals(sdom)
  oklocal2 = typeof(localaniso) <: LocalAnisotropy
  @assert oklocal1 "number of local anisotropies must match domain points"
  @assert oklocal2 "wrong format of local anisotropies"

  # determine bounded search method
  ref_searcher = if isnothing(sneigh)
    refdist = local_search ? Mahalanobis(Symmetric(qmat(localaniso, 1))) : Euclidean()
    KNearestSearch(domain(sdat), maxneighbors; metric=refdist)
  else
    refneigh = if local_search
      N = ndims(localaniso)
      angs = quat_to_dcm(rotation(localaniso, 1))[SOneTo(N), SOneTo(N)]'
      ranges = radii(sneigh)
      ranges = length(ranges) < 2 ? Tuple(ranges .* magnitude(localaniso, 1)) : ranges
      m = MetricBall(ranges, angs)
      sneigh isa MetricBall ? m : (@set sneigh.ball = m)
    else
      sneigh
    end
    KBallSearch(domain(sdat), maxneighbors, refneigh)
  end

  if get_only_weights
    predict_weights(
      sdat,
      sdom,
      smodel,
      localaniso,
      local_search,
      inds,
      sneigh,
      ref_searcher,
      neighbors,
      minneighbors,
      point
    )
  else
    predict_variables(
      sdat,
      sdom,
      dom,
      smodel,
      localaniso,
      local_search,
      inds,
      sneigh,
      ref_searcher,
      neighbors,
      minneighbors,
      point,
      predfun
    )
  end
end

function predict_variables(
  sdat,
  sdom,
  dom,
  smodel,
  localaniso,
  local_search,
  inds,
  sneigh,
  ref_searcher,
  neighbors,
  minneighbors,
  point,
  predfun
)
  # predict variables
  cols = Tables.columns(values(sdat))
  vars = Tables.columnnames(cols)
  pred = @inbounds map(inds) do ind
    center = centroid(sdom, ind)

    ## modify search with local anisotropy
    searcher = make_local_searcher(ind, ref_searcher, sneigh, localaniso, local_search)

    nneigh = search!(neighbors, center, searcher)
    if nneigh ≥ minneighbors
      ninds = view(neighbors, 1:nneigh)
      samples = view(sdat, ninds)
      fmodel = local_fit(smodel, samples, i=ind, m=ninds)
      geom = point ? center : sdom[ind]
      vals = predfun(fmodel, vars, geom)
    else
      # missing prediction
      vals = fill(missing, length(vars))
    end
    (; zip(vars, vals)...)
  end

  # convert to original table type
  georef(pred |> Tables.materializer(values(sdat)), dom)
end

function predict_weights(
  sdat,
  sdom,
  smodel,
  localaniso,
  local_search,
  inds,
  sneigh,
  ref_searcher,
  neighbors,
  minneighbors,
  point
)
  # predict weights
  accu_weights = zeros(Float64, length(domain(sdat)))
  for ind in inds
    center = centroid(sdom, ind)

    ## modify search with local anisotropy
    searcher = make_local_searcher(ind, ref_searcher, sneigh, localaniso, local_search)

    nneigh = search!(neighbors, center, searcher)
    if nneigh ≥ minneighbors
      ninds = view(neighbors, 1:nneigh)
      samples = view(sdat, ninds)
      fmodel = local_fit(smodel, samples, i=ind, m=ninds)

      # save weights
      geom = point ? center : sdom[ind]
      wgts = get_weights(fmodel, geom)
      accu_weights[ninds] .+= abs.(wgts)
    end
  end

  accu_weights ./= sum(accu_weights)
  accu_weights
end

get_weights(m, g) =
  m isa FittedIDW ? ustrip.(GeoStatsModels.weights(m, g)) :
  m isa LocalFittedKriging ? weights(m, g).λ : GeoStatsModels.weights(m, g).λ

function make_local_searcher(ind, ref_searcher, sneigh, localaniso, local_search)
  if !local_search
    ref_searcher
  elseif isnothing(sneigh)
    Qi = qmat(localaniso, ind)
    @set ref_searcher.tree.metric = Mahalanobis(Symmetric(Qi))
  else
    N = ndims(localaniso)
    angs = quat_to_dcm(rotation(localaniso, ind))[SOneTo(N), SOneTo(N)]'
    ranges = radii(sneigh)
    ranges = length(ranges) < 2 ? Tuple(ranges .* magnitude(localaniso, ind)) : ranges
    localsneigh = MetricBall(ranges, angs)
    searcher_ = @set ref_searcher.ball = localsneigh
    @set searcher_.tree.metric = metric(localsneigh)
  end
end
