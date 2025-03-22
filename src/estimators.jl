# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsTransforms.jl
# ------------------------------------------------------------------

"""
    LocalInterpolate(params ...)
"""
struct LocalInterpolate{D<:Domain,GM<:GeoStatsModel,P,N,M} <: TableTransform
  domain::D
  model::GM
  path::P
  point::Bool
  prob::Bool
  minneighbors::Int
  maxneighbors::Int
  neighborhood::N
  distance::M
end

LocalInterpolate(
  domain::Domain;
  model,
  path=LinearPath(),
  point=true,
  prob=false,
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing,
  distance=Euclidean()
) = LocalInterpolate(domain, model, path, point, prob, minneighbors, maxneighbors, neighborhood, distance)

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
    neighbors=true,
    minneighbors=transform.minneighbors,
    maxneighbors=transform.maxneighbors,
    neighborhood=transform.neighborhood,
    distance=transform.distance
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
  neighbors=true,
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing,
  distance=Euclidean()
)
  # point or volume support
  pdat = point ? GeoStatsModels._pointsupport(dat) : dat

  if model isa KCModels && isnothing(model.hdlocalaniso)
    model = initmodel(model, dat, dom)
  end

  # scale objects for numerical stability
  smodel, sdat, sdom, sneigh = GeoStatsModels._scale(model, pdat, dom, neighborhood)

  nobs = nrow(sdat)
  if maxneighbors > nobs || maxneighbors < 1
    maxneighbors = nobs
  end
  if minneighbors > maxneighbors || minneighbors < 1
    minneighbors = 1
  end

  # determine bounded search method
  searcher = if isnothing(sneigh)
    # nearest neighbor search with a metric
    KNearestSearch(domain(sdat), maxneighbors; metric=distance)
  else
    # neighbor search with ball neighborhood
    KBallSearch(domain(sdat), maxneighbors, sneigh)
  end

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, maxneighbors)

  # traverse domain with given path
  inds = traverse(sdom, path)

  # prediction function
  predfun = prob ? GeoStatsModels._marginals ∘ predictprob : predict

  # check pars
  localaniso = smodel.localaniso
  oklocal1 = nvals(localaniso) == nvals(dom)
  oklocal2 = typeof(localaniso) <: LocalAnisotropy
  @assert oklocal1 "number of local anisotropies must match domain points"
  @assert oklocal2 "wrong format of local anisotropies"

  # predict variables
  cols = Tables.columns(values(sdat))
  vars = Tables.columnnames(cols)
  pred = @inbounds map(inds) do ind
    # centroid of estimation
    center = centroid(sdom, ind)

    # # modify search with local anisotropy
    # searcher_ = if isnothing(sneigh)
    #     Qi = qmat(localaniso, ind)
    #     anisodistance = Mahalanobis(Symmetric(Qi))
    #     KNearestSearch(sdom, maxneighbors; metric = anisodistance)
    # else
    #     # change angles here later
    #     searcher
    # end

    # find neighbors with data
    nneigh = search!(neighbors, center, searcher)

    # predict if enough neighbors
    if nneigh ≥ minneighbors
      # final set of neighbors
      ninds = view(neighbors, 1:nneigh)

      # view neighborhood with data
      samples = view(sdat, ninds)

      # fit model to samples
      fmodel = local_fit(smodel, samples, i=ind, m=ninds)

      # save prediction
      geom = point ? center : sdom[ind]
      vals = predfun(fmodel, vars, geom)
    else
      # missing prediction
      vals = fill(missing, length(vars))
    end
    (; zip(vars, vals)...)
  end

  # convert to original table type
  georef(pred |> Tables.materializer(values(sdat)), sdom)
end
