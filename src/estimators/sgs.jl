# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsProcesses.jl
# ------------------------------------------------------------------

"""
    LocalSGS(; [paramaters])
    LocalSGS(localaniso; [paramaters])
    LocalSGS(method, localaniso; [paramaters])

The sequential process method introduced by Gomez-Hernandez 1993.
It traverses all locations of the geospatial domain according to a path,
approximates the conditional Gaussian distribution within a neighborhood
using simple Kriging, and assigns a value to the center of the neighborhood
by sampling from this distribution.

## Parameters

* `path`         - Process path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `36`)
* `neighborhood` - Search neighborhood (default to `:range`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)
* `init`         - Data initialization method (default to `NearestInit()`)

For each location in the process `path`, a maximum number of
neighbors `maxneighbors` is used to fit the conditional Gaussian
distribution. The neighbors are searched according to a `neighborhood`.

The `neighborhood` can be a `MetricBall`, the symbol `:range` or `nothing`.
The symbol `:range` is converted to `MetricBall(range(γ))` where `γ` is the
variogram of the Gaussian process. If `neighborhood` is `nothing`, the nearest
neighbors are used without additional constraints.

## References

* Gomez-Hernandez & Journel 1993. [Joint Sequential Simulation of
  MultiGaussian Fields](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8)

### Notes

* This method is very sensitive to the various parameters.
  Care must be taken to make sure that enough neighbors
  are used in the underlying Kriging model.
"""
@kwdef struct LocalSGS{P,N,D} <: FieldSimulationMethod
  method::Symbol = :MovingWindows
  localaniso::LocalAnisotropy
  path::P = LinearPath()
  minneighbors::Int = 1
  maxneighbors::Int = 36 # 6x6 grid cells
  neighborhood::N = :range
  distance::D = Euclidean()
end

LocalSGS(localaniso::LocalAnisotropy; kwargs...) = LocalSGS(; localaniso, kwargs...)
LocalSGS(method::Symbol, localaniso::LocalAnisotropy; kwargs...) = LocalSGS(; method, localaniso, kwargs...)

function GeoStatsProcesses.preprocess(rng, process, meth::LocalSGS, init, domain, data)
  (; path, minneighbors, maxneighbors, neighborhood, distance) = meth
  method_ = SEQSIM(path, minneighbors, maxneighbors, neighborhood, distance)
  GeoStatsProcesses.preprocess(rng, process, method_, init, domain, data)
end

#function GeoStatsProcesses.preprocess(rng, process, meth::LocalSGS, init, domain, data)
#  (; path, minneighbors, maxneighbors, neighborhood, distance) = meth
#  method = SEQSIM(path, minneighbors, maxneighbors, neighborhood, distance)
#
#  func = process.func
#
#  ## scale objects for numerical stability
#  #dom, dat, fun, neigh = _scale(domain, data, func, method.neighborhood)
#  # or ignore that for now
#  dom, dat, fun, neigh = domain, data, func, method.neighborhood
#
#  # determine search method and min/max neighbors
#  path, searcher, nmin, nmax = GeoStatsProcesses._search(dom, neigh, method)
#
#  # determine probability model
#  model, prior = GeoStatsProcesses._probmodel(process, fun)
#
#  # transform process and data
#  sdom, sdat, cache = GeoStatsProcesses._transform(process, dom, dat)
#
#  (; path, searcher, nmin, nmax, model, prior, sdom, sdat, cache, init)
#end

function GeoStatsProcesses.randsingle(rng, process, meth::LocalSGS, domain, data, preproc)
  # retrieve parameters
  (; path, searcher, nmin, nmax, model, prior, sdom, sdat, cache, init) = preproc
  (; method, localaniso) = meth

  # initialize realization and mask
  real, mask = GeoStatsProcesses.randinit(process, sdom, sdat, init)

  # realization in matrix form for efficient updates
  realization = ustrip.(stack(Tables.rowtable(real)))

  # save units of columns to restore later
  units = isnothing(sdat) ? GeoStatsProcesses._units(process) : GeoStatsProcesses._units(sdat)

  # locations with all variables already simulated
  simulated = map(all, Tables.rowtable(mask))

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, nmax)

  # retrieve variable names
  vars = keys(real)

  # variables passed to probability model
  mvars = length(vars) == 1 ? first(vars) : vars

  # restart searcher using local aniso
  ref_searcher = if searcher isa KNearestSearch
    Qi = qmat(localaniso, 1)
    anisodistance = Mahalanobis(Symmetric(Qi))
    KNearestSearch(searcher.domain, searcher.k; metric=anisodistance)
  else
    N = ndims(localaniso)
    angs = quat_to_dcm(rotation(localaniso, 1))[SOneTo(N), SOneTo(N)]'
    ranges = radii(searcher.ball)
    ranges = length(ranges) < 2 ? Tuple(ranges .* magnitude(localaniso, 1)) : ranges
    localsneigh = MetricBall(ranges, angs)
    KBallSearch(searcher.domain, searcher.k, localsneigh)
  end

  # simulation loop
  @inbounds for ind in traverse(domain, path)
    if !simulated[ind]
      # center of target location
      center = centroid(sdom, ind)

      # buffer at target location
      buffer = view(realization, :, ind)

      # modify local search
      local_searcher = if searcher isa KNearestSearch
        Qi = qmat(localaniso, ind)
        @set ref_searcher.tree.metric = Mahalanobis(Symmetric(Qi))
      else
        N = ndims(localaniso)
        angs = quat_to_dcm(rotation(localaniso, ind))[SOneTo(N), SOneTo(N)]'
        ranges = radii(searcher.ball)
        ranges = length(ranges) < 2 ? Tuple(ranges .* magnitude(localaniso, ind)) : ranges
        localsneigh = MetricBall(ranges, angs)
        searcher_ = @set ref_searcher.ball = localsneigh
        @set searcher_.tree.metric = metric(localsneigh)
      end

      # search neighbors with simulated data
      n = search!(neighbors, center, local_searcher, mask=simulated)

      if n < nmin
        # draw from prior
        GeoStatsProcesses._draw!(rng, prior, buffer)
      else
        # neighborhood with data
        neigh = let
          inds = view(neighbors, 1:n)
          ndom = view(sdom, inds)
          nmat = view(realization, :, inds)
          ntab = (; zip(vars, eachrow(nmat))...)
          georef(ntab, ndom)
        end

        # fit probability model
        hdlocalaniso = neighs_localaniso(localaniso, sdom, neigh; method)
        local_model = local_probmodel(model, localaniso, hdlocalaniso)

        # fit distribution probmodel
        fitted = local_fit(local_model, neigh, i=ind, m=1:nrow(neigh))

        # draw from conditional
        conditional = if local_status(fitted)
          GeoStatsProcesses._conditional(process, fitted, mvars, center)
        else
          prior
        end
        GeoStatsProcesses._draw!(rng, conditional, buffer)
      end

      # mark location as simulated and continue
      simulated[ind] = true
    end
  end

  # convert back to table format
  @inbounds for (i, var) in enumerate(vars)
    real[var] .= realization[i, :] * units[i]
  end

  # undo data transformations
  rdat = GeoStatsProcesses._bwdtransform(process, georef(real, sdom), cache)

  # return realization values
  values(rdat)
end

function local_probmodel(probmodel, localaniso, hdlocalaniso)
  isnothing(hdlocalaniso) ? MW_SKModel(localaniso, probmodel.fun, probmodel.mean) :
  KC_SKModel(localaniso, probmodel.fun, probmodel.mean, hdlocalaniso)
end

function Meshes._pboxes(::Type{𝔼{N}}, points) where {N}
  p = first(points)
  ℒ = lentype(p)
  cmin = fill(typemax(ℒ), N)
  cmax = fill(typemin(ℒ), N)

  for p in points
    c = getfield(coords(p), :coords)
    for i in 1:N
      cmin[i] = min(c[i], cmin[i])
      cmax[i] = max(c[i], cmax[i])
    end
  end
  Box(withcrs(p, Tuple(cmin)), withcrs(p, Tuple(cmax)))
end
