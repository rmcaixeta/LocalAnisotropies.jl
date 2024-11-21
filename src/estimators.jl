# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsTransforms.jl
# ------------------------------------------------------------------


"""
    LocalInterpolate(params ...)
"""
struct LocalInterpolate{D<:Domain,N,M} <: TableTransform
    domain::D
    selectors::Vector{ColumnSelector}
    models::Vector{GeoStatsModel}
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
    minneighbors = 1,
    maxneighbors = 10,
    neighborhood = nothing,
    distance = Euclidean(),
    point = true,
    prob = false,
) = LocalInterpolate(
    domain,
    collect(ColumnSelector, selectors),
    collect(GeoStatsModel, models),
    minneighbors,
    maxneighbors,
    neighborhood,
    distance,
    point,
    prob,
)

LocalInterpolate(domain, pairs::Pair{<:Any,<:GeoStatsModel}...; kwargs...) =
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
        localfitpredict(
            model,
            data,
            domain,
            point,
            prob,
            minneighbors,
            maxneighbors,
            neighborhood,
            distance,
            path,
        )
    end

    newgeotable = reduce(hcat, interps)

    newgeotable, nothing
end

function localfitpredict(
    model::GeoStatsModel,
    geotable::AbstractGeoTable,
    pdomain::Domain,
    point = true,
    prob = false,
    minneighbors = 1,
    maxneighbors = 10,
    neighborhood = nothing,
    distance = Euclidean(),
    path = LinearPath(),
)

    table = values(geotable)
    ddomain = domain(geotable)
    vars = Tables.schema(table).names

    # adjust data
    data = if point
        pset = PointSet(centroid(ddomain, i) for i = 1:nelements(ddomain))
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
        KNearestSearch(ddomain, maxneighbors; metric = distance)
    else
        # neighbor search with ball neighborhood
        KBallSearch(ddomain, maxneighbors, neighborhood)
    end

    # prediction order
    inds = traverse(pdomain, path)

    # predict function
    predfun = prob ? predictprob : predict

    # check pars
    localaniso = model.localaniso
    oklocal1 = length(localaniso.rotation) == nvals(pdomain)
    oklocal2 = typeof(localaniso) <: LocalAnisotropy
    @assert oklocal1 "number of local anisotropies must match domain points"
    @assert oklocal2 "wrong format of local anisotropies"

    if model isa LocalKrigingModel
        okmeth = model.method in [:MovingWindows, :KernelConvolution]
        @assert okmeth "method must be :MovingWindows or :KernelConvolution"
        # add KC info
        if model.method == :KernelConvolution
            if isnothing(model.hdlocalaniso)
                hdlocalaniso = grid2hd_qmat(geotable, pdomain, model.localaniso)
            else
                hdlocalaniso = toqmat(model.hdlocalaniso)
            end
            model = LocalKrigingModel(
                model.method,
                model.localaniso,
                model.γ,
                model.skmean,
                hdlocalaniso,
            )
        end
    end

    # predict variable values
    function pred(var)
        tmap(inds) do ind
            # centroid of estimation
            center = centroid(pdomain, ind)

            # find neighbors with data
            ninds = search(center, searcher)
            nneigh = length(ninds)

            # predict if enough neighbors
            if nneigh ≥ minneighbors
                # view neighborhood with data
                samples = view(data, ninds)

                # fit model to samples
                fmodel = local_fit(model, samples, i = ind, m = ninds)

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
