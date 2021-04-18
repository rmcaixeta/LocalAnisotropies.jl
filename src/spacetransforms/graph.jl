# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    graph(domain, localaniso, metric, searcher)
    graph(domain, localaniso, metric, refvariogram, searcher)
	graph(samples, domain, localaniso, metric, searcher)
    graph(samples, domain, localaniso, metric, refvariogram, searcher)

Create a graph connecting `domain` and `samples` points locally. Number of
edges/neighbors are defined by the `searcher` object associated to the domain.
Distance between points are based on local anisotropies `localaniso` and the
desired `metric` to calculate distance between two points. Available metrics:

* `AnisoDistance()` - averaged anisotropic distance
* `LocalVariogram()`  - variogram distance with averaged anisotropy

A reference variogram `refvariogram` is necessary if metric is `LocalVariogram()`.
"""
function graph(obj::SpatialData, lpar::LocalAnisotropy, metric::LocalMetric,
	searcher::NeighborSearchMethod)

	D = LocalGeoData(obj,lpar)
	graph!(D, metric, searcher)
end

function graph(obj::SpatialData, lpar::LocalAnisotropy, metric::LocalMetric,
	refvario::Variogram, searcher::NeighborSearchMethod)

	D = LocalGeoData(obj, lpar, refvario)
	graph!(D, metric, searcher)
end

function graph(hd::SpatialData, obj::SpatialData, lpar::LocalAnisotropy, metric::LocalMetric,
	searcher::NeighborSearchMethod)

	D = LocalGeoData(hd, obj,lpar)
	graph!(D, metric, searcher)
end

function graph(hd::SpatialData, obj::SpatialData, lpar::LocalAnisotropy, metric::LocalMetric,
	refvario::Variogram, searcher::NeighborSearchMethod)

	D = LocalGeoData(hd, obj, lpar, refvario)
	graph!(D, metric, searcher)
end

function graph!(D::LocalGeoData, metric::LocalMetric, searcher::NeighborSearchMethod)
	O = searcher.domain
	n = nelements(O)
	sources, dest, wgts = Vector{Int}(), Vector{Int}(), Vector{Float64}()

	for i in 1:n
		icoord = centroid(O,i)
		idxs   = search(icoord, searcher)
		for j in idxs
			push!(sources, i)
			push!(dest, j)
			push!(wgts, evaluate(D, metric, i, j))
		end
	end

	if sdata(D)
		h1, hn = n+1, nall(D)
		for i in h1:hn
			j = spars(D,i)
			push!(sources, i)
			push!(dest, j)
			push!(wgts, evaluate(D, metric, i, j))
		end
	end

	g = SimpleWeightedGraph(sources, dest, wgts)
	comps = length(connected_components(g))
	@assert comps==1 "There are $comps subgroups of isolated vertices. Need to increase number of neighbors"
	D = @set D.graph = g
end

cleangraph!(D::LocalGeoData) = (D = @set D.graph = nothing)
