
function addgraph(obj::GeoData, lpar::LocalParameters, metric::LocalMetric,
	searcher::NeighborSearchMethod)

	D = LocalGeoData(obj,lpar)
	addgraph!(D, metric, searcher)
end

function addgraph!(D::LocalGeoData, metric::LocalMetric, searcher::NeighborSearchMethod)
	O = searcher.object
	n = nelms(O)
	sources, dest, wgts = Vector{Int}(), Vector{Int}(), Vector{Float64}()

	for i in 1:n
		icoord = coordinates(O,i)
		idxs = search(icoord, searcher)
		for j in idx
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

	# make with hd, id being i+n

	g = SimpleWeightedGraph(sources, dest, wgts)
	comps = length(connected_components(g))
	@assert comps==1 "There are $comps subgroups of isolated vertices. Need to increase number of neighbors"
	D = @set D.graph = g
end

cleangraph!(D::LocalGeoData) = (D = @set D.graph = nothing)
