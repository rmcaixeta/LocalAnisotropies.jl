
function addgraph(obj::SpatialData, lpar::LocalParameters, metric::LocalMetric,
	searcher::NeighborSearchMethod)

	D = LocalSpatialData(obj,lpar)
	addgraph!(D, metric, searcher)
end

function addgraph!(D::LocalSpatialData, metric::LocalMetric, searcher::NeighborSearchMethod)
	n = nelms(searcher.object)
	sources, dest, wgts = Vector{Int}(), Vector{Int}(), Vector{Float64}()

	for i in 1:n
		idxs = search(i, searcher)
		for j in idx
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

cleangraph!(D::LocalSpatialData) = (D = @set D.graph = nothing)
