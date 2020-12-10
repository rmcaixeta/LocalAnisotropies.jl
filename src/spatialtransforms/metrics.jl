



function makegraph(ref_coords::AbstractArray{<:Number,2},neigh_type,neigh_val,anchor)

	# use weights for setanchors
	nb_points = size(ref_coords,2)
	@assert neigh_type in ["knn","inrange"] "Invalid neighborhood type"
	idxs, dists = _get_neighbors(ref_coords, neigh_type, neigh_val)
	sources,destinations,weights = [Int64[],Int64[],Float64[]]

	for i in 1:nb_points
		for j in idxs[i][idxs[i].!=i]
			push!(sources, i)
			push!(destinations, j)
			push!(weights, euclidean(view(ref_coords,:,i),view(ref_coords,:,j)))
		end
	end

	g = SimpleWeightedGraph(convert(Array{Int64,1}, sources), convert(Array{Int64,1}, destinations), weights)
	comps = length(connected_components(g))
	@assert comps==1 string("There are ",comps," subgroups of isolated vertices. Need to increase number of neighbors")
	g
end
