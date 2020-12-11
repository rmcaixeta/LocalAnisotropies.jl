
function colwise(D::LocalSpatialData, ::Type{LocalAnisotropy}, i::Int, xj::Vector{Int})
	Qi = qmat(rotation(D,i),magnitude(D,i))
	d = Vector{Float64}()
	for j in xj
		Qj = qmat(rotation(D,j),magnitude(D,j))
		Qij = (Qi+Qj)/2
		local_d = Mahalanobis(Qij)
		push!(d,evaluate(local_d, xi, xj))
	end
	d
end

function colwise(D::LocalSpatialData, ::Type{LocalVariogram}, i::Int, xj::Vector{Int})
	Qi, xi = qmat(rotation(D,x),magnitude(D,i)), coordinates(obj(D),i)
	d = Vector{Float64}()
	for j in xj
		Qj, xj = qmat(rotation(D,x),magnitude(D,j)), coordinates(obj(D),j)
		dγ = sill(γ) - kcvario(D.refγ, xi, xj, Qi, Qj)
		push!(d,dγ)
	end
	d
end

function colwise(D::LocalSpatialData, ::Type{GraphDistance}, i::Int, j::Ints)
	d = dijkstra_shortest_paths(D.graph,i).dists
	view(d,j)
end


function evaluate(D::LocalSpatialData, ::Type{LocalAnisotropy}, i::Int, j::Int)
	Qi, Qj = [qmat(rotation(D,x),magnitude(D,x)) for x in (i,j)]
	Qij = (Qi+Qj)/2
	local_d = Mahalanobis(Qij)
	evaluate(local_d, xi, xj)
end

function evaluate(D::LocalSpatialData, ::Type{LocalVariogram}, i::Int, j::Int)
	Qi, Qj = [qmat(rotation(D,x),magnitude(D,x)) for x in (i,j)]
	xi, xj = [coordinates(obj(D),x) for x in (i,j)]
	sill(γ) - kcvario(D.refγ, xi, xj, Qi, Qj)
end

evaluate(D::LocalSpatialData, ::Type{GraphDistance}, i::Int, j::Ints) = colwise(GraphDistance, D, i, j)
