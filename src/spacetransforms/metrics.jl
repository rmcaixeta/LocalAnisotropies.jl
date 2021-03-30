
function colwise(D::LocalGeoData, ::LocalAnisotropy, i::Int, xj::Vector{Int})
	Qi, ix = qmat(rotation(D,i),magnitude(D,i)), coords(D,i)
	d = Vector{Float64}()
	for j in xj
		Qj = qmat(rotation(D,j),magnitude(D,j))
		Qij = (Qi+Qj)/2
		jx = coords(D,j)
		push!(d,Distances.evaluate(Mahalanobis(Qij), ix, jx))
	end
	d
end

function colwise(D::LocalGeoData, ::LocalVariogram, i::Int, xj::Vector{Int})
	Qi, ix = qmat(rotation(D,i),magnitude(D,i)), centro(D,i)
	d = Vector{Float64}()
	for j in xj
		Qj, jx = qmat(rotation(D,j),magnitude(D,j)), centro(D,j)
		dγ = abs(sill(D.refvario) - kccov(D.refvario, ix, jx, Qi, Qj))
		push!(d,dγ)
	end
	d
end

function colwise(D::LocalGeoData, ::GraphDistance, i::Int, j::Ints)
	d = dijkstra_shortest_paths(D.graph,i).dists
	view(d,j)
end


function evaluate(D::LocalGeoData, ::LocalAnisotropy, i::Int, j::Int)
	Qi, Qj = [qmat(rotation(D,x),magnitude(D,x)) for x in (i,j)]
	ix, jx = [coords(D,x) for x in (i,j)]
	Qij = (Qi+Qj)/2
	Distances.evaluate(Mahalanobis(Qij), ix, jx)
end

function evaluate(D::LocalGeoData, ::LocalVariogram, i::Int, j::Int)
	Qi, Qj = [qmat(rotation(D,x),magnitude(D,x)) for x in (i,j)]
	ix, jx = [centro(D,x) for x in (i,j)]
	abs(sill(D.refvario) - kccov(D.refvario, ix, jx, Qi, Qj))
end

evaluate(D::LocalGeoData, ::GraphDistance, i::Int, j::Int) = colwise(GraphDistance, D, i, j)[1]
