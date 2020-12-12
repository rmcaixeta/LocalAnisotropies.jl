
function colwise(D::LocalGeoData, ::LocalAnisotropy, i::Int, xj::Vector{Int})
	Qi, ix = qmat(rotation(D,i),magnitude(D,i)), coords(D,i)
	d = Vector{Float64}()
	for j in xj
		Qj = qmat(rotation(D,j),magnitude(D,j))
		Qij = (Qi+Qj)/2
		local_d = Mahalanobis(Qij)
		jx = coords(D,j)
		push!(d,Distances.evaluate(local_d, ix, jx))
	end
	d
end

function colwise(D::LocalGeoData, ::LocalVariogram, i::Int, xj::Vector{Int})
	Qi, ix = qmat(rotation(D,x),magnitude(D,i)), coords(D,i)
	d = Vector{Float64}()
	for j in xj
		Qj, jx = qmat(rotation(D,x),magnitude(D,j)), coords(D,j)
		dγ = sill(γ) - kcvario(D.refvario, ix, jx, Qi, Qj)
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
	local_d = Mahalanobis(Qij)
	Distances.evaluate(local_d, ix, jx)
end

function evaluate(D::LocalGeoData, ::LocalVariogram, i::Int, j::Int)
	Qi, Qj = [qmat(rotation(D,x),magnitude(D,x)) for x in (i,j)]
	ix, jx = [coords(D,x) for x in (i,j)]
	sill(γ) - kcvario(D.refγ, ix, jx, Qi, Qj)
end

evaluate(D::LocalGeoData, ::GraphDistance, i::Int, j::Ints) = colwise(GraphDistance, D, i, j)
