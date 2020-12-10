# spatial deformation process

# dissimilarity == anchor nbs
#               == anisodist or KC
#               == graph or nograph
#               == nb_connections (if graph)
#               == max out dims


function lmds(coords, dmetric; anchors=1500, maxoutdim=10, weights=nothing)

	n = size(coords,2)
	n<anchors && (anchors = n)

	anchor_ids = setanchors(coords,anchors,weights)
	ADM = Array{Float64}(undef,(anchor,anchor))
	dissmatrix!(ADM,coords,dmetric)

	if n > anchors
		other_ids = setdiff(1:n,anchor_ids)
		anchor_coords, M1, M3 = anchors_mds(ADM,maxoutdim)
		dim = size(anchor_coords,1)
		# give variance explained by that
		other_coords = Array{Float64}(undef,(dim,length(other_ids)))
		Threads.@threads for i in 1:length(other_ids)
			out = triangulation(i,dmetric,anchor_ids,other_ids,M1,M3)
			other_coords[:,i] .= out#[:,1]
		end
		transf_ref_coords = Array{Float64}(undef, size(coords))
		transf_ref_coords[:,anchor_ids] .= anchor_coords
		transf_ref_coords[:,other_ids] .= other_coords
	else
		M = fit(MDS, ADM, maxoutdim=maxoutdim, distances=true)
    	transf_ref_coords = transform(M)
	end
	transf_ref_coords
end


function setanchors(coords,anchors,weights)

	n = size(coords,2)
	n == anchors && (return 1:n)

	# maybe need to convert GeoStatsBase weights to StatsBase format
	args = weights == nothing ? (1:n, anchors) : (1:n, wgt, anchors)
	anchor_ids = sample(args..., replace=false)
	sort!(anchor_ids)
end

function dissmatrix!(ids, dmetric)
	n = length(ids)
	for i in 1:n
		d = dmetric(ids[i]) # dijkstra format; maybe need to generalize it
		for j in i:n
			ADM[i,j] = d[ids[j]]
			ADM[j,i] = d[ids[j]]
		end
	end
end

function anchors_mds(ADM,maxoutdim)
	nx = size(ADM,1)
	G = dmat2gram(ADM)
    F = eigen(Symmetric(G))
	λ = F.values
	maxdim = minimum([maxoutdim, nx-1, sum(λ .> 0)])
	sortidx = sortperm(F.values,rev=true)
	sortλ = λ[sortidx]
    EM = (F.vectors[:,sortidx])[:,1:maxdim]
	println("Explained variance: $(sum(sortλ[1:maxdim])/sum(sortλ[λ .> 0]))")
    sq_eigenvals = sortλ[1:maxdim].^0.5
    AM = Diagonal(sq_eigenvals)
    anchor_coords = permutedims(EM*AM)

	# to use later for another points allocations
	M1 = permutedims(EM)
	M1 ./= reshape(sq_eigenvals,1,maxdim)
	M3 = Array{Float64}(undef,(nx,1))
	mean!(M3,ADM.^2)

	anchor_coords, M1, M3
end

function triangulation(i,dmetric,anchor_ids,other_ids,M1,M3)
	d = dmetric(other_ids[i])
	M2 = mapreduce(x->d[x]^2,vcat,anchor_ids)
	-0.5*M1*(M2-M3)
end
