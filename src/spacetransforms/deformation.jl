
function deformspace(obj::SpatialData, lpar::LocalParameters, metric::LocalMetric;
	anchors=1500, maxoutdim=10, weights=nothing)
	D = LocalSpatialData(obj,lpar)
	deformation(D, metric, anchors=anchors, maxoutdim=maxoutdim, weights=weights)
end

# metric = LocalVariogram, LocalAnisotropy, GraphDistance
function deformspace(D::LocalSpatialData, metric::LocalMetric;
	anchors=1500, maxoutdim=10, weights=nothing)

	#@assert graph exists if GraphDistance
	dim, n = ndims(D), nvals(D)
	n < anchors && (anchors = n)

	ianchors = setanchors(n, anchors, weights)
	ADM = Array{Float64}(undef, (anchors, anchors))
	dissmatrix!(ADM, D, metric, ianchors)

	if n > anchors
		iothers = setdiff(1:n, ianchors)
		atcoords, M1, M3 = anchors_mds(ADM, maxoutdim)
		dim = size(atcoords, 1)
		otcoords = Array{Float64}(undef, (dim, length(iothers)))
		Threads.@threads for i in 1:length(iothers)
			otcoords[:,i] .= triangulation(D,metric,iothers[i],ianchors,M1,M3)
		end
		tcoords = Array{Float64}(undef, (dim,n))
		tcoords[:,ianchors] .= atcoords
		tcoords[:,iothers] .= otcoords
	else
		M = fit(MDS, ADM, maxoutdim=maxoutdim, distances=true)
    	tcoords = transform(M)
		#println("Explained variance: $(sum(sortλ[1:maxdim])/sum(sortλ[λ .> 0]))")
	end
	georef(obj(D).table, tcoords)
end

function setanchors(n,anchors,weights)
	n == anchors && (return 1:n)
	# maybe need to convert GeoStatsBase weights to StatsBase format
	args = weights == nothing ? (1:n, anchors) : (1:n, weights, anchors)
	ianchors = sample(args..., replace=false)
	sort!(ianchors)
end

function dissmatrix!(ADM, D::LocalSpatialData, metric::LocalMetric, ia::Vector{Int})
	n = nvals(D)
	for i in 1:n
		dcols = colwise(D, metric, i, ia[i:end])
		for (d,j) in zip(dcols,i:n)
			ADM[i,j] = ADM[j,i] = d
		end
	end
end

function anchors_mds(ADM, maxoutdim)
	nx = size(ADM,1)
	G = dmat2gram(ADM)
    F = eigen(Symmetric(G))
	λ = F.values
	maxdim = minimum([maxoutdim, nx-1, sum(λ .> 0)])
	sorti = sortperm(F.values,rev=true)
	sortλ = λ[sorti]
    EM = (F.vectors[:,sorti])[:,1:maxdim]
	println("Explained variance: $(sum(sortλ[1:maxdim])/sum(sortλ[λ .> 0]))")
    sq_eigenvals = sortλ[1:maxdim].^0.5
    AM = Diagonal(sq_eigenvals)
    atcoords = permutedims(EM*AM)

	# to use later for another points allocations
	M1 = permutedims(EM)
	M1 ./= reshape(sq_eigenvals,1,maxdim)
	M3 = Array{Float64}(undef,(nx,1))
	mean!(M3,ADM.^2)

	atcoords, M1, M3
end

function triangulation(D::LocalSpatialData, metric::LocalMetric,i,j,M1,M3)
	M2 = colwise(D, metric, i, j) .^ 2
	# reshape/transpose it?
	-0.5*M1*(M2-M3)
end
