"""
    deformspace(domain, localpars, metric; kwargs)
    deformspace(domain, localpars, metric, refvariogram; kwargs)
	deformspace(samples, domain, localpars, metric; kwargs)
    deformspace(samples, domain, localpars, metric, refvariogram; kwargs)
    deformspace(graphobject, metric=GraphDistance(); kwargs)
	kwargs = (anchors=1500, maxoutdim=10, weights=nothing)

Spatial deformation is a method that transforms coordinates to an isotropic
space in high dimension where the local anisotropies information are embedded
into the new distances. Traditional estimations and simulations can be performed
in this new data configuration. It's necessary to inform the `domain` object
(and the `samples` if applied). The local parameters are passed via `localpars`
and `metric` define how the local parameters are used to calculate the new
data distances. Available metrics to calculate distance between two points are:

* `LocalAnisotropy()` - averaged anisotropic distance
* `LocalVariogram()`  - variogram distance with averaged anisotropy
* `GraphDistance()`   - geodesic distances of a informed graph. More details
  of the graph construction in [`graph`](@ref) docstring.

A reference variogram `refvariogram` is necessary if metric is `LocalVariogram()`.
Additional keyword arguments are the number of `anchors` to perform landmark MDS
if the number of data is too big for a traditional MDS. `maxoutdim` is the
maximum number of dimensions to retain after MDS. `weights` are the declustering
weights of the domain data (+ samples) in order to help selecting declustered
anchor points for the transformations.
"""
function deformspace(obj::SpatialData, lpar::LocalParameters, metric::LocalMetric;
	anchors=1500, maxoutdim=10, weights=nothing)
	D = LocalGeoData(obj,lpar)
	deformspace(D, metric, anchors=anchors, maxoutdim=maxoutdim, weights=weights)
end

function deformspace(obj::SpatialData, lpar::LocalParameters, metric::LocalMetric,
	refvario::Variogram; anchors=1500, maxoutdim=10, weights=nothing)
	D = LocalGeoData(obj,lpar,refvario)
	deformspace(D, metric, anchors=anchors, maxoutdim=maxoutdim, weights=weights)
end

function deformspace(hd::SpatialData, obj::SpatialData, lpar::LocalParameters,
	metric::LocalMetric; anchors=1500, maxoutdim=10, weights=nothing)
	D = LocalGeoData(hd,obj,lpar)
	deformspace(D, metric, anchors=anchors, maxoutdim=maxoutdim, weights=weights)
end

function deformspace(hd::SpatialData, obj::SpatialData, lpar::LocalParameters,
	metric::LocalMetric, refvario::Variogram; anchors=1500, maxoutdim=10, weights=nothing)
	D = LocalGeoData(hd,obj,lpar,refvario)
	deformspace(D, metric, anchors=anchors, maxoutdim=maxoutdim, weights=weights)
end

# metric = LocalVariogram, LocalAnisotropy, GraphDistance
function deformspace(D::LocalGeoData, metric::LocalMetric=GraphDistance();
	anchors=1500, maxoutdim=10, weights=nothing)

	#@assert graph exists if GraphDistance
	dim, n = ndims(D), nall(D)
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

	outobj(D, tcoords)
end

function setanchors(n,anchors,weights)
	n == anchors && (return collect(1:n))
	# maybe need to convert GeoStatsBase weights to StatsBase format
	args = weights == nothing ? (1:n, anchors) : (1:n, weights, anchors)
	ianchors = sample(args..., replace=false)
	sort!(ianchors)
end

function dissmatrix!(ADM, D::LocalGeoData, metric::LocalMetric, iax::Vector{Int})
	n = size(ADM,1)
	for (i,ia) in enumerate(iax)
		dcols = colwise(D, metric, ia, iax[i:n])
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
	#println("Explained variance: $(sum(sortλ[1:maxdim])/sum(sortλ[sortλ .> 0]))")
    sq_eigenvals = sortλ[1:maxdim].^0.5
    AM = Diagonal(sq_eigenvals)
    atcoords = permutedims(EM*AM)

	# to use later for another points allocations
	M1 = permutedims(EM)
	M1 ./= reshape(sq_eigenvals,maxdim,1)
	M3 = Array{Float64}(undef,nx) #Array{Float64}(undef,(nx,1))
	mean!(M3,ADM.^2)

	atcoords, M1, M3
end

function triangulation(D::LocalGeoData, metric::LocalMetric,i,j,M1,M3)
	M2 = colwise(D, metric, i, j) .^ 2
	-0.5*M1*(M2-M3)
end


function outobj(D, coord)
	dom = values(obj(D))
	out = if sdata(D)
		n, hn = nvals(D), nall(D)
		h1 = n+1
		dc   = view(coord,:,1:n)
		dobj = dom isa Domain ? PointSet(dc) : georef(dom, dc)
		georef(values(sobj(D)), view(coord,:,h1:hn)), dobj
	else
		dobj = dom isa Domain ? PointSet(coord) : georef(dom, coord)
		dobj
	end
	out
end

"""
    to3d(geodata)

The spatial deformation method often return an output in dimensions higher
than 3-D. This function create a new spatial object retaining only the
first three dimensions in order to help plotting of the data.
"""
function to3d(s)
	dom = s.domain
	c = reduce(hcat,[coordinates(centroid(dom,x))[1:3] for x in 1:nelements(dom)])
	georef(values(s),PointSet(c))
end
