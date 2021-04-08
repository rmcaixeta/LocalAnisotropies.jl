"""
    localparameters(Geometric(), searcher; reg=false)

Extract `LocalParameters` from a searcher object. Based on the spatial data and
the parameters from the searcher object, a group of neighbor points is collected
at each place. PCA is applied to their spatial coordinates, returning eigenvectors
that represent a local ellipse/ellipsoid aligned with the directions with more
clustered points. The magnitude is the ratio of the PCA eigenvalues. It's useful
to extract local parameters from 3-D surface points. If `reg = true` and data is
in 3-D, only the minimum direction is extracted from PCA; the other axes are
calculated (main direction will be always along max dip direction and the
intermediate axis will be orthogonal to the others).

## Example

```julia
searcher = KNearestSearch(pts3dsurface, 10)
prelpars = localparameters(Geometric(), searcher) # extract local pars at surface
lpars = idwpars(prelpars, searcher, grid) # interpolate local pars to a grid
```
"""
function localparameters(::Geometric, searcher::NeighborSearchMethod;
	reg::Bool=false)
	X = coordinates(centroid(searcher.object))
	N, len = size(X)

    quat = Array{Quaternion}(undef,len)
    m = Array{Vector}(undef,len)

    Threads.@threads for i in 1:len
        icoord = view(X,:,i)
		neighids = search(icoord, searcher)
		λ, v = pca(view(X,:,neighids),reg)

		det(v) < 0 && (v = Diagonal(SVector{3}([-1,1,1])) * v)

        q = dcm_to_quat(v)
        quat[i] = q
        m[i] = λ/λ[1]
    end

	# deal with -1 eigvals
	m = reduce(hcat,vec(m))
	posm = m .> 0
	if sum(posm) < N*len
		for d in 2:N
			posd = view(posm,d,:)
			sum(posd) == len && continue
			minx = minimum(view(m,d,findall(posd)))
			m[d,findall(.!posd)] .= minx
		end
	end

    LocalParameters(quat, m)
end

function pca(X,reg)
	N = size(X,1)
	M = fit(PCA, X, maxoutdim=N, pratio=1)
	λ = principalvars(M)
	nv = length(λ)
	v = N == 3 ? projection(M) : vcat(projection(M),[0 0 1][1:nv])

	if nv == 1
		vx = [-v[2,1]; v[1,1]; v[3,1]]
		v = hcat(v,vx,cross(v[:,1],vx))
		append!(λ,[1.,1.])
	elseif nv == 2 && N==3
		v = hcat(v,cross(v[:,1],v[:,2]))
		push!(λ,1.)
	end

	# if reg && N==3
	# 	# use maxdip vector and cross(maxdip,normal) as the main vectors
	# end

	λ[1:N], SMatrix{3,3}(v)
end
