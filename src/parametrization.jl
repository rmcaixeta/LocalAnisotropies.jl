"""
 Extract local anisotropy from some reference input
 - Gradients 2D and 3D
 - PCA for scattered points

 Store parameters as quaternions + ratios

- Need to resolve:
- best design to scale ratios easily?
    fix reference vario and scale between ]0,1]?
    fix vario range to one and scale between [min,max] ranges?
- how to assign local pars differently to each structures? mayber better to only
allow if the ratios should rescale all or given list of structures



"""



struct LocalParameters
  rotation::Quaternion
  magnitude::AbstractVector
end

function gridneighbors(img, i::CartesianIndex, window::Int)
    dim = Size(img)
    minid(d) = max(1,i[d]-window)
    maxid(d) = min(dim[d],i[d]+window)
    idx = Tuple([minid(d):maxid(d) for d in 1:length(dim)])
    CartesianIndices(idx)
end

function gradients(preimg, prop, window)

    dims = preimg.domain.dims
    N = length(dims)

    # extract gradients (maybe need an extra method to deal with big datasets)
    img = reshape(preimg.table[prop], Size(dims))
    g = imgradients(img, KernelFactors.sobel, "replicate")

    pars = Array{LocalParameters}(undef,size(img)) # make some better way to store it

    Threads.@threads for i in CartesianIndices(img)
        tensor = zeros(Float64,N,N)
        for ng in gridneighbors(img, i, window)
            for j in CartesianIndices((1:N,1:N))
                tensor[j] += g[j[1]][ng]*g[j[2]][ng]
            end
        end

        T = eigen(SMatrix{N,N}(tensor))
        V = T.vectors[:, sortperm(T.values)]'
        if N==3
            eigv = V
        else
            eigv = zeros(Float64,3,3)
            eigv[1:2,1:2] = V
            eigv[3,3] = 1.0
            eigv = SMatrix{3,3}(eigv)
        end
        q = dcm_to_quat(eigv)
        λ = sort(T.values, rev=true)
        pars[i] = LocalParameters(q,λ/λ[end])
    end
    vec(pars)

end


# function geometry(ref_pt, pars)
#     # PCA on mesh points
# end


# function EmpiricalVariogram(ref_pt, pars)
#     # local variography
# end


## Interpolate LocalParameters: NN and IDW
