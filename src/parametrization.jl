"""
 Extract local anisotropy from some reference input
 - Gradients 2D and 3D
 - PCA for scattered points

 Store parameters as quaternions + ratios

 # https://math.stackexchange.com/questions/2816986/plotting-an-ellipsoid-using-eigenvectors-and-eigenvalues
 # world frame rotation example
 e, x = np.linalg.eig(E)
 rot_matrix = np.array([x[:, 0],
                       x[:, 1],
                       x[:, 2]]) # or .T
 dcm_to_quat(rot_matrix) # need to test

"""


# cast points as grid or special method for pts

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
        pars[i] = LocalParameters(q,sort(T.values, rev=true))
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
