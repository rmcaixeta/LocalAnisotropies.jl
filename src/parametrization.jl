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


# cast points as grid

struct LocalAnisotropies
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

function Gradient(img, prop)

    #img = reshape(preimg.table[:Z], Size(300,260))
    g = imgradients(img, KernelFactors.sobel, "replicate")

    w = 5
    N = length(g)
    out = Array{LocalAnisotropies}(undef,300,260) # make some better way to store it

    Threads.@threads for i in CartesianIndices(img)
        tensor = zeros(Float64,N,N)
        for ng in gridneighbors(img, i, 10)
            for j in CartesianIndices((1:N,1:N))
                tensor[j] += g[j[1]][ng]*g[j[2]][ng]
            end
        end

        T = eigen(SMatrix{N,N}(tensor))
        if N==3
            eigv = T.vectors
        else
            eigv = zeros(Float64,3,3)
            eigv[1:2,1:2] = T.vectors
            eigv[3,3] = 1.0
            eigv = SMatrix{3,3}(eigv)
        end
        q = dcm_to_quat(eigv')
        out[i] = LocalAnisotropies(q,T.values)
    end

end
