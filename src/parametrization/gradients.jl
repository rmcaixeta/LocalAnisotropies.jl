
function localparameters(::Gradients, obj, prop, window)
    # get dimensions
    dims = obj.domain.dims
    N = length(dims)

    # extract gradients (maybe need an extra method to deal with big datasets)
    factor = maximum(obj.table[!,prop]) - minimum(obj.table[!,prop])
    factor = factor < 1000 ? 100 : 1
    img = round.(Int,(factor .* reshape(obj.table[!,prop], Size(dims)) )) # temp solution
    g = imgradients(img, KernelFactors.sobel, "replicate")

    quat = Array{Quaternion}(undef,size(img)) # make some better way to store it
    m = Array{Vector}(undef,size(img)) # make some better way to store it

    Threads.@threads for i in CartesianIndices(img)
        tensor = zeros(Float64,N,N)
        for ng in gridneighbors(img, i, window)
            for j in CartesianIndices((1:N,1:N))
                tensor[j] += g[j[1]][ng]*g[j[2]][ng]
            end
        end

        T = eigen(Symmetric(SMatrix{N,N}(tensor)))
        V = T.vectors[:, sortperm(T.values)]'
        if N==3
            eigv = V
			det(V) < 0 && (eigv = Diagonal(SVector{3}([-1,1,1])) * eigv)
        else
            eigv = zeros(Float64,3,3)
            eigv[1:2,1:2] = V
			eigv[1,1] ≉ eigv[2,2] && (eigv[1,:] .*= -1)
            eigv[3,3] = 1.0
            eigv = SMatrix{3,3}(eigv)
        end
        q = dcm_to_quat(eigv)
        λ = sort(T.values, rev=true)
        quat[i] = q
        m[i] = λ/λ[1]
    end

    LocalParameters(vec(quat), reduce(hcat,vec(m)))
end
