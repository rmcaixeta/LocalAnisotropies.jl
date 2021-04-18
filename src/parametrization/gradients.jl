# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    localanisotropies(Gradients(), grid, propname, w)

Extract `LocalAnisotropy` from a reference scenario. The `propname` variable
from the cartesian `grid` object is scanned and local gradients are extracted
from it. These gradients are smoothed within a squared/cubic window sized
2`w` x 2`w` (x 2`w`) in order to return an ellipse/ellipsoid along the
directions with more/less variability (`w` must be an `Int` value and corresponds
to the number of grid cells). A preliminar magnitude is also extracted from it
and can be later rescaled to proper limits using [`rescale_magnitude`](@ref)
function.

## Example

```julia
lpars = localanisotropies(Gradients(), grid, :CO2, 5)
```
"""
function localanisotropies(::Gradients, obj, prop, window)
    # get dimensions
    dims  = obj.domain.dims
    N     = length(dims)
	propv = Tables.getcolumn(Tables.columns(values(obj)), prop)

    # extract gradients (maybe need an extra method to deal with big datasets)
    factor = maximum(propv) - minimum(propv)
    factor = factor < 1000 ? 100 : 1
    img = round.(Int,(factor .* reshape(propv, Size(dims)) )) # temp solution
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

    LocalAnisotropy(vec(quat), reduce(hcat,vec(m)))
end
