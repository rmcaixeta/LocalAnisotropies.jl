# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    localanisotropies(Gradients, grid, propname, w)

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
lpars = localanisotropies(Gradients, grid, :CO2, 5)
```
"""
function localanisotropies(::Type{Gradients}, obj, prop, window)
    # get dimensions
    subd = domain(obj) isa SubDomain
    dom = subd ? parent(domain(obj)) : domain(obj)
    dims = size(dom)
    N = length(dims)
    propv = rescale_min_max(Tables.getcolumn(Tables.columns(values(obj)), prop))

    # extract gradients
    img = if subd
        im = fill(NaN, dims)
        im[parentindices(domain(obj))] .= propv
        im
    else
        reshape(propv, Size(dims))
    end

    g = imgradients(img, KernelFactors.sobel, "replicate")

    quat = Array{Quaternion}(undef, size(img)) # make some better way to store it
    m = Array{Vector}(undef, size(img)) # make some better way to store it

    Threads.@threads for i in CartesianIndices(img)
        isnan(img[i]) && continue
        tensor = zeros(Float64, N, N)
        for ng in gridneighbors(img, i, window)
            for j in CartesianIndices((1:N, 1:N))
                val = g[j[1]][ng] * g[j[2]][ng]
                !isnan(val) && (tensor[j] += val)
            end
        end

        T = eigen(Symmetric(SMatrix{N,N}(tensor)))
        V = T.vectors[:, sortperm(T.values)]'
        if N == 3
            eigv = V
            det(V) < 0 && (eigv = Diagonal(SVector{3}([-1, 1, 1])) * eigv)
        else
            eigv = zeros(Float64, 3, 3)
            eigv[1:2, 1:2] = V
            eigv[1, 1] ≉ eigv[2, 2] && (eigv[1, :] .*= -1)
            eigv[3, 3] = 1.0
            eigv = SMatrix{3,3}(eigv)
        end
        q = dcm_to_quat(DCM(eigv))
        λ = sort(T.values, rev = true)
        quat[i] = q
        m[i] = λ / λ[1]
    end

    outid = findall(x -> !isnan(x), img)
    LocalAnisotropy(vec(quat[outid]), reduce(hcat, vec(m[outid])))
end


function rescale_min_max(v)
    # Extract non-NaN values
    valid_vals = filter(!isnan, v)

    # If there are no valid values, return the original vector
    if isempty(valid_vals)
        return v
    end

    # Compute min and max of non-NaN values
    min_val = minimum(valid_vals)
    max_val = maximum(valid_vals)

    # Rescale the vector, keeping NaNs unchanged
    return [isnan(x) ? NaN : (x - min_val) / (max_val - min_val) for x in v]
end
