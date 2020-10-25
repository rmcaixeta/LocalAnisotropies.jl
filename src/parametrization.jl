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
  rotation::AbstractVector{Quaternion}
  magnitude::AbstractArray{AbstractFloat,2}
end

function gridneighbors(img, i::CartesianIndex, window::Int)
    dim = Size(img)
    minid(d) = max(1,i[d]-window)
    maxid(d) = min(dim[d],i[d]+window)
    idx = Tuple([minid(d):maxid(d) for d in 1:length(dim)])
    CartesianIndices(idx)
end

function gradients(preimg, prop, window)
    # get dimensions
    dims = preimg.domain.dims
    N = length(dims)

    # extract gradients (maybe need an extra method to deal with big datasets)
    img = reshape(preimg.table[prop], Size(dims))
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
        quat[i] = q
        m[i] = λ/λ[1]
    end

    LocalParameters(vec(quat), reduce(hcat,vec(m)))
end


# function geometry(ref_pt, pars)
#     # PCA on mesh points
# end


function rescale_magnitude(lp::LocalParameters, bounds1, bounds2=nothing)
    N = bounds2==nothing ? 2 : size(lp.magnitude,1)
    m = lp.magnitude
    bounds = [bounds1,bounds2]

    for i in 2:N
        mi = m[i,:]
        qx = quantile(mi, [0.05,0.95])
        mi[mi .< qx[1]] .= qx[1]
        mi[mi .> qx[2]] .= qx[2]
        r = (mi .- qx[1]) ./ (qx[2]-qx[1])
        typeof(bounds[i-1]) <: Number && (bounds[i-1]=(bounds[i-1],bounds[i-1]))
        b1, b2 = bounds[i-1]
        m[i,:] .= (b1 .+ (b2-b1) .* r)
    end

    LocalParameters(lp.rotation,m)
end

"""
Generic function that cross-validate multiple scenarios and return an optimal local field

searchtests(pdata, maxneighbors::AbstractVector{Int}; metric=Euclidean()) =
    [KNearestSearcher(pdata, n, metric=metric) for n in maxneighbors]


struct TestSet
  refimgs::AbstractVector=[(img,prop1),(img,prop2)]
  smooth::AbstractVector{Number}
  meths::AbstractVector{Symbol}=[:MovingWindows,:KernelConvolution]
  search::AbstractArray{Search}=[KNN(15),KNN(30),KNN(50)]
  variogram::AbstractVector=[(:Y,γ1),(:X,γ2)]
  ratio1::AbstractArray=[(0.2,1.0),(0.5,1.0)]
  ratio2::AbstractArray=[nothing]

  folds::Int=10
  forcegradients::Bool=false

end


function localpars(problem::EstimationProblem, t::TestSet)

    vars = [v for (v,V) in variables(problem)]
    for (img,w) in Iterators.product(t.refimgs, t.smooth)
        M = typeof(img[1]) <: RegularGrid || forcegradients
        lpars = M ? gradients(img[1],img[2],w) : geometry(img,w)

        for (r1,r2) in Iterators.product(t.ratio1, t.ratio2)
            lx = rescale_magnitude(lpars, r1, r2)

            for (v,s,e,γ) in Iterators.product(vars, t.search, t.meths, t.variogram)
                p = (variogram=γ, localpars=lx, method=e, neighborhood=s)
                solver = LocalKriging(v => p)
                cv = CrossValidation(t.folds)
                ε = cverror(solver, problem, cv)
            end
        end
    end

    # store ε and pars
    # get optimal pars
    # generic call for search or is it constrained to original pdata in case o cv?

    (variogram=(:Z,γ), neighborhood=neigh, localpars=optlocalpars, method=m)
end

"""


## Interpolate LocalParameters: NN and IDW
