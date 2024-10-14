# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# grid neighborhood for gradients
function gridneighbors(img, i::CartesianIndex, window::Int)
    dim = Size(img)
    minid(d) = max(1, i[d] - window)
    maxid(d) = min(dim[d], i[d] + window)
    idx = Tuple([minid(d):maxid(d) for d = 1:length(dim)])
    CartesianIndices(idx)
end

# nearest grid data id to some other sample data
function grid2hd_ids(pdata, pdomain)
    hd = [ustrip.(to(centro(pdata, i))) for i = 1:nvals(pdata)]
    grid = [ustrip.(to(centro(pdomain, i))) for i = 1:nvals(pdomain)]

    tree = KDTree(grid)
    idxs, dists = nn(tree, hd)

    [i for i in idxs]
end

# nearest grid data mahalanobis matrix to some other sample data
function grid2hd_qmat(pdata, pdomain, localaniso)
    idxs = grid2hd_ids(pdata, pdomain)
    [qmat(rotation(localaniso, i), magnitude(localaniso, i)) for i in idxs]
end
