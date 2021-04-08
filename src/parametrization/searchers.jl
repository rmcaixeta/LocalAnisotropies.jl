
# grid neighborhood for gradients
function gridneighbors(img, i::CartesianIndex, window::Int)
    dim = Size(img)
    minid(d) = max(1,i[d]-window)
    maxid(d) = min(dim[d],i[d]+window)
    idx = Tuple([minid(d):maxid(d) for d in 1:length(dim)])
    CartesianIndices(idx)
end

# nearest grid data id to some other sample data
function grid2hd_ids(pdata,pdomain)
  hd   = [coordinates(centroid(pdata, i)) for i in 1:nelements(pdata)]
  grid = [coordinates(centroid(pdomain, i)) for i in 1:nelements(pdomain)]

  tree = KDTree(grid)
  idxs, dists = nn(tree, hd)

  [i for i in idxs]
end

# nearest grid data mahalanobis matrix to some other sample data
function grid2hd_qmat(pdata,pdomain,localpars)
  idxs = grid2hd_ids(pdata,pdomain)
  [qmat(rotation(localpars,i),magnitude(localpars,i)) for i in idxs]
end
