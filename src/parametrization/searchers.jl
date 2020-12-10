

function gridneighbors(img, i::CartesianIndex, window::Int)
    dim = Size(img)
    minid(d) = max(1,i[d]-window)
    maxid(d) = min(dim[d],i[d]+window)
    idx = Tuple([minid(d):maxid(d) for d in 1:length(dim)])
    CartesianIndices(idx)
end

function grid2hd(pdata,pdomain,localpars)
  hd = coordinates(pdata)
  grid = coordinates(pdomain)

  tree = KDTree(grid)
  idxs, dists = nn(tree, hd)

  [qmat(localpars.rotation[i],localpars.magnitude[:,i]) for i in idxs]
end

function grid2hd_ids(pdata,pdomain)
  hd = coordinates(pdata)
  grid = coordinates(pdomain)

  tree = KDTree(grid)
  idxs, dists = nn(tree, hd)

  [i for i in idxs]
end
