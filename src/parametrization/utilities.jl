

function rescale_magnitude(lp::LocalParameters, r1, r2=nothing; clip=[0.05,0.95])
    N = size(lp.magnitude,1)
    m = lp.magnitude
    bounds = [r1,r2]

    for i in 2:N
		bounds[i-1] == nothing && continue
        mi = m[i,:]
        qx = quantile(mi, clip)
        mi[mi .< qx[1]] .= qx[1]
        mi[mi .> qx[2]] .= qx[2]
        r = (mi .- qx[1]) ./ (qx[2]-qx[1])
        typeof(bounds[i-1]) <: Number && (bounds[i-1]=(bounds[i-1],bounds[i-1]))
        b1, b2 = bounds[i-1]
        m[i,:] .= (b1 .+ (b2-b1) .* r)
    end

    LocalParameters(lp.rotation,m)
end


function localpars2vtk(vtkfile,coords,lpars; dir=1,magnitude=:r1)
	dim = size(coords,1)
	n = size(coords,2)
	xyz = dim == 2 ? vcat(coords, zeros(Float64,1,n)) : coords
	verts = [MeshCell( VTKCellTypes.VTK_VERTEX, [i]) for i in 1:n]
	ijk = zeros(Float64, 3, n, 1, 1)
	for x in 1:n
		dcm = quat_to_dcm(lpars.rotation[x])
        f = dcm[3,3] < 0 ? -1 : 1
		ijk[1, x, 1, 1] = f*dcm[dir,1]
		ijk[2, x, 1, 1] = f*dcm[dir,2]
		ijk[3, x, 1, 1] = f*dcm[dir,3]
	end
	scale = magnitude == :r1 ? 2 : 3
	outfiles = vtk_grid(vtkfile,xyz,(verts)) do vtk
		vtk["orient"] = ijk
		vtk["scale"] = (1 + eps()) .- lpars.magnitude[scale,:]
	end
end

toqmat(lp) = [qmat(lp.rotation[i],lp.magnitude[:,i]) for i in 1:length(lp.rotation)]
