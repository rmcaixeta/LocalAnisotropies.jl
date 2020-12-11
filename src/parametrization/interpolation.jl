

# Later adapt it to interpolation using IDW weights
function smooth(lpars, searcher::NeighborSearchMethod)
	D = searcher.object
	N, len = ncoords(D), nelms(D)

    quat = Array{Quaternion}(undef,len)
    m = Array{Vector}(undef,len)

    Threads.@threads for i in 1:len
        icoords = coordinates(D,i)
		neighids = search(icoords, searcher)
		quat[i] = quatavg(rotation(lpars,neighids))
    end

    LocalParameters(quat, lpars.magnitude)
end


# Rerference: Markley, F. Landis, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman.
# "Averaging quaternions." Journal of Guidance, Control, and Dynamics 30,
# no. 4 (2007): 1193-1197.
# Link: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872.pdf

quat2vector(q) = [q.q0; q.q1; q.q2; q.q3]
tensor(q) = q' .* q
wtensor(q) = q[1] .* (q[2]' .* q[2])

function quatavg(qarr, warr=[])
    Q = mapreduce(quat2vector, hcat, qarr)
    n = size(Q,2)
    #@assert in(length(warr),[n,0])
    W = length(warr) > 0

    A = W ? mapreduce(wtensor, +, zip(warr, eachcol(Q))) : mapreduce(tensor, +, eachcol(Q))

    # scale
    Nw = W ? sum(warr) : n
    A ./= Nw

    # compute eigenvalues and -vectors
    T = eigen(Symmetric(SMatrix{4,4}(A)))
    V = view(T.vectors[:, sortperm(T.values,rev=true)],:,1)
    Quaternion(V...)
end
