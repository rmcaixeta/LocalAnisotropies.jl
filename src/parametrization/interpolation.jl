
# adapt to LocalGeoData
# mean instead median for magnitude? -> mean(magnitude(lpars,neighids),dims=2)

function interpolate(lpars, searcher::NeighborSearchMethod, domain=nothing;
	power::Float64, metric=Euclidean())
	D = searcher.object
	N = ncoords(D)
	len = domain==nothing ? nelms(D) : nelms(domain)
	@assert nelms(D)==length(lpars) "searcher domain must match number of local parameters"

    quat = Array{Quaternion}(undef,len)
    m    = domain==nothing && power==0 ? lpars.magnitude : Array{Vector}(undef,len)

    Threads.@threads for i in 1:len
        icoords  = domain==nothing ? coordinates(D,i) : coordinates(domain,i)
		neighids = search(icoords, searcher)

		if length(neighids) == 0
			throw(ErrorException("zero neighbors at some location; adjust searcher"))
			# or accept missing in LocalParameters
		elseif power==0.0
			quat[i] = quatavg(rotation(lpars,neighids))
			if domain!=nothing
				mi   = magnitude(lpars,neighids)
				m[i] = mapreduce(x->quantile(view(mi,x,:),0.5), vcat, 1:N)
			end
		else
			xcoords = map(x->coordinates(D,x), neighids)
			weigths = 1 ./ (eps() .+ colwise(metric, icoords, xcoords)) .^ power
			weights = Weights(weights ./ sum(weights))
			quat[i] = quatavg(rotation(lpars,neighids), weights)
			mi      = magnitude(lpars,neighids)
			m[i]    = mapreduce(x->StatsBase.quantile(view(mi,x,:),weights,0.5), vcat, 1:N)
		end
    end

    LocalParameters(quat, m)
end

IDWpars(lpars, searcher::NeighborSearchMethod, domain; power=2.0, metric=Euclidean()) =
	interpolate(lpars, searcher, domain, power=power)

smoothpars(lpars, searcher::NeighborSearchMethod; power=0.0, metric=Euclidean()) =
	interpolate(lpars, searcher, power=power)


# Rerference: Markley, F. Landis, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman.
# "Averaging quaternions." Journal of Guidance, Control, and Dynamics 30,
# no. 4 (2007): 1193-1197.
# Link: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872.pdf

quat2vector(q) = [q.q0; q.q1; q.q2; q.q3]
tensor(q) = q' .* q
wtensor(q) = q[1] .* (q[2]' .* q[2])

function quatavg(qarr, warr=[])
	length(qarr) == 1 && (return qarr)
	#@assert in(length(warr),[n,0])
	Q = mapreduce(quat2vector, hcat, qarr)
    n = size(Q,2)
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
