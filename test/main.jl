using StanLogDensityProblems, LogDensityProblems, PosteriorDB

function get_reasonable_point(prb, d)
	x = randn(d)
	γ = 1e-6
	for i=1:100
		f, g = LogDensityProblems.logdensity_and_gradient(prb, x)
		xp = x + γ*g
		fp = LogDensityProblems.logdensity(prb, xp)
		if fp >= f + 0.5*γ*sum(g.^2)
			x = xp
			γ *= 2.0
		else
			γ /= 2.0
		end
	end
	return x
end
		

function get_subsample_size(pdb, subsampled_posterior)
	modelnamesub = PosteriorDB.info(subsampled_posterior)["model_name"]
	datanamesub = PosteriorDB.info(subsampled_posterior)["data_name"]
	modelsub = PosteriorDB.model(pdb, modelnamesub)
	datasub = PosteriorDB.load(PosteriorDB.dataset(pdb, datanamesub))
	subsampleszstr = PosteriorDB.info(modelsub)["subsample_size"]
	subsamplesznms = split(subsampleszstr, "*")
	subsample_sz = 1
	for var in subsamplesznms
		subsample_sz *= datasub[var]
	end
	return subsample_sz
end

function main()
	pdb = PosteriorDB.database()
	posteriors = PosteriorDB.posterior_names(pdb)
	for subpostnm in posteriors
		# only look at subsampled posteriors
		if !occursin("_subsampled", subpostnm)
			continue
		end
		postnm = subpostnm[1:end-11]
		println("Running test for $postnm")
		post = PosteriorDB.posterior(pdb, postnm)
		prb = StanProblem(post, "stan")
		postsub = PosteriorDB.posterior(pdb, subpostnm)
		prbsub = StanProblem(postsub, "stan")
		subsample_sz = get_subsample_size(pdb, postsub)
		d = LogDensityProblems.dimension(prb)
		z = get_reasonable_point(prb, d)
        ll = LogDensityProblems.logdensity(prb, z)
		lls = 0.0
        for i = 1:subsample_sz
			lls += LogDensityProblems.logdensity(prbsub, vcat(z,i))
		end
		lls /= subsample_sz
		println(postnm*": ll = $(round(ll,sigdigits=2)) lls = $(round.(lls,sigdigits=2)) err = $(round(abs(ll-lls),sigdigits=2)) relerr = $(round(abs(ll-lls)/abs(ll),sigdigits=2))")

		ll, gll = LogDensityProblems.logdensity_and_gradient(prb, z)
		glls = zeros(d)
		lls = 0.0
        for i = 1:subsample_sz
			lli, glli = LogDensityProblems.logdensity_and_gradient(prbsub, vcat(z,i))
			glls += glli[1:end-1]
			lls += lli
		end
		glls /= subsample_sz
		lls /= subsample_sz
		idcs = rand(1:length(gll), 5)
		println(postnm*": gll = $(round.(gll[idcs],sigdigits=2)) glls = $(round.(glls[idcs],sigdigits=2)) err = $(round(sqrt(sum((gll-glls).^2)),sigdigits=2)) relerr = $(round(sqrt(sum((gll-glls).^2))/sqrt(sum(gll.^2)),sigdigits=2))")
		println(postnm*": ll = $(round(ll,sigdigits=2)) lls = $(round.(lls,sigdigits=2)) err = $(round(abs(ll-lls),sigdigits=2)) relerr = $(round(abs(ll-lls)/abs(ll),sigdigits=2))")
	end
end

main()
