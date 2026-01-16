using StanLogDensityProblems, LogDensityProblems, PosteriorDB

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
        z = randn(d)
        ll = LogDensityProblems.logdensity(prb, z)
		lls = 0.0
        for i = 1:subsample_sz
			lls += LogDensityProblems.logdensity(prbsub, vcat(z,i))
		end
		lls /= subsample_sz
		println(postnm*": ll = $(round(ll,digits=2)) lls = $(round.(lls,digits=2)) err = $(round(abs(ll-lls),digits=2)) relerr = $(round(abs(ll-lls)/abs(ll),digits=2))")

		_, gll = LogDensityProblems.logdensity_and_gradient(prb, z)
		glls = zeros(d)
        for i = 1:subsample_sz
			_, glli = LogDensityProblems.logdensity_and_gradient(prbsub, vcat(z,i))
			glls += glli[1:end-1]
		end
		glls /= subsample_sz
		idcs = rand(1:length(gll), 5)
		println(postnm*": gll = $(round.(gll[idcs],digits=2)) glls = $(round.(glls[idcs],digits=2)) err = $(round(sqrt(sum((gll-glls).^2)),digits=2)) relerr = $(round(sqrt(sum((gll-glls).^2))/sqrt(sum(gll.^2)),digits=2))")
	end
end

main()
