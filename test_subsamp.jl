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
	for postnm in posteriors
		if "_subsampled" in postnm
			continue
		end
		println("Running test for $postnm")
		post = PosteriorDB.posterior(pdb, postnm)
		prb = StanProblem(post, "stan")
		postsub = PosteriorDB.posterior(pdb, postnm*"_subsampled")
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
		println(postnm*": ll = $ll lls = $lls")

		_, gll = LogDensityProblems.logdensity_and_gradient(prb, z)
		glls = 0.0
                for i = 1:subsample_sz
			_, glli = LogDensityProblems.logdensity_and_gradient(prbsub, vcat(z,i))
			glls += glli[1:end-1]
		end
		glls /= subsample_sz
		println(postnm*": gll = $gll glls = $glls")
	end
end

main()
