using StanLogDensityProblems, LogDensityProblems, PosteriorDB
using Plots

function get_reasonable_point(prb, d)
	x = randn(d)
	γ = 1e-6
	for i=1:30
		f, g = LogDensityProblems.logdensity_and_gradient(prb, x)
		xp = x + γ*g
		fp = LogDensityProblems.logdensity(prb, xp)
		if fp >= f + 0.5*γ*sum(g.^2)
			x = xp
			γ *= 2.0
		else
			γ /= 2.0
		end
		x += 1e-2*randn(d) # add a bit of noise
	end
	return x
end
		

function get_subsample_size(pdb, subsampled_posterior)
	modelnamesub = PosteriorDB.info(subsampled_posterior)["model_name"]
	datanamesub = PosteriorDB.info(subsampled_posterior)["data_name"]
	modelsub = PosteriorDB.model(pdb, modelnamesub)
	datasub = PosteriorDB.load(PosteriorDB.dataset(pdb, datanamesub))
	subsampleszstr = PosteriorDB.info(modelsub)["subsample_size"]
	subsample_sz = -1
	if occursin("*", subsampleszstr)
		subsamplesznms = split(subsampleszstr, "*")
		subsample_sz = 1
		for var in subsamplesznms
			subsample_sz *= datasub[var]
		end
	elseif occursin("+", subsampleszstr)
		subsamplesznms = split(subsampleszstr, "+")
		subsample_sz = 0
		for var in subsamplesznms
			subsample_sz += datasub[var]
		end
	elseif occursin("-", subsampleszstr)
		subsamplesznms = split(subsampleszstr, "-")
		subsample_sz = datasub[subsamplesznms[1]] - datasub[subsamplesznms[2]]
	else
		subsample_sz = datasub[subsampleszstr]
	end
	return subsample_sz
end

function main()
	pdb = PosteriorDB.database()
	posteriors = PosteriorDB.posterior_names(pdb)
	ii = 0
	trels = []
	tgradrels = []
	for subpostnm in posteriors[1:10]
		ii += 1
		# only look at subsampled posteriors
		if !occursin("_subsampled", subpostnm)
			continue
		end
		# create logprob objects
		postnm = subpostnm[1:end-11]
		println("Running test $ii for $postnm")
		post = PosteriorDB.posterior(pdb, postnm)
		prb = StanProblem(post, "stan")
		postsub = PosteriorDB.posterior(pdb, subpostnm)
		prbsub = StanProblem(postsub, "stan")
		subsample_sz = get_subsample_size(pdb, postsub)
		d = LogDensityProblems.dimension(prb)
		# create two test points (compute ll diffs, stan doesn't guarantee constants are the same)
		z1 = get_reasonable_point(prb, d)
		z2 = get_reasonable_point(prb, d)
		# run functions one time to ensure compilation before timing
		LogDensityProblems.logdensity(prb, z1)
		LogDensityProblems.logdensity(prbsub, vcat(z1,1))
		LogDensityProblems.logdensity_and_gradient(prb, z1)
		LogDensityProblems.logdensity_and_gradient(prbsub, vcat(z1,1))
		# compare logdensity and subsampled average with timing
		t_full = time_ns()
		ll = LogDensityProblems.logdensity(prb, z1) - LogDensityProblems.logdensity(prb, z2)
		t_full = (time_ns() - t_full)/1e9
		t_sub = time_ns()
		lls = 0.0
		ll2s = 0.0
		for i = 1:subsample_sz
			lldiff = LogDensityProblems.logdensity(prbsub, vcat(z1,i)) - LogDensityProblems.logdensity(prbsub, vcat(z2,i))
			lls += lldiff
			ll2s += lldiff^2
		end
		lls /= subsample_sz
		ll2s /= subsample_sz
		t_sub = (time_ns() - t_sub)/1e9
                abserr = abs(ll-lls)
		relerr = abs(ll-lls)/abs(ll)
		if (abserr > 1e-6 || relerr > 1e-6)
			println(postnm*": ll = $(round(ll,sigdigits=2)) lls = $(round(lls,sigdigits=2)) var = $(round(ll2s-lls^2,sigdigits=2)) t_sub/t_full = $(round(t_sub/t_full,sigdigits=2)) err = $(round(abserr,sigdigits=2)) relerr = $(round(relerr,sigdigits=2))")
		end
		push!(trels, t_sub/t_full)

		t_full = time_ns()
		_, gll = LogDensityProblems.logdensity_and_gradient(prb, z1)
		t_full = (time_ns() - t_full)/1e9
		t_sub = time_ns()
		glls = zeros(d)
		gll2s = 0.0
		for i = 1:subsample_sz
			_, glli = LogDensityProblems.logdensity_and_gradient(prbsub, vcat(z1,i))
			glls += glli[1:end-1]
			gll2s += sum(glli[1:end-1].^2)
		end
		glls /= subsample_sz
		gll2s /= subsample_sz
		t_sub = (time_ns() - t_sub)/1e9
		idcs = rand(1:length(gll), 5)
		abserr = sqrt(sum((gll-glls).^2))
		relerr = sqrt(sum((gll-glls).^2))/sqrt(sum(gll.^2))
		if (abserr > 1e-6 || relerr > 1e-6)
			println(postnm*": gll = $(round.(gll[idcs],sigdigits=2)) glls = $(round.(glls[idcs],sigdigits=2)) var = $(round(gll2s-sum(glls.^2),sigdigits=2))  t_sub/t_full = $(round(t_sub/t_full,sigdigits=2)) err = $(round(abserr,sigdigits=2)) relerr = $(round(relerr,sigdigits=2))")
		end
		push!(tgradrels, t_sub/t_full)
	end

	p1 = histogram(trels, xscale=:log10, bins=30, xlabel="Relative Time Increase when Subsampling", ylabel="Count", legend=false, title="Relative Time (sum log p_i) / log p")
	p2 = histogram(tgradrels, xscale=:log10, bins=30, xlabel="Relative Time Increase when Subsampling", ylabel="Count", legend=false, title ="Relative Time (sum ∇log p_i) / ∇log p")
	savefig(p1, "reltime_logp.png", dpi=300)
	savefig(p2, "reltime_gradlogp.png", dpi=300)
end

main()
