# Subsampled PosteriorDB

This repository contains data-subsampled versions of most of the posteriors in PosteriorDB.
These enable sub-indexing of posterior log densities. More specifically, if the standard
log density function for a posterior named `posterior_name` from PosteriorDB is `logp(x)`,
this repository provides a posterior named `posterior_name_subsampled` 
that returns a log density 
```
logp(x,i) where logp(x) = (1/N)(logp(x,1) + logp(x,2) + ... + logp(x,N))
```
as well as a method to compute the appropriate "data size" `N` for each posterior.
Note the scaling by `N`; this ensures that 
```
I = rand(1:N)
lphat = logp(subsampled_posterior, (x, I))
```
provides an unbiased estimate of
```
lp = logp(full_posterior, x)
```

**Important:** The implementation of a single call to `logp(x,i)` is in some cases just as slow as `logp(x)`.
Do not use these posteriors to compare subsampled methods to full-data methods. This repository is
meant to be used for comparing different methods on the same posterior, whether subsampled or not.

## Basic Usage

You can use these posteriors from any interface to Stan.
For example, from Julia, first make sure there is a directory named `stan/` in the working directory to hold compiled models.
Then to work with the full-data `dogs-dogs` posterior one could run
```
using StanLogDensityProblems, PosteriorDB, LogDensityProblems

# initialize the DB
pdb = PosteriorDB.database()

# compile the model if not already compiled, otherwise just load it
post = PosteriorDB.posterior(pdb, "dogs-dogs")
prob = StanProblem(post, "stan")

# compute the log density and gradient at a random point z
d = LogDensityProblems.dimension(prob)
z = randn(d)
LogDensityProblems.logdensity(prob, z)
LogDensityProblems.logdensity_and_gradient(prob, z)
```
and one could run very similar code to work with its subsampled version, `dogs-dogs_subsampled`:
```
postsub = PosteriorDB.posterior(pdb, "dogs-dogs_subsampled")
probsub = StanProblem(postsub, "stan")
subsample_sz = get_subsample_size(pdb, postsub)

# data index to use
i = rand(1:subsample_sz)

# append i to the end of the original state to create a valid state for the subsampled model
LogDensityProblems.logdensity(probsub, vcat(z,i))
LogDensityProblems.logdensity_and_gradient(probsub, vcat(z,i))
```

## Installation

These instructions are provided for a computer running linux; if you're on Mac or 
Windows you'll need to adjust some of the steps below but it should all still be doable.

First, make sure `julia` is installed. If you don't yet have `PosteriorDB.jl` installed, open a REPL and run
```
]
activate --temp
add PosteriorDB
```

This should install the `PosteriorDB.jl` package. You now will need to modify it manually to add the 
subsampled posteriors as well as challenge problems.
Navigate to `~/.julia/artifacts` and run

```
find . -name "posteriordb*"
```
This should return something like (where `A_REALLY_LONG_HASH` is a bunch of letters and numbers):
```
./A_REALLY_LONG_HASH/posteriordb-1.0.0
```

Open up `move.sh` in this repository, and replace `A_REALLY_LONG_HASH` with the numbers and letters from your own output.
Make sure the directory structure in that script correctly points to your installation of PosteriorDB.
Then run `./move.sh` (you may need to `chmod u+x move.sh` first if it isn't executable) to move the set of
stan models and posteriors into your installation.

To check your installation, you can navigate to the `test/` folder in this repo and run `julia --project=. main.jl`.

## How It Works
Since Stan requires real-valued inputs, 
we add one parameter `SUBIDX` to each model, and loop over indices until the value of `SUBIDX` matches that index
(and this is why the implementation is generally fairly inefficient).
For example, for the `earn_height` posterior from PosteriorDB, we create the `earn_height_subsampled` 
posterior as follows:

```
// Original model code
data {
  int<lower=0> N;
  vector[N] earn;
  vector[N] height;
}
parameters {
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  earn ~ normal(beta[1]+beta[2]*height, sigma);
}

// New subsampled model code
data {
  int<lower=0> N;
  vector[N] earn;
  vector[N] height;
}
parameters {
  vector[2] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
	for (i in 1:N){
		if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
		  target += N*normal_lpdf(earn[i] | beta[1]+beta[2]*height[i], sigma);
		  break;
		}
	}
}
```
