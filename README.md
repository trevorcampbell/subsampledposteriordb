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

## Important Caveat
Because of how the subsampled posteriors are implemented (see below), 
`logp(x)` can be 10-10,000x faster than adding up all the individual `logp(x,i)` terms.
For that reason, **do not use these Stan posteriors to compare subsampled methods to full-data methods.** This repository is
meant to be used for comparing different subsampling methods on a single subsampled posterior, or different full-data methods on
a single full-data posterior.

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
and one could run very similar code to work with its subsampled version, `dogs-dogs_subsampled`,
where the `get_subsample_size` function in `test/main.jl` is used to get the total "data size" for this posterior:
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

This should install the `PosteriorDB.jl` package. You now will need to modify it to add the 
subsampled posteriors. Navigate to `~/.julia/artifacts` and run

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

On Mac/Windows, you may not be able to use `move.sh`, so you'll need to manually copy the files 
in `models/` and `posteriors/` to their correct locations in the julia artifacts cache.

## How It Works
Since Stan requires real-valued inputs, 
we add one parameter `SUBIDX` to each model, and loop over indices until the value of `SUBIDX` matches that index.
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

The maximum subsampling index ("dataset size") is found in each model's information JSON file in the `subsample_size` entry.
For example, for the `earn_height_subsampled` model, the `subsample_size` is `N`, one of the variables from the model code.

```
{
  "name": "earn_height_subsampled",
  "keywords": ["ARM", "Ch. 4", "stan_examples"],
  "title": "One Predictor Linear Model",
  "subsample_size": "N",
  "prior": {
    "keywords": "stan_recommended_35dbfe6"
  },
  "description": "earn ~ height",
  "urls": "https://raw.githubusercontent.com/stan-dev/example-models/master/ARM/Ch.4/earn_height.stan",
  "model_implementations": {
    "stan": {
      "model_code": "models/stan/earn_height_subsampled.stan",
      "stan_version": ">=2.26.0"
    }
  },
  "references": "gelman2006data",
  "added_date": "2020-01-17",
  "added_by": "Oliver Järnefelt",
  "licence": "BSD3"
}
```

Some other models have a subsample size that is a simple algebraic composite of two variables, like `N+K` or `N*K` or `N-K`; processing code must
be able to handle this properly. The provided `get_subsample_size` function in `test/main.jl` shows how to do this in Julia, but you'll need
to implement it yourself for other languages.

The loop over `1:N` is why this implementation is very inefficient and should not be used to compare subsampling to full-data methods.
The two figures below show the histogram of relative time to compute `logp(x)` (full data, original posterior) versus summing
over all the individual `logp(x,i)` terms (one pass over the whole dataset via subsampling):

<img width="450" alt="reltime_logp" src="https://github.com/user-attachments/assets/6379cc54-23b5-4105-87be-666e9b8660bb" />
<img width="450" alt="reltime_gradlogp" src="https://github.com/user-attachments/assets/983d99ab-e63d-4fe4-be93-def767a0d8c8" />

## Other implementation ideas

The best alternative option is to set the index as an entry in `data{...}` in the Stan model, and update that to change which data point is being considered. This would indeed be nice; it would circumvent having to loop over the entire data index set to find `SUBIDX`. Unfortunately, as far as I can tell, Stan/BridgeStan bake the data into the C++ model object directly and don't expose a reference to it, meaning that any changes to the data involve reconstructing that object from the JSON data string. This ends up being much more expensive than running a loop with mostly no-ops in search of the right index.

Another is to try to convert the real-valued `SUBIDX` to an integer in the `model{...}` block; this doesn't seem possible currently.
