# Subsampled PosteriorDB

This repository contains modified versions of most of the posteriors in PosteriorDB 
that enable sub-indexing of posterior log densities. Specifically, if the standard
log density function for a posterior from PosteriorDB is $$\log p(x)$$, then this 
repository provides $$\log p(x, i)$$ where $$\log p(x) = \sum_{i=1}^N \log p(x,i)$$.



You can use these from any interface to Stan;
for example, from Julia one can run the following code to work with the `dogs-dogs` posterior
and its subsampled version, `dogs-dogs_subsampled`:

```
post = PosteriorDB.posterior(pdb, "dogs-dogs")
prob = StanProblem(post, "stan")
postsub = PosteriorDB.posterior(pdb, subpostnm)
prbsub = StanProblem(postsub, "stan")
subsample_sz = get_subsample_size(pdb, postsub)
d = LogDensityProblems.dimension(prb)


but the instructions in this repository will pertain specifically to using them from 
the Julia interface

of Stan's
many interany interface to 
Stan, from any language,
but the instructions in this repost


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

To check


