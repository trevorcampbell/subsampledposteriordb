#!/bin/bash

echo "You need to edit this script to replace A_REALLY_LONG_HASH with the location of your PosteriorDB.jl install." >&2
echo "Once that's done you can remove this printout and the exit 1 that follows." >&2
exit 1

cp -r models/* ~/.julia/artifacts/A_REALLY_LONG_HASH/posteriordb-1.0.0/posterior_database/models
cp -r posteriors/* ~/.julia/artifacts/A_REALLY_LONG_HASH/posteriordb-1.0.0/posterior_database/posteriors
cp -r data/* ~/.julia/artifacts/A_REALLY_LONG_HASH/posteriordb-1.0.0/posterior_database/posteriors

