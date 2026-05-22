data {
    int<lower=1> I; // data dimension (96 for mution data)
    int<lower=1> J; // num observations 
    int<lower=1> K; // num signatures (try 25)
    array[I, J] int<lower=0> X;
    real<lower=0> alpha0; // default = 1
    real<lower=0> a;      // default = 1
    real<lower=0> J0;     // default = 10
    real<lower=0> eps;    // default = 0.001
    real<lower=0, upper=1> lik_power; 
    // lik_power values for the mutation datasets: 
    // Lung: 0.4
    // Stomach: 0.2
    // Skin: 0.01
    // Ovary: 0.6
    // Breast: 0.2 
    // Liver: 0.6
}

transformed data {
    vector<lower=0>[I] alpha_array = rep_vector(alpha0, I);
    real<lower=0> a0 = J0 * a + 1;
    real<lower=0> b0 = eps * (a0 - 1);
}

parameters {
    matrix<lower=0>[K, J] theta;
    array[K] simplex[I] r;
    vector<lower=0>[K] mu;
}

model {
    real mutation_rate;

    for (k in 1:K) {
        mu[k] ~ inv_gamma(a0, b0);
        theta[k] ~ gamma(a, a / mu[k]);
        r[k] ~ dirichlet(alpha_array);
    }

    for (i in 1:I) {
        for (j in 1:J) {
            mutation_rate = 0;
            for (k in 1:K)
                mutation_rate += r[k,i]*theta[k,j];

            target += lik_power*poisson_lpmf(X[i,j] | mutation_rate);
        }
    }
}
