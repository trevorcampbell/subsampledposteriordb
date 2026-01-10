data {
  int<lower=0> I;
  int<lower=0> J;
  array[I, J] int<lower=0, upper=1> y;
}
parameters {
  real<lower=0> sigma_theta;
  vector[J] theta;
  
  real<lower=0> sigma_a;
  vector<lower=0>[I] a;
  
  real mu_b;
  real<lower=0> sigma_b;
  vector[I] b;
  real SUBIDX;
}
model {
  sigma_theta ~ cauchy(0, 2);
  theta ~ normal(0, sigma_theta);
  
  sigma_a ~ cauchy(0, 2);
  a ~ lognormal(0, sigma_a);
  
  mu_b ~ normal(0, 5);
  sigma_b ~ cauchy(0, 2);
  b ~ normal(mu_b, sigma_b);
  
  int ii = 1;
  for (i in 1 : I) {
    for (j in 1 : J){
        if (ii-0.5 <= SUBIDX && ii+0.5 >= SUBIDX){
           target += I*J*bernoulli_logit_lpmf(y[i,j] | a[i]*(theta[j] - b[i]));
        }
        ii += 1;
    }
  }
}


