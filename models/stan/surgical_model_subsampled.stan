data {
  int<lower=0> N;
  array[N] int r;
  array[N] int n;
}
parameters {
  real mu;
  real<lower=0> sigmasq;
  array[N] real b;
  real SUBIDX;
}
transformed parameters {
  real<lower=0> sigma;
  array[N] real<lower=0, upper=1> p;
  sigma = sqrt(sigmasq);
  for (i in 1 : N) {
    p[i] = inv_logit(b[i]);
  }
}
model {
  mu ~ normal(0.0, 1000.0);
  sigmasq ~ inv_gamma(0.001, 0.001);
  b ~ normal(mu, sigma);
  for (i in 1:N){
    if (i - 0.5 <= SUBIDX && i+0.5 >= SUBIDX){
      target += N*binomial_logit_lpmf(r[i]| n[i], b[i]);
      break;
    }
  }
}
generated quantities {
  real pop_mean;
  pop_mean = inv_logit(mu);
}


