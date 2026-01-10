data {
  int<lower=0> N;
  int<lower=0> n_groups;
  int<lower=0> n_scenarios;
  array[N] int<lower=1, upper=n_groups> group_id;
  array[N] int<lower=1, upper=n_scenarios> scenario_id;
  vector[N] y;
}
parameters {
  vector[n_groups] a;
  vector[n_scenarios] b;
  real mu_a;
  real mu_b;
  real<lower=0, upper=100> sigma_a;
  real<lower=0, upper=100> sigma_b;
  real<lower=0, upper=100> sigma_y;
  real SUBIDX;
}
transformed parameters {
  vector[N] y_hat;
  
  for (i in 1 : N) {
    y_hat[i] = a[group_id[i]] + b[scenario_id[i]];
  }
}
model {
  mu_a ~ normal(0, 1);
  a ~ normal(10 * mu_a, sigma_a);
  
  mu_b ~ normal(0, 1);
  b ~ normal(10 * mu_b, sigma_b);
  
  for (n in 1:N){
    if (n-0.5 <= SUBIDX && n+0.5 >= SUBIDX){
      target += N*normal_lpdf(y[n] | y_hat[n], sigma_y);
      break;
    }
  }
}


