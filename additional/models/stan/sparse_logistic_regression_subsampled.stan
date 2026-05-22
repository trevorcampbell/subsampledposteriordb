data {
  int<lower=1> N; // Number of data
  int<lower=1> M; // Number of covariates
  matrix[N, M] X;
  array[N] int<lower=0, upper=1> y;
}

transformed data {
  real sigma = 2;    // effective likelihood standard deviation
  real m0 = 3;           // Expected number of large slopes
  real slab_scale = 2;    // Scale for large slopes
  real slab_scale2 = square(slab_scale);
  real slab_df = 4;      // Effective degrees of freedom for large slopes
  real half_slab_df = 0.5 * slab_df;
}

parameters {
  vector[M] beta_tilde;
  vector<lower=0>[M] lambda;
  real<lower=0> c2_tilde;
  real<lower=0> tau_tilde;
  real alpha;
  real SUBIDX;
}

transformed parameters {
  vector[M] beta;
  {
    real tau0 = (m0 / (M - m0)) * (sigma / sqrt(1.0 * N));
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0)

    // c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scale2)
    // Implies that marginally beta ~ student_t(slab_df, 0, slab_scale)
    real c2 = slab_scale2 * c2_tilde;

    vector[M] lambda_tilde =
      sqrt( c2 * square(lambda) ./ (c2 + square(tau) * square(lambda)) );

    // beta ~ normal(0, tau * lambda_tilde)
    beta = tau * lambda_tilde .* beta_tilde;
  }
}

model {
  beta_tilde ~ normal(0, 1);
  lambda ~ cauchy(0, 1);
  tau_tilde ~ cauchy(0, 1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  alpha ~ normal(0, 10);
  for (i in 1:N){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += N*bernoulli_lpmf(y[i] | inv_logit(X[i,:]*beta + alpha));
  		break;
  	}
  }
}
