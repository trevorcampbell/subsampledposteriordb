data {
  int<lower=1> N;
  array[N] real x;
  array[N] int k;
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  vector[N] f_tilde;
  real SUBIDX;
}
transformed parameters {
  vector[N] f;
  {
    matrix[N, N] cov = gp_exp_quad_cov(x, alpha, rho)
                       + diag_matrix(rep_vector(1e-10, N));
    matrix[N, N] L_cov = cholesky_decompose(cov);
    f = L_cov * f_tilde;
  }
}
model {
  rho ~ gamma(25, 4);
  alpha ~ normal(0, 2);
  f_tilde ~ normal(0, 1);
  
  for (n in 1:N){
  	if (n-0.5 <= SUBIDX && n+0.5 >= SUBIDX){
  		target += N*poisson_log_lpmf(k[n] | f[n]);
  		break;
  	}
  }
}


