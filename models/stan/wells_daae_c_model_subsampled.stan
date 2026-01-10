data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> switched;
  vector[N] dist;
  vector[N] arsenic;
  vector[N] assoc;
  vector[N] educ;
}
transformed data {
  vector[N] c_dist100 = (dist - mean(dist)) / 100.0;
  vector[N] c_arsenic = arsenic - mean(arsenic);
  vector[N] da_inter = c_dist100 .* c_arsenic;
  vector[N] educ4 = educ / 4.0;
  matrix[N, 5] x = [c_dist100', c_arsenic', da_inter', assoc', educ4']';
}
parameters {
  real alpha;
  vector[5] beta;
  real SUBIDX;
}
model {
  for (i in 1:N){
    if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
    	target += N*bernoulli_logit_glm_lpmf(switched[i] | x[i, :], alpha, beta);
    	break;
    } 
  }
}


