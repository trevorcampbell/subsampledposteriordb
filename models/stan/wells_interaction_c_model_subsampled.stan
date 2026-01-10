data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> switched;
  vector[N] dist;
  vector[N] arsenic;
}
transformed data {
  // centering
  vector[N] c_dist100 = (dist - mean(dist)) / 100.0;
  vector[N] c_arsenic = arsenic - mean(arsenic);
  // interaction
  vector[N] inter = c_dist100 .* c_arsenic;
  matrix[N, 3] x = [c_dist100', c_arsenic', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
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


