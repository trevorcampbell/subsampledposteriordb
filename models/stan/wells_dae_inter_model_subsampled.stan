data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> switched;
  vector[N] dist;
  vector[N] arsenic;
  vector[N] educ;
}
transformed data {
  // centering
  vector[N] c_dist100 = (dist - mean(dist)) / 100.0;
  vector[N] c_arsenic = arsenic - mean(arsenic);
  vector[N] c_educ4 = (educ - mean(educ)) / 4.0;
  // interactions
  vector[N] da_inter = c_dist100 .* c_arsenic;
  vector[N] de_inter = c_dist100 .* c_educ4;
  vector[N] ae_inter = c_arsenic .* c_educ4;
  matrix[N, 6] x = [c_dist100', c_arsenic', c_educ4', da_inter', de_inter',
                    ae_inter']';
}
parameters {
  real alpha;
  vector[6] beta;
  real SUBIDX;
}
model {
  for (i in 1:N){
    if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
    	target += N*bernoulli_logit_glm_lpmf(switched[i] | [x[i, :]], alpha, beta);
    	break;
    } 
  }
}


