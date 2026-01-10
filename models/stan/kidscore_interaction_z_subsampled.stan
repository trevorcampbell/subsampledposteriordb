data {
  int<lower=0> N;
  vector[N] kid_score;
  vector[N] mom_hs;
  vector[N] mom_iq;
}
transformed data {
  // standardizing
  vector[N] z_mom_hs;
  vector[N] z_mom_iq;
  vector[N] inter;
  z_mom_hs = (mom_hs - mean(mom_hs)) / (2 * sd(mom_hs));
  z_mom_iq = (mom_iq - mean(mom_iq)) / (2 * sd(mom_iq));
  inter = z_mom_hs .* z_mom_iq;
}
parameters {
  vector[4] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  for (i in 1:N){
    if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
    	target += N*normal_lpdf(kid_score[i] | beta[1] + beta[2] * z_mom_hs[i] + beta[3] * z_mom_iq[i]
                     + beta[4] * inter[i], sigma);
        break;
    }
  }
}


