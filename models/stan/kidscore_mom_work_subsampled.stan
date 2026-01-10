data {
  int<lower=0> N;
  vector[N] kid_score;
  array[N] int mom_work;
}
transformed data {
  vector[N] work2;
  vector[N] work3;
  vector[N] work4;
  for (i in 1 : N) {
    work2[i] = mom_work[i] == 2;
    work3[i] = mom_work[i] == 3;
    work4[i] = mom_work[i] == 4;
  }
}
parameters {
  vector[4] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  for (i in 1:N){
    if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
    	target += N*normal_lpdf(kid_score[i] | beta[1] + beta[2] * work2[i] + beta[3] * work3[i]
                     + beta[4] * work4[i], sigma);
        break;
    }
  }
}


