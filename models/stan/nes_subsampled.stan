data {
  int<lower=0> N;
  vector[N] partyid7;
  vector[N] real_ideo;
  vector[N] race_adj;
  vector[N] educ1;
  vector[N] gender;
  vector[N] income;
  array[N] int age_discrete;
}
transformed data {
  vector[N] age30_44; // age as factor
  vector[N] age45_64;
  vector[N] age65up;
  
  for (n in 1 : N) {
    age30_44[n] = age_discrete[n] == 2;
    age45_64[n] = age_discrete[n] == 3;
    age65up[n] = age_discrete[n] == 4;
  }
}
parameters {
  vector[9] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  // vectorization
  for (n in 1:N){
    if (n-0.5 <= SUBIDX && n+0.5 >= SUBIDX){
      target += N*normal_lpdf(partyid7[n] | beta[1] + beta[2] * real_ideo[n] + beta[3] * race_adj[n]
                    + beta[4] * age30_44[n] + beta[5] * age45_64[n]
                    + beta[6] * age65up[n] + beta[7] * educ1[n] + beta[8] * gender[n]
                    + beta[9] * income[n], sigma);
      break;
    }
  }
}


