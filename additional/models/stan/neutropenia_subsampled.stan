
functions{
  vector twoCptNeutModelODE(real t, vector y, 
                            real CL, real Q, real VC, real VP, real ka,
                            real mtt, real circ0, real gamma, real alpha){
    real ktr      = 4.0 / mtt;
    real prol     = y[4] + circ0;
    real transit1 = y[5] + circ0;
    real transit2 = y[6] + circ0;
    real transit3 = y[7] + circ0;
    real circ     = fmax(machine_precision(), y[8] + circ0);
    real conc     = y[2] / VC;
    real EDrug    = alpha * conc;

    vector[8] dydt;

    // PK equations
    dydt[1] = - ka * y[1];
    dydt[2] = ka * y[1] - (CL / VC + Q / VC) * y[2] + Q / VP * y[3];
    dydt[3] = Q / VC * y[2] - Q / VP * y[3];

    // PD equations
    dydt[4] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dydt[5] = ktr * (prol - transit1);
    dydt[6] = ktr * (transit1 - transit2);
    dydt[7] = ktr * (transit2 - transit3);
    dydt[8] = ktr * (transit3 - circ);

    return dydt;
  }
  
  matrix solve_events(data int nEvent,
                      data array[] real time,
                      data array[] real amt,
                      data array[] int evid,
                      data array[] int cmt,
                      real CL, real Q, real VC, real VP, real ka,
                      real mtt, real circ0, real gamma, real alpha,
                      data real rtol, data real atol, data int max_num_steps,
                      data real eps_jitter) {
    int nCmt = 8;
    matrix[nCmt, nEvent] mass;
    vector[nCmt] y;
    real t_prev;
    
    y = rep_vector(0.0, nCmt);
    t_prev = time[1];
    
    for (i in 1:nEvent) {
      
      // Integrate from previous time to current time
      if (i > 1 && time[i] > t_prev) {
        array[1] real ts = { time[i] };
        
        array[1] vector[nCmt] sol = 
          ode_rk45_tol(twoCptNeutModelODE, y, t_prev, ts,
                       rtol, atol, max_num_steps,
                       CL, Q, VC, VP, ka,
                       mtt, circ0, gamma, alpha);
        
        y = sol[1] + rep_vector(eps_jitter, nCmt);
      }
      
      // Apply bolus dose AFTER integration
      if (evid[i] == 1) y[cmt[i]] += amt[i];
      
      mass[, i] = y;
      t_prev = time[i];
    }
    
    return mass;
  }
}

data{
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  
  array[nObsPK] int<lower=1> iObsPK;
  array[nObsPD] int<lower=1> iObsPD;
  array[nt] real<lower=0> amt;
  array[nt] int<lower=0> cmt;
  array[nt] int<lower=0> evid;
  array[nt] real<lower=0> time;
  
  vector<lower=0>[nObsPK] cObs;
  vector<lower=0>[nObsPD] neutObs;

  // data for priors
  real<lower = 0> CLPrior;
  real<lower = 0> QPrior;
  real<lower = 0> V1Prior;
  real<lower = 0> V2Prior;
  real<lower = 0> kaPrior;
  real<lower = 0> CLPriorCV;
  real<lower = 0> QPriorCV;
  real<lower = 0> V1PriorCV;
  real<lower = 0> V2PriorCV;
  real<lower = 0> kaPriorCV;
  real<lower = 0> circ0Prior;
  real<lower = 0> circ0PriorCV;
  real<lower = 0> mttPrior;
  real<lower = 0> mttPriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
  real<lower = 0> alphaPrior;
  real<lower = 0> alphaPriorCV;

  // control parameters for ODE
  real<lower = 0> rtol;
  real<lower = 0> atol;
}

transformed data{
  int<lower = 0> max_num_steps = 1000000;  // 1e6pmx
  int nCmt = 8;
  real eps_jitter = 1e-7;  // jittering to correct negative ODE solution
  
  // Event schedule passed in data is already augmented, 
  // so no need to expand in transformed data.
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  //  real<lower = 0> ka; // ka unconstrained
  real<lower = (CL / V1 + Q / V1 + Q / V2 +
		sqrt((CL / V1 + Q / V1 + Q / V2)^2 -
		     4 * CL / V1 * Q / V2)) / 2> ka; // ka > lambda_1
  real<lower = 0> mtt;
  real<lower = 0> circ0;
  real<lower = 0> alpha;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;
  real SUBIDX;
}

transformed parameters{
  array[9] real<lower = 0> parms
    = {CL, Q, V1, V2, ka, mtt, circ0, gamma, alpha};
  
  matrix[nCmt, nt] x
    = solve_events(nt, time, amt, evid, cmt,
                   CL, Q, V1, V2, ka,
                   mtt, circ0, gamma, alpha,
                   rtol, atol, max_num_steps,
                   eps_jitter);

  vector<lower=0>[nt] cHat = x[2, ]' / V1;
  vector<lower=0>[nt] neutHat = x[8, ]' + circ0;

  vector<lower=0>[nObsPK] cHatObs = cHat[iObsPK];
  vector<lower = 0>[nObsPD] neutHatObs = neutHat[iObsPD];
}

model{
  CL ~ lognormal(log(CLPrior), CLPriorCV);
  Q ~ lognormal(log(QPrior), QPriorCV);
  V1 ~ lognormal(log(V1Prior), V1PriorCV);
  V2 ~ lognormal(log(V2Prior), V2PriorCV);
  ka ~ lognormal(log(kaPrior), kaPriorCV);
  sigma ~ normal(0, 1);

  mtt ~ lognormal(log(mttPrior), mttPriorCV);
  circ0 ~ lognormal(log(circ0Prior), circ0PriorCV);
  alpha ~ lognormal(log(alphaPrior), alphaPriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaNeut ~ normal(0, 1);

  int ii = 1;
  for (i in 1:nObsPK){
  	if (ii-0.5 <= SUBIDX <= ii+0.5){
  		target += (nObsPK+nObsPD)*lognormal_lpdf(cObs[i] | log(cHatObs[i]), sigma); // observed data likelihood
  	}
  	ii += 1;
  }
  for (i in 1:nObsPD){
  	if (ii-0.5 <= SUBIDX <= ii+0.5){
  		target += (nObsPK+nObsPD)*lognormal_lpdf(neutObs[i] | log(neutHatObs[i]), sigmaNeut);
  	}
  	ii += 1;
  }
}
