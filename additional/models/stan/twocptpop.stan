
functions {

  matrix solve_two_cpt_events(int nEvent,
                              array[] real time,
                              array[] real amt,
                              array[] int evid,
                              array[] int cmt,
                              real CL, real Q, real VC, real VP, real ka) {
    matrix[3, nEvent] mass;
    vector[3] y;
    real t_prev;
    real dt;
    
    // Construct matrix for linear ODE
    real k10 = CL / VC;
    real k12 = Q  / VC;
    real k21 = Q  / VP;

    matrix[3,3] K;
    K[1,1] = -ka;      K[1,2] = 0;                K[1,3] = 0;
    K[2,1] =  ka;      K[2,2] = -(k10 + k12);     K[2,3] =  k21;
    K[3,1] =  0;       K[3,2] =  k12;             K[3,3] = -k21;
    
    y = rep_vector(0.0, 3);
    t_prev = time[1];
    
    for (i in 1:nEvent) {
      
      // Integrate from previous time to current time
      if (i > 1 && time[i] > t_prev) {
        // array[1] real ts = { time[i] };
        
        dt = time[i] - t_prev;
        y = matrix_exp(K * dt) * y + rep_vector(1e-12, 3); 
        
        // array[1] vector[3] sol = 
        //   ode_rk45(two_cpt, y, t_prev, ts, CL, Q, VC, VP, ka);
        
        // y = sol[1];
      }
      
      // Apply bolus dose AFTER integration
      if (evid[i] == 1) y[cmt[i]] += amt[i];
      
      mass[, i] = y;
      t_prev = time[i];
    }
    
    return mass;
  }
}

data {
  int<lower=1> nEvent;
  int<lower=1> nObs;
  array[nObs] int<lower=1> iObs;
  int<lower=1> nSubjects;
  array[nSubjects] int<lower=1> start;
  array[nSubjects] int<lower=1> end;
  vector<lower=0>[nSubjects] weight;

  // Event schedule
  array[nEvent] int<lower=1> cmt;
  array[nEvent] int evid;
  array[nEvent] int addl;
  array[nEvent] int ss;
  array[nEvent] real amt;
  array[nEvent] real time;
  array[nEvent] real rate;
  array[nEvent] real ii;

  vector<lower=0>[nObs] cObs;
}

transformed data {
  int nTheta = 5;
  int nCmt = 3;
  int nIIV = 5;

  array[nIIV] real prior_sd = {0.25, 0.5, 0.25, 0.5, 1};
  
  ////////////////////////////////////////////////////////////////////////
  // Expand event schedule to make each dosing event explicit
  
  array[nObs] int iObs_e;
  int obs_pos = 1;
  
  // --- compute total size ---
  int nEvent_expanded = 0;
  for (j in 1:nSubjects) {
    for (i in start[j]:end[j]) {
      nEvent_expanded += 1;
      if (evid[i] == 1)
        nEvent_expanded += addl[i];
    }
  }

  // --- allocate ---
  array[nEvent_expanded] real time_e;
  array[nEvent_expanded] real amt_e;
  array[nEvent_expanded] int evid_e;
  array[nEvent_expanded] int cmt_e;

  array[nSubjects] int start_new;
  array[nSubjects] int end_new;

  int pos = 1;

  // --- main loop ---
  for (j in 1:nSubjects) {

    start_new[j] = pos;

    int i = start[j];  // pointer in original data

    // track current ADDL sequence
    int addl_left = 0;
    real next_addl_time = 0;
    real addl_amt = 0;
    int addl_cmt = 0;
    real addl_ii = 0;

    while (i <= end[j] || addl_left > 0) {

      real t_orig = positive_infinity();
      if (i <= end[j])
        t_orig = time[i];

      real t_addl = positive_infinity();
      if (addl_left > 0)
        t_addl = next_addl_time;

      // --- choose next event (merge step) ---
      if (t_addl < t_orig) {
        // ADDL dose comes next

        time_e[pos]  = t_addl;
        amt_e[pos]   = addl_amt;
        evid_e[pos]  = 1;
        cmt_e[pos]   = addl_cmt;

        pos += 1;

        addl_left -= 1;
        if (addl_left > 0)
          next_addl_time += addl_ii;

      } else {
        // ORIGINAL event comes next

        time_e[pos]  = time[i];
        amt_e[pos]   = amt[i];
        evid_e[pos]  = evid[i];
        cmt_e[pos]   = cmt[i];
        
        if (evid[i] == 0) {
          iObs_e[obs_pos] = pos;
          obs_pos += 1;
        }

        pos += 1;

        // initialize ADDL sequence if needed
        if (evid[i] == 1 && addl[i] > 0) {
          addl_left = addl[i];
          addl_ii   = ii[i];
          addl_amt  = amt[i];
          addl_cmt  = cmt[i];
          next_addl_time = time[i] + ii[i];
        }

        i += 1;
      }
    }

    end_new[j] = pos - 1;
  }
  
  if (obs_pos - 1 != nObs) {
    reject("Mismatch in number of observations after expansion");
  }
}

parameters {
  // Population parameters
  real<lower=0> CL_pop;
  real<lower=0> Q_pop;
  real<lower=0> VC_pop;
  real<lower=0> VP_pop;

  // ka constrained to be > lambda_1
  real<lower=(CL_pop / VC_pop + Q_pop / VC_pop + Q_pop / VP_pop +
              sqrt((CL_pop / VC_pop + Q_pop / VC_pop + Q_pop / VP_pop)^2 -
                   4 * CL_pop / VC_pop * Q_pop / VP_pop)) / 2> ka_pop;

  // Inter-individual variability
  vector<lower=0>[nIIV] omega;
  array[nSubjects, nTheta] real<lower=0> theta;

  real<lower=0> sigma;
}

transformed parameters {
  vector<lower=0>[nTheta] theta_pop
    = to_vector({CL_pop, Q_pop, VC_pop, VP_pop, ka_pop});

  row_vector<lower=0>[nEvent_expanded] concentration;
  matrix<lower=0>[nCmt, nEvent_expanded] mass;

  // Individual parameters
  vector<lower=0>[nSubjects] CL =
    to_vector(theta[, 1]) .* exp(0.75 * log(weight / 70));
  vector<lower=0>[nSubjects] Q  =
    to_vector(theta[, 2]) .* exp(0.75 * log(weight / 70));
  vector<lower=0>[nSubjects] VC =
    to_vector(theta[, 3]) .* (weight / 70);
  vector<lower=0>[nSubjects] VP =
    to_vector(theta[, 4]) .* (weight / 70);
  vector<lower=0>[nSubjects] ka =
    to_vector(theta[, 5]);

  for (j in 1:nSubjects) {
    mass[, start_new[j]:end_new[j]] =
      solve_two_cpt_events(end_new[j] - start_new[j] + 1,
                           time_e[start_new[j]:end_new[j]],
                           amt_e[start_new[j]:end_new[j]],
                           evid_e[start_new[j]:end_new[j]],
                           cmt_e[start_new[j]:end_new[j]],
                           CL[j], Q[j], VC[j], VP[j], ka[j]);
    
    concentration[start_new[j]:end_new[j]] =
      mass[2, start_new[j]:end_new[j]] / VC[j];
  }
}

model {
  // priors
  CL_pop ~ lognormal(log(10), prior_sd[1]);
  Q_pop  ~ lognormal(log(15), prior_sd[2]);
  VC_pop ~ lognormal(log(35), prior_sd[3]);
  VP_pop ~ lognormal(log(105), prior_sd[4]);
  ka_pop ~ lognormal(log(2.5), prior_sd[5]);
  sigma  ~ normal(0, 0.5);
  omega  ~ normal(0, 0.5);

  // interindividual variability
  for (j in 1:nSubjects) {
    theta[j, ] ~ lognormal(log(theta_pop), omega);
  }

  // likelihood
  cObs ~ lognormal(log(concentration[iObs_e]), sigma);
}
