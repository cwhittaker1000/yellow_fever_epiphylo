data {
  int<lower=1> N0;         // number of days for which to impute infections
  int<lower=1> N2;         // days of observed data 
  array[N2] int cases;     // reported cases
  array[N2] int deaths;    // reported deaths -- the rows with i > N contain -1 and should be ignored
  vector[N2] death_delay;  // death delay distribution
  int EpidemicStart;       // time at which point you start fitting to deaths
  real pop;                // population 
  array[N2] real SI;       // fixed pre-calculated SI using emprical data from Neil
}

transformed data {
  vector[N2] SI_rev;              // SI in reverse order
  vector[N2] death_delay_rev;     // death delay distribution in reversed order
  
  for(i in 1:N2)                  // reversing the SI and death delay distribution
    SI_rev[i] = SI[N2 - i + 1];       // done for convenience (see code below)
  for(i in 1:N2) {
    death_delay_rev[i] = death_delay[N2 - i + 1];
  }
}


parameters {
  real<lower=0> mu;             // basic reproduction number
  real<lower=0> kappa;          // to do with prior on the reproduction number
  real<lower=0> y;              // initial infections
  real<lower=0> tau;            // prior on initial infection dispersion 
  real<lower=0.001> phi;            // overdispersion in the deaths data observational model
  // real<lower=0> ifr;            // infection fatality ratio
}

transformed parameters {
  vector<lower=0>[N2] infections = rep_vector(0, N2);   // number of infections at each timepoint
  vector<lower=0>[N2] E_deaths  = rep_vector(0, N2);    // number of deaths at each timepoint
  vector[N2] Rt = rep_vector(0, N2);           // Reproduction number over time (no adjustment)
  vector[N2] Rt_adj = Rt;                      // Adjusted Reproduction number over time (adjustment for susceptible depletion)
  
  {
    vector[N2] cumm_sum = rep_vector(0, N2);            // cumulative incidence of cases
    infections[1:N0] = rep_vector(y, N0);               // learn the number of cases in the first N0 days
    cumm_sum[2:N0] = cumulative_sum(infections[2:N0]);
    for (i in 1:N2) {
      Rt[i] = mu;
    }
    Rt_adj[1:N0] = Rt[1:N0];
    
    // Calculating infections over time
    for (i in (N0 + 1):N2) {
      real convolution = dot_product(infections[1:(i - 1)], tail(SI_rev, i - 1)); 
      cumm_sum[i] = cumm_sum[i - 1] + infections[i - 1];
      Rt_adj[i] = ((pop - cumm_sum[i]) / pop) * Rt[i];
      infections[i] = Rt_adj[i] * convolution;
    }
    // Calculating deaths over time
    E_deaths[1]= 1e-15 * infections[1];
    for (i in 2:N2){
      // E_deaths[i] = ifr * dot_product(infections, tail(death_delay_rev, i - 1));
      E_deaths[i] = dot_product(infections[1:(i - 1)], tail(death_delay_rev, i - 1));
    }
  }
}
model {
  tau ~ exponential(0.1);
  y ~ exponential(1 / tau);
  phi ~ normal(0, 5);
  kappa ~ normal(0, 0.8);
  mu ~ normal(3.28, kappa); 
  // ifr ~ normal(1, 0.1);  // change this
  deaths[EpidemicStart:N2] ~ neg_binomial_2(E_deaths[EpidemicStart:N2], phi);
}

generated quantities {
  
  vector[N2] infections0 = rep_vector(0, N2);
  vector[N2] E_deaths0  = rep_vector(0, N2);
  vector[N2] Rt_adj0 = rep_vector(0, N2);                    
  
  {
    vector[N2] cumm_sum0 = rep_vector(0, N2);
    for (i in 2:N0){
      cumm_sum0[i] = cumm_sum0[i - 1] + y; 
    }
    infections0[1:N0] = rep_vector(y, N0); 
    for (i in (N0 + 1):N2) {
      real convolution0 = dot_product(infections0[1:(i - 1)], tail(SI_rev, i - 1)); 
      cumm_sum0[i] = cumm_sum0[i - 1] + infections0[i - 1];
      Rt_adj0[i] = ((pop - cumm_sum0[i]) / pop) * mu;
      infections0[i] = Rt_adj0[i] * convolution0;
    }
    E_deaths0[1]= 1e-15 * infections0[1];
    for (i in 2:N2){
      // E_deaths0[i] = ifr * dot_product(infections0, tail(death_delay_rev, i - 1));
      E_deaths0[i] = dot_product(infections0[1:(i - 1)], tail(death_delay_rev, i - 1));
    }
  }
}
