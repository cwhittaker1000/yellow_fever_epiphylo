data{
	int n;                       // number of cases in the linelist 
	int q;                       // number of not imported cases in the linelist
	vector[n] t;                 // timing of symptom onsets relative to first case
	array[n, n] real<lower=0> d; // matrix of distances between cases
	int NI[q];                   // index in the linelist of not imported cases
	int z;                       // number of potential infector-infectee pairs to consider
	
	real<lower=0> epsilon_mean;        // mean on epsilon edges
  // real<lower=0> spatial_scale;  // mean on spatial scale parameter
	real<lower=0> A_mean;              // serial interval mean
	real<lower=0> A_sd;                // serial interval SD
	real<lower=0> A_trunc_low;         // serial interval lower bound
	real<lower=0> A_trunc_high;        // serial interval upper bound
	real<lower=0> spatial_scale_mean;  // mean on spatial scale parameter
}

parameters {
	real<lower=A_trunc_low, upper=A_trunc_high> a;
	real<lower=0> epsilon_edge;
	real<lower=0> spatial_scale;
} 

transformed parameters {
  vector<lower=A_trunc_low, upper=A_trunc_high>[z] A;
  for (i in 1:z) {
    A[i] = a;
  }
}
	
model {
	real H=0; 		 
	real S=0;
	// for (i in 1:z) {
	// 	A[i] ~ normal(A_mean, A_sd) T[A_trunc_low, A_trunc_high];
	// }
  a ~ normal(A_mean, A_sd) T[A_trunc_low, A_trunc_high];
	epsilon_edge ~ exponential(epsilon_mean);
	spatial_scale ~ exponential(spatial_scale_mean);	
	
	int counter=1;	
	for (i in 1:q) {
		H = 0;
		for(j in 1:(NI[i] - 1)){
      H = H + A[counter] * (t[NI[i]] - t[j]) * exp(-spatial_scale * d[i, j]); // larger spatial scale makes long-distance transmission less likely
			counter = counter + 1;
		}
    target += log(H + epsilon_edge);
  }
  counter=1;
	for (i in 1:q) {
	  for(k in 1:(NI[i]-1)){
		  target += (-A[counter] * ((t[NI[i]] - t[k])^2) * 0.5);
			counter = counter+1;
		}
	}
	target += -epsilon_edge; 
}

generated quantities {
  real beta;
  beta = exponential_rng(spatial_scale_mean); // 1 divided by spatial scale mean is the mean
}
