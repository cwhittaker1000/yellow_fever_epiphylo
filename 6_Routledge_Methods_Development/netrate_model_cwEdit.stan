data{
	int n;           // number of cases in the linelist 
	int q;           // number of not imported cases in the linelist
	vector[n] t;     // timing of symptom onsets relative to first case
	int NI[q];       // index in the linelist of not imported cases
	int z;           // number of potential infector-infectee pairs to consider
	
	real<lower=0> epsilon_mean;  // Mean on epsilon edges
	real<lower=0> A_mean;        // Serial interval mean
	real<lower=0> A_sd;          // Serial interval SD
	real<lower=0> A_trunc_low;   // Serial interval lower bound
	real<lower=0> A_trunc_high;  // Serial interval upper bound
}

parameters {
  vector<lower=A_trunc_low, upper=A_trunc_high>[z] A; // need to check this upper constraint
                                                      // also can I just replace with a single A?
	real<lower=0> epsilon_edge;
} 
	
model {
	real H=0; 		 
	real S=0;
	for (i in 1:z) {
		A[i] ~ normal(A_mean, A_sd) T[A_trunc_low, A_trunc_high];
	}
	epsilon_edge ~ exponential(epsilon_mean);
	
	int counter=1;	
	for (i in 1:q) {
		H = 0;
		for(j in 1:(NI[i] - 1)){
			H = H + A[counter] * (t[NI[i]] - t[j]); 
			counter = counter+1;
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
