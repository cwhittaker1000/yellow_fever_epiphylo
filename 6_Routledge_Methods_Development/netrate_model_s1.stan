data{
	int n;
	int q;
	vector[n] t;
	int NI[q];
	int z;
}

parameters {
  vector<lower=0, upper=0.01>[z] A; // need to check this upper constraint
	real<lower=0> epsilon_edge;
} 
	
model {
	real H=0; 		 
	real S=0;
	for (i in 1:z) {
		A[i] ~ normal(0.003,0.01) T[0,0.01];
	}
	epsilon_edge ~ exponential(200);
	
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
