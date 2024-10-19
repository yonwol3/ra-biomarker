data {

	int<lower=1> N; 	//total number of observations
	int<lower=1> M;        // total number of unique individuals
	real Y[N];            //Outcome 
	vector[N] t;      // time 
        vector[N] g;    // group-indicator for each observation
        int id[N];     // indexing unique individual
        real a;       // prior mean for beta1:3 and mu parameters
        real r;       // prior sd for beta1:3 and mu parameters 
}

Parameters  {

        real alpha[M];  // random intercepts for each unique individual
        real mu;   // mean of random intercepts
	real beta1;  // coefficient for time 
        real beta2;  // coefficient for RA-group
        real beta3;   // coefficient for the interaction term 
        real kappa;  // change point time 
        real sigma_e;    // standard deviation for error
        real sigma_0;     // standard deviation for random-intercept

 }

model {
           real eta[N];
	// Random_intercept
          alpha ~ normal (mu, sigma_0);
          
	// Likelihood
         for ( i in N ) {
             eta[i]= alpha[id[i]] + beta1*t[i] +beta2*g[i] +beta3*g[i]*fdim(t[i], kappa);

            } 
	 Y~ normal(eta, sigma_e); 

         // Priors 

          mu~normal(a,r)
	  beta1~ normal (a,r);
          beta2~ normal (a,r);
          beta3~ normal (a,r);
          kappa ~uniform (-20, 5); 
          sigma_e ~ cauchy(0,10);
          sigma_0 ~ cauchy (0,10);

}  



