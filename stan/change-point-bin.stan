
data {
  
  // dimensions
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> K;
  int id[N];
  
  // data
  vector[N] t;
  vector[N] g;
  int<lower=0, upper=1> Y[N, K];

  // hyperparameters
  vector[K] a;
  vector[K] b;
  cov_matrix[K] R;
  cov_matrix[K] S;
  
}

parameters {
  
  // mean function
  vector[K] mu;
  vector[K] alpha[M];
  vector[K] beta1;
  vector[K] beta2;
  vector[K] gamma;
  
  // covariance
  vector<lower=0>[K] sigma_0;
  vector<lower=-20,upper=10>[K] kappa;

}

transformed parameters {

  cov_matrix[K] Sigma_0 = diag_matrix(square(sigma_0));

}

model {

  // Priors
  mu ~ multi_normal(a, R);
  beta1 ~ multi_normal(b, S);
  beta2 ~ multi_normal(b, S);
  gamma ~ multi_normal(b, S);
  
  sigma_0 ~ cauchy(0, 5);
  kappa ~ uniform(-20, 10);

  // Random Intercept
  alpha ~ multi_normal(mu, Sigma_0);
  
  // Likelihood
  for (i in 1:N) {
    
    for (k in 1:K) {
      
      // change point 
      Y[i,k] ~ bernoulli_logit(alpha[id[i],k] + beta1[k]*g[i] + beta2[k]*t[i] + gamma[k]*g[i]*fdim(t[i],kappa[k]));
      
    }

    
  }


}
