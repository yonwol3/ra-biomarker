
data {
  
  // dimensions
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> K;
  int id[N];
  
  // data
  vector[K] D[N];
  vector[K] L[N];
  vector[K] U[N];
  vector[N] t;
  vector[N] g;
  vector<lower=L>[K] Y[N];

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
  vector<lower=0>[K] sigma_e;
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
  
  sigma_e ~ cauchy(0, 5);
  sigma_0 ~ cauchy(0, 5);
  kappa ~ uniform(-20, 10);

  // Random Intercept
  alpha ~ multi_normal(mu, Sigma_0);
  real eps = 1e-10;  // small nudge to prevent log(0)
  
  // Likelihood
  for (i in 1:N) {

    vector[K] eta_i;
    
    for (k in 1:K) {

      eta_i[k] = alpha[id[i],k] + beta1[k]*g[i] + beta2[k]*t[i] + gamma[k]*g[i]*fdim(t[i],kappa[k]);

      if (D[i,k] == 0) {

        target += normal_lpdf(Y[i,k] | eta_i[k], sigma_e[k]);

      } else if (D[i,k] == 1) {

        real z = (U[i,k] - eta_i[k]) / sigma_e[k];
        real p = 1.0 - Phi(z);
        target += log(fmax(p, 1e-10));  // avoid log(0)

      }
      
    }
    
  }

}
