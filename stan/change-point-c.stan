
data {

  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> K;

  matrix[N,K] Y;
  vector[N] t;
  int id[N];

  vector[K] a;
  vector[K] b_1;
  vector[K] b_2;

  cov_matrix[K] R;
  cov_matrix[K] S_1;
  cov_matrix[K] S_2;
  cov_matrix[K] V_e;
  cov_matrix[K] V_0;

  int<lower=K + 1> u_e;
  int<lower=K + 1> u_0;

}

parameters {

  vector[K] mu;
  matrix[M,K] alpha;
  matrix[2,K] beta;
  vector<lower=-20,upper=5>[K] kappa;
  cov_matrix[K] Sigma_e;
  cov_matrix[K] Sigma_0;

} 

model {
  
  int idx[N,K];
  matrix[N,K] eta;

  // Random Effects
  for (j in 1:M) {

    alpha[j,1:K] ~ multi_normal(mu, Sigma_0);

  }

  // Likelihood
  for (i in 1:N) {

    for (k in 1:K) {

      // change point 
      eta[i,k] = alpha[id[i],k] + beta[1,k]*t[i] + beta[2,k]*fdim(t[i],kappa[k]);

    }

    Y[i,1:K] ~ multi_normal(eta[i,1:K], Sigma_e);
    
  }

  // Priors
  mu[1:K] ~ multi_normal(a, R);
  beta[1,1:K] ~ multi_normal(b_1, S_1);
  beta[2,1:K] ~ multi_normal(b_2, S_2);

  kappa[1] ~ uniform(-20, 10);
  kappa[2] ~ uniform(-20, 10);
  kappa[3] ~ uniform(-20, 10);
  kappa[4] ~ uniform(-20, 10);
  kappa[5] ~ uniform(-20, 10);
  kappa[6] ~ uniform(-20, 10);

  Sigma_e ~ wishart(u_e, V_e);
  Sigma_0 ~ wishart(u_0, V_0);

}
