data {
  int<lower=1> N; //number of observations
  int<lower=1> Dx; //number of predictors + first column = intercept alpha_0
  matrix[N,Dx] x; //covariate matrix
  int<lower=0> y[N];//response, whitefish counts
  
  int<lower=1> Np; //number of observations
  matrix[Np,Dx] xp; //covariate matrix
  int<lower=0> yp[Np];//response, whitefish counts
}

parameters {
  vector [Dx] beta;//regression parameters
  real<lower=0> phi; // the overdispersion parameter
}

transformed parameters {
  vector [Np] eta; //linear predictor
  eta = multiply(xp,beta); //use the log link
  //phi = 1. / reciprocal_phi;
}

model {
  //priors for beta and phi coefficients 

  //reciprocal_phi ~ cauchy(0.,5);
  phi ~ cauchy(0,3);
  
  for (i in 1:Dx)
    beta[i] ~ normal(0, 10);
    
  //the likelihood  - log alternatove reparametrization
  y ~   neg_binomial_2_log(multiply(x,beta), phi) ;
  
}

//define lik to be usd in k-fold CV 
generated quantities {
  vector[Np] Lik;
  vector[Np] mu;
  
  mu = exp(eta);

  for (n in 1:Np) {
    //log probability mass function 
    Lik[n] = exp(neg_binomial_2_log_lpmf(yp[n] | eta[n], phi));
  }
}

