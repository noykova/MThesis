data {
  int<lower=1> N; //number of observations
  int<lower=1> Dx; //number of predictors + first column = intercept alpha_0
  matrix[N,Dx] x; //covariate matrix
  int<lower=0> y[N];
}

parameters {
  vector [Dx] beta;
}

model {
  //priors for alpha and beta coefficients 
  //all alpha and beta coefficients are involved in beta column vector
  //beta =[alpha0, alpha1, ..., alpha6,beta1,...,beta5]
  //Set week priors for intercept and all other covariates
  for (i in 1:Dx)
    beta[i] ~ normal(0, 10);
    
  //the likelihood  
  //y ~ bernoulli(inv_logit(x*beta)); - does not work; not stable numercally
  //y ~ bernoulli_logit((x*beta));
  y ~ bernoulli_logit(multiply(x,beta));
}

//define log_lik to be usd in LOO CV 
generated quantities {
  vector[N] lik;
  for (n in 1:N) {
    //The log Bernoulli probability mass of y given chance of success  inv_logit(x[n] * beta) 
    //take exp in order to get only likelihood
    lik[n] = exp(bernoulli_logit_lpmf(y[n] | x[n] * beta));
  }
}

