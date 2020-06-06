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
  //y ~ poisson_log(log_theta), log_theta = x*beta
  //y ~ poisson_log((x*beta));
    y ~ poisson_log(multiply(x,beta));
  //y ~ poisson(exp(multiply(x,beta)));
}

//define log_lik to be usd in LOO CV 
generated quantities {
  vector[N] lik;
  vector[N] y_rep;
  for (n in 1:N) {
    //probability mass function 
    //Increment target log probability density with poisson_log_lpmf( n | log_theta) 
    //dropping constant additive terms.
    lik[n] = exp(poisson_log_lpmf(y[n] | (x[n] * beta)));
    y_rep[n] = poisson_log_rng((x[n] * beta)); 
  }
}

