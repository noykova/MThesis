data {
  int<lower=1> N; //number of observations
  int<lower=1> Dx; //number of predictors + first column = intercept alpha_0
  matrix[N,Dx] x; //covariate matrix
  int<lower=0> y[N];//response, whitefish counts
}

parameters {
  vector [Dx] beta;//regression parameters
  real<lower=0> phi; // the overdispersion parameter
}

transformed parameters {
  vector [N] eta; //linear predictor
  eta = multiply(x,beta); //use the log link
  //phi = 1. / reciprocal_phi;
}

model {
  //priors for beta and phi coefficients 

  //reciprocal_phi ~ cauchy(0.,5);
  phi ~ cauchy(0,3);
  
  for (i in 1:Dx)
    beta[i] ~ normal(0, 10);
    
  //the likelihood  - log alternatove reparametrization
  y ~   neg_binomial_2_log(eta, phi) ;
  
}

//define lik to be usd in k-fold CV 
generated quantities {
  vector[N] y_rep;
  vector[N] lik;
  vector[N] mu;
  
  mu = exp(eta);

  for (n in 1:N) {
    //log probability mass function 
    lik[n] = exp(neg_binomial_2_log_lpmf(y[n] | eta[n], phi));
    
    //random generator to check accuracy 
    //y_rep[n] = neg_binomial_2_rng(mu[n], phi);
    // eta[n] = log(mean mu), mu=log(eta)
    if (eta[n] > 15) { 
            //15<29*log(2)=20.10127
            // To avoid erros like the below during the warmup. 
            // neg_binomial_2_rng: Random number that came from gamma distribution is 3.02668e+39, but must be less than 1.07374e+09  
          // https://groups.google.com/forum/#!topic/stan-users/4g2hbwtRELQ 
     y_rep[n] = -1;
     } else { 
       y_rep[n] = neg_binomial_2_rng(mu[n], phi); 
       //y_rep[n] = neg_binomial_2_log_rng(eta[n], phi); 
     }

  }
}

