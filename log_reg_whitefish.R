
setwd("D:/Neli/MThesis/Ex_7_trials")

if (.Platform$OS.type != "windows") {
  if (Sys.getenv("RSTUDIO") == "1") {
    windows <- function( ... ) X11( ... )
  } else {
    windows <- dev.new
  }
}


library("rstan")
library(matrixStats)
options(mc.cores = parallel::detectCores())

library("raster")
library(Rcpp)


# Load the data
whitefish.dat = read.table("data_whitefish.txt", header=TRUE, sep="\t")
whitefish.raster = read.table("predraster_whitefish.txt", header=TRUE, sep="\t")

## The coordinate system in the data is ETRS89 (http://www.euref.eu/euref_egrs.html) which is 
## a 3D cartesian coordinate system constructed specifically to European areas. Hence, we can use
## the coordinates as such. 
##  In whitefish.dat the coordinates are N_etrs89 E_etrs89
##  In whitefish.raster the coordinates are X Y                 

# First examine the data by plotting few of the raster layers 

# !!!!!!!!! Correct
# Visualize few environmental covariates
e <- extent(cbind(whitefish.raster$X,whitefish.raster$Y))
r <- raster(e, ncol=length(unique(whitefish.raster$X)), nrow=length(unique(whitefish.raster$Y)))
# Visualize the study area
windows()
par(mfrow=c(1,2))
z <- rasterize(cbind(whitefish.raster$X,whitefish.raster$Y), r, whitefish.raster$DIST20M, fun=mean)
plot(z, xlim=cbind(min(whitefish.raster$X),max(whitefish.raster$X)),
     ylim=cbind(min(whitefish.raster$Y),max(whitefish.raster$Y)), main="Distance to 20m deep water")
# Plot the locations of sampling sites
points(whitefish.dat$E_etrs89,whitefish.dat$N_etrs89)
z <- rasterize(cbind(whitefish.raster$X,whitefish.raster$Y), r, whitefish.raster$ICELAST09, fun=mean)
plot(z, xlim=cbind(min(whitefish.raster$X),max(whitefish.raster$X)),
     ylim=cbind(min(whitefish.raster$Y),max(whitefish.raster$Y)), main="Last ice cover day in 2009")
# Plot the locations of sampling sites
points(whitefish.dat$E_etrs89,whitefish.dat$N_etrs89)


## START the analysis. The environmental covariates that we are interested in are
#   BOTTOMCLS - Bottom type classifiction, a categorical variable with classes
#               0 = not shallow
#               1 = open water
#               2 = other
#               3 = sand
#               4 = sand/mud
#               5 = sand/stone
#   DIS_SAND  - distance to sandy shore, continuous variable
#   FE300ME   - The average fetch (opennes/exposure) over all directions, continuous variable
#   ICELAST09 - The last ice cover date in winter 2009-10, continuous variable
#   RIVERS    - Influence of rivers (~weighted average distance to river mouths), continuous variable
#   SAL910WIN - Winter salinity in 2009-2010, continuous variable

## The end variable is
#   WHISUM - the number of whitefishes caught in sampling occasion. 

## An "off-set" variable / covariate for likelihood
#   VOLUME - the volume of water sampled

# Set up data
s = as.matrix(cbind(whitefish.dat$E_etrs89,whitefish.dat$N_etrs89)) / 1000   # spatial coordinates in km
x = matrix(0,nrow=nrow(s),ncol=12)      # intercept + 6 BOTTOMCLS classes + 5 continues covariates
x[,1] = 1                             # Set the column corresponding to intercept to 1
x[whitefish.dat$BOTTOMCLS==0,2] = 1   # Set the elements corresponding to BOTTOMCLS = 0 to 1
x[whitefish.dat$BOTTOMCLS==1,3] = 1   # Set the elements corresponding to BOTTOMCLS = 1 to 1
x[whitefish.dat$BOTTOMCLS==2,4] = 1   # Set the elements corresponding to BOTTOMCLS = 2 to 1
x[whitefish.dat$BOTTOMCLS==3,5] = 1   # Set the elements corresponding to BOTTOMCLS = 3 to 1
x[whitefish.dat$BOTTOMCLS==4,6] = 1   # Set the elements corresponding to BOTTOMCLS = 4 to 1
x[whitefish.dat$BOTTOMCLS==5,7] = 1   # Set the elements corresponding to BOTTOMCLS = 5 to 1
xcont = as.matrix(cbind(whitefish.dat$DIS_SAND,
                        whitefish.dat$FE300ME,
                        whitefish.dat$ICELAST09,
                        whitefish.dat$RIVERS,
                        whitefish.dat$SAL910WIN))
stdxcont = apply(xcont, 2, sd)
mxcont = apply(xcont, 2, mean)
x[,8:12] = t( apply( t(apply(xcont,1,'-',mxcont)),1,'/',stdxcont) )    # "standardize" the continuous covariates

# End variable
#transform WHISUM to binary: 0 if 0, 1 otherwise
whitefish.dat$WHISUM[whitefish.dat$WHISUM > 0] <- 1
y = whitefish.dat$WHISUM               # number of counted fish larvae
#head(y)
#tail(y)

# Prediction variables
spred = as.matrix(cbind(whitefish.raster$X,whitefish.raster$Y)) / 1000  # spatial coordinates in km
xpred = matrix(0,nrow=nrow(spred),ncol=12)      # intercept + 6 BOTTOMCLS classes + 5 continues covariates
xpred[,1] = 1                             # Set the column corresponding to intercept to 1
xpred[whitefish.raster$BOTTOMCLS==0,2] = 1   # Set the elements corresponding to BOTTOMCLS = 0 to 1
xpred[whitefish.raster$BOTTOMCLS==1,3] = 1   # Set the elements corresponding to BOTTOMCLS = 1 to 1
xpred[whitefish.raster$BOTTOMCLS==2,4] = 1   # Set the elements corresponding to BOTTOMCLS = 2 to 1
xpred[whitefish.raster$BOTTOMCLS==3,5] = 1   # Set the elements corresponding to BOTTOMCLS = 3 to 1
xpred[whitefish.raster$BOTTOMCLS==4,6] = 1   # Set the elements corresponding to BOTTOMCLS = 4 to 1
xpred[whitefish.raster$BOTTOMCLS==5,7] = 1   # Set the elements corresponding to BOTTOMCLS = 5 to 1
xpredcont = as.matrix(cbind(whitefish.raster$DIS_SAND,
                            whitefish.raster$FE300ME,
                            whitefish.raster$ICELAST09,
                            whitefish.raster$RIVERS,
                            whitefish.raster$SAL910WIN))
xpred[,8:12] = t( apply( t(apply(xpredcont,1,'-',mxcont)),1,'/',stdxcont) )    # "standardize" the continuous covariates


######################################################
## START STAN inference
#############################################################

# Test with smaller data first
#seq(from = 1, to = 1, by = ((to - from)/(length.out - 1)),
#length.out = NULL, along.with = NULL, ...)
#by == step; length == # samples

 x = x[seq(1,217,length=60),]
 #s = s[seq(1,217,length=60),]
 y = y[seq(1,217,length=60)]
 #V = V[seq(1,217,length=60)]
# x = x[seq(1,217,length=120),]
# s = s[seq(1,217,length=120),] 
# y = y[seq(1,217,length=120)] 
# V = V[seq(1,217,length=120)] 


whitefish_dat <- list(Dx = ncol(x),
                      N = nrow(x),
                      x = x,
                      y = y)

#set initial values of the chains
#rnorm(#samples, mu,sigma)
set.seed(124)
#init1 <-list(beta=as.vector(rep(0,ncol(x))))
init1 <- list(beta=as.vector(rnorm(ncol(x),0,10)))
init2 <- list(beta=as.vector(rnorm(ncol(x),0,10)))
## initializing with init2 does not produce samples in Stan
## ?? How to choose initial values of the chains?
init3 <- list(beta=as.vector(rnorm(ncol(x),0,10)))
init4 <- list(beta=as.vector(rnorm(ncol(x),0,10)))
inits = list(init1,init2,init3,init4)


# Compile the STAN model
#try it with init2 and init3 separately, chain=1
library(rstudioapi)
fit <- stan(file = "Log_reg_whitefish.stan", data = whitefish_dat, warmup=200,
            iter = 1000, init=inits , control = list(adapt_delta = 0.99), pars=c("beta"))
beta <- as.matrix(fit)

windows()
print(fit)              
#lp_ ->log density up to a constant 
# https://www.jax.org/news-and-insights/jax-blog/2015/october/lp-in-stan-output
#log density = log likelihood of the observations conditioned on the posterior parameters: p(y | p_post). 

summary(fit)            # summary
stan_trace(fit, pars = c("beta"))         # traceplot
#stan_trace(fit, pars = c("beta") )         #traceplotstan_ac(fit,inc_warmup = FALSE, lags = 25)   # autocorrelation between Markov chain samples
#posterior densties and histograms
stan_dens(fit, pars='beta')
stan_hist(fit, pars='beta', bins=30)
#show credible intervals
plot(fit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon", pars='beta')
#pairwise correlations
pairs(fit)

#plot autocorrelation function of all beta parameters & constant term:
acf(beta)



# posterior density of parameters beta
quietgg(h <- stan_hist(fit, pars = "beta"))


# scatter plot
stan_scat(fit, pars = c("beta[1]", "beta[2]"), color = "black", size = 3)

beta_thinned = as.matrix(beta[seq(1,800,4),])
#mean estimates
means_b = colMeans(beta_thinned)
means_b

#

#POSTERIOR OF THETA (LOG(THETA/1-THETA)=BETA*X. 
#POSTERIOR PREDICTIVE CHECKING
#MCMC for prediction: https://medium.com/@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed

#best thing: to use x_test and y-test data. 
#for this toy example I use one point of x = mean(x_train)


#beta_thinned[,13] - coefficent lp_
#take only 12 coefficients alpha,,,beta and 200 samples from x
#sum (alpha*1[i] + ... + beta1*x[i,j]+ ... beta5*x[i,12])  
mean_lin_sums<-rowSums(beta_thinned[,1:12]*x[1:200,])

theta = 1/(1+exp(-mean_lin_sums))
m_theta = mean(theta)
m_theta
#compute quintiles instead of std for 95% conf interval. 
#q_lower = quantile(theta,0.025)
#q_lower
#q_upper = quantile(theta,0.975)
#q_upper

#simulate y_test
y_Rep <- rbinom(200, 1, theta)
y_Rep
#Compute the probability y=0 and y=1 
p_y0 = sum(y_Rep==0)/length(y_Rep)
p_y0
p_y1 = sum(y_Rep==1)/length(y_Rep)
p_y1


acf(beta)

#I do not know how to compute the accuracy of the estimates; need test data. 

#We can compare STAN estimates to estimates obtained in R, 
#outcomeModel <- glm(y ~ ., data = whitefish_dat, family = binomial(link = "logit"))
#summary(outcomeModel)
# Frequentist
#tableone::ShowRegTable(outcomeModel, exp = FALSE)

#I do not understand tha last plots in the Ex.7.2 example. 