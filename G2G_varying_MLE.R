#Code for shifted geometric distribution with gamma randomized intercept
#and complementary log-log link function on time-varying covariates
#https://brucehardie.com/notes/037/time-varying_covariates_in_BG.pdf

library(tidyverse)
library(survival)

#####################Simulate data######################

#parameters
N <- 5000
shape=2
rate=3
b1 <- (-.1)
b2 <- (.1)

#To simulate time-varying covariates
#Use a reasonably long enough M for the number of periods
M <- 50

x1 <- matrix( rnorm(N*M,mean=0,sd=.5), N, M) 
x2 <- matrix( rbinom(N*M,1,.6), N, M)

#gamma randomized intercept
lam <- rgamma(N,shape=shape,rate=rate)

#create individual, time-varying p with covariates, p is a NxM matrix
p = 1 - exp(-lam*exp(x1*(b1)+x2*(b2)))

#Simulation cannot use geometric because p is changing 
#Bernoulli is used instead
#z is a binary NxM matrix  
z <- matrix( rbinom(length(p) ,1, p), nrow=dim(p)[1])

#Look for the first 1 in z for each row to know when the event occurs
y=rep(0,N)
for (i in 1:N) {
  y[i] <- min(which(z[i,]==1));
}

#Check to see if y is reasonable.
#It is ok to see Inf if the event has not occurred within M periods
summary(y)

#create a censoring point since real application can have censoring
#using uniform censoring here, the only requirement is that it is < M
censor_point=10

#get the id's whose value is above above the censoring point
censor_id <-which(y>censor_point)

#censor some rows to mimic real life application
y2 <-y
y2[censor_id] = censor_point

#create a binary variable to indicate censoring
d<-rep(1,N)
d[censor_id] = 0


#create a data frame version to mimic real life application
#But there is no covariate in this data frame
data_df<-data.frame(time = y2,
                    status = d)
#Instead use a 3D array for K covariates, NxMxK  
#a 2D dataframe is too cumbersome
X = array(c(x1,x2),dim=c(dim(x1)[1],dim(x1)[2],2))

###########Simulation Complete###############################

#################Estimation function#################


#Create a logdiffexp function to avoid computation error
logdiffexp <- function (a,b) {
  
  c = pmax(a,b);
  return (c + log(exp(a-c)-exp(b-c))) ;
  
}

#model likelihood function
G2G <- function(par,y,X,id) {
#par: parameters
#y: duration for each person (uncensored and censored)
#X: all the covariates in a 3D array
#id: the index of the censored rows
  
  r=par[1];
  alpha=par[2];
  coeff=par[-(1:2)];
  
  
  #uncensored piece of likelihood 

  uncen=y[-id];
  X_uncen=X[-id, , ];
  LL_uncen = 0;
  C_u=matrix(0,nrow=length(uncen),ncol=max(uncen));

  for (t in 1:max(uncen)) {

    C_u[, t] = exp( X_uncen[,t,] %*% coeff );
    
    index_t = which(uncen==t);
    
    if (t==1) {
      
      LL_uncen = sum( logdiffexp(0,-r*log(1+C_u[index_t,1]/alpha) ) );
      
    }
    
    else if (t==2) {
      
      LL_uncen= LL_uncen + sum (logdiffexp(-r*log(1+(C_u[index_t,1])/alpha), 
                                           -r*log(1+rowSums(C_u[index_t,1:2])/alpha) ) 
      ); 
      
    }
    
    else {
      
      LL_uncen= LL_uncen + sum (logdiffexp(-r*log(1+rowSums(C_u[index_t,1:t-1])/alpha), 
                                           -r*log(1+rowSums(C_u[index_t,1:t])/alpha) ) 
                              ); 
    }
  }
  
  #censored piece of likelihood  
  cen=y[id];
  X_cen = X[id, ,];
  LL_cen = 0;
  C_c=matrix(0,nrow=length(cen),ncol=max(cen));
  
  for (t in 1:max(cen)) {
    
    C_c[, t] = exp( X_cen[,t,] %*% coeff);
    
    index_t = which(cen==t);
    
    if (length(index_t) > 0) {
      
        if (t == 1) {
          LL_cen = LL_cen + sum ( -r*log(1+(C_c[index_t,1]/alpha) ));
        }
        else { 
          LL_cen= LL_cen + sum (-r*log(1+rowSums(C_c[index_t,1:t]/alpha) )); 
        }
    }
  }  
  
  return ( -(LL_uncen)-(LL_cen) );

}


#######Estimation Function Complete##############################################


#######Parameter recovery####################################################
#To use the estimation function to obtain the parameters for the simulated data frame
#we use the standard 'optim' routine from R


solution=optim(par=c(1,2,rep(0,dim(X)[3])),fn=G2G,
               y = data_df$time,
               id = which(data_df$status==0),
               X=X,
               method="L-BFGS-B",
               lower=c(1e-5,1e-5,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf),
               hessian = TRUE)
#standard error and 95% CI
solution$par_stderr<-sqrt(diag(solve(solution$hessian)))
solution$par_upper<-solution$par+1.96*solution$par_stderr
solution$par_lower<-solution$par-1.96*solution$par_stderr
solution$par_true<-c(shape,rate,b1,b2)

#showing the estimated parameters are within the 95% CI
solution$par
solution$par_upper
solution$par_lower
solution$par_true

#visualize the density of the intercept p
ran <- data.frame(xx = rgamma(10000,solution$par[1],solution$par[2]))
ran$pp <- 1-exp(-ran$xx)
ggplot(ran, aes(x=pp)) + 
  geom_density()

