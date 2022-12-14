#Code for shifted geometric distribution with gamma randomized intercept
#and complementary log-log link function on static covariates
#https://brucehardie.com/notes/037/time-varying_covariates_in_BG.pdf

library(tidyverse)
library(survival)

#####################Simulate data######################

#parameters
N <- 5000
shape <- 2
rate <- 1
b1 <- -.1
b2 <- -.2

#two static covariates
x1 <- rnorm(N, 7, 2)
x2 <- rbinom(N,1,.8)

#gamma randomized intercept
lam <- rgamma(N,shape=shape,rate=rate)

#create individual p with covariates 
p = 1 - exp(-lam*exp(x1*b1+x2*b2))

#simulate from shifted geometric distribution for each p
y <- rgeom(N, p) + 1

#create a censoring point since real application can have censoring
#using uniform censoring here
censor_point=10

#get the id's whose value is above the censoring point
censor_id <-which(y>censor_point)

#censor some rows to mimic real life application
y2 <-y
y2[censor_id] = censor_point

#create a binary variable to indciate censoring
d<-rep(1,N)
d[censor_id] = 0

#create a data frame version to mimic real life application
data_df<-data.frame(time = y2,
                    status = d ,
                    x1=x1,
                    x2=x2)

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
#X: all the covariates in a 2D format
#id: the index of the censored rows
  
  r=par[1];
  alpha=par[2];
  beta=par[-(1:2)];
  
  #uncensored piece of likelihood 
  uncen=y[-id];
  C_u = exp(X[-id,] %*% beta);
  LL_uncen=logdiffexp( -r*log(1+C_u*(uncen-1)/alpha), 
                       -r*log(1+C_u*(uncen)/alpha) );
    
  #censored piece of likelihood  
  cen=y[id];
  C_c = exp(X[id,] %*% beta);
  LL_cen = -r*log(1+C_c*(cen)/alpha);
  
  return (-sum(LL_uncen)-sum(LL_cen));
}

#######Estimation Function Complete##############################################


#######Parameter recovery####################################################
#To use the estimation function to obtain the parameters for the simulated data frame
#we use the standard 'optim' routine from R

solution=optim(par=c(1.5,1.5,rep(0,dim(X)[2]-2)),fn=G2G,
               y = data_df$time,
               X = as.matrix(data_df[,3:dim(data_df)[2]]),
               id = which(data_df$status==0),
               method="L-BFGS-B",
               lower=c(1e-5,1e-5,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf),
               hessian = TRUE)
#standard error and 95% CI
solution$par_stderr<-sqrt(diag(solve(solution$hessian)))
solution$par_upper<-solution$par+1.96*solution$par_stderr
solution$par_lower<-solution$par-1.96*solution$par_stderr
solution$par_true <- c(shape,rate,b1,b2)

#showing the estimated parameters are within the 95% CI
solution$par
solution$par_lower
solution$par_upper
solution$par_true

#visualize the density of the intercept p using gamma parameters
ran <- data.frame(xx = rgamma(10000,solution$par[1],solution$par[2]))
ran$pp <- 1-exp(-ran$xx)
ggplot(ran, aes(x=pp)) + 
  geom_density()


#######Parameter recovery complete##############################################




##########################Following can be ignored##############################



##model survival function
G2G_surv <- function(par,X,t) {
  
  r=par[1];
  alpha=par[2];
  beta=par[-(1:2)];
  
  surv=matrix(0,nrow = dim(X)[1],ncol=length(t));
  
  C_x = exp(X %*% beta);
  
  for (i in (1:length(t)) ) {
    surv[,i] = exp(-r*log(1+C_x*t[i]/alpha));
  }
  
  S=colSums(surv)/dim(X)[1];
  
  return ( data.frame(time=t,
                      surv=S) );
  
}


#It is easier to create X to auto-determine the #covariates
X = as.matrix(data_df[,3:dim(data_df)[2]])

G2G_df<-G2G_surv(solution$par,X,c(seq(0, 20, by = 1)))

#kaplan-meier
fit.km = survfit( Surv(data_df$time, data_df$status) ~ 1, conf.int=F)
km_df <- data.frame(fit.km$time,fit.km$surv)


#Plot them side by side
ggplot() + 
  geom_line(data = G2G_df, aes(x=time, y=surv, color = "Model")) +
  geom_line(data = km_df, aes( x=fit.km.time, y=fit.km.surv, color = "K-M")) + 
  scale_color_manual(name = "", 
                     values = c("Model"="blue", "K-M" = "red"))


###########################################################
#Actual survival curve (without using K-M)
S=rep(0,21);
for (i in 0:20) {
  S[i+1] = sum(y>i)/N;
}

G2G_df$S<-S


#Plot them side by side
ggplot(G2G_df, aes(time)) + 
  geom_line(aes(y = surv, colour = "Model")) + 
  geom_point(aes(y = S, colour = "Actual")) +
  scale_color_manual(name = "", 
                     values = c("Model"="blue", "Actual" = "red"))

