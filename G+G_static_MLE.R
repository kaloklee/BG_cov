#https://brucehardie.com/notes/037/time-varying_covariates_in_BG.pdf

library(tidyverse)
library(survival)

#Simulate from shifted geometric distribution with gamma randomized intercept
#and complementary log-log link function on (static) covariates

N <- 5000

x1 <- rnorm(N, 7, 2)

x2 <- rbinom(N,1,.8)

lam <- rgamma(N,shape=2,rate=1)

p = 1 - exp(-lam*exp(x1*(-.1)+x2*(-.2)))

X = as.matrix(cbind(x1,x2))

y <- rgeom(N, p) + 1

censor_point=10
censor_id <-which(y>censor_point)

y2 <-y
y2[censor_id] = censor_point
d<-rep(1,N)
d[censor_id] = 0

#create a data frame version
data_df<-data.frame(time = y2,
                    status = d ,
                    x1=x1,
                    x2=x2)

#Create a logdiffexp function to avoid computation error
logdiffexp <- function (a,b) {
  
  c = pmax(a,b);
  return (c + log(exp(a-c)-exp(b-c))) ;
  
}

#model likelihood function
G2G <- function(par,y,X,id) {
  
  r=par[1];
  alpha=par[2];
  beta=par[-(1:2)];
  
  
  #uncen 
  uncen=y[-id];
  C_u = exp(X[-id,] %*% beta);
  LL_uncen=logdiffexp( -r*log(1+C_u*(uncen-1)/alpha), 
                 -r*log(1+C_u*(uncen)/alpha) );
    
  #cen  
  cen=y[id];
  C_c = exp(X[id,] %*% beta);
  LL_cen = -r*log(1+C_c*(cen)/alpha);
  
  return (-sum(LL_uncen)-sum(LL_cen));
}

#If someone is using a real data_df and not from simulation, it's best to create X
X = as.matrix(data_df[,3:dim(data_df)[2]])

solution=optim(par=c(1.5,1.5,rep(0,dim(X)[2])),fn=G2G,
               #y=y2, 
               #X=X,  
               #id=censor_id, 
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

solution$par

#visualize the density of the intercept p
ran <- data.frame(xx = rgamma(10000,solution$par[1],solution$par[2]))
#ggplot(ran, aes(x=xx)) + 
#  geom_density()
ran$pp <- 1-exp(-ran$xx)
ggplot(ran, aes(x=pp)) + 
  geom_density()


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

