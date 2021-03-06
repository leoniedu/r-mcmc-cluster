## EXAMPLE 1 - 8 schools - rjags and snow package
library(snow)
library(rjags)
n.chains <- 5
if (!exists("cl")) cl <- makeCluster(n.chains, type = "SOCK")
model.file <- system.file(package="R2WinBUGS", "model", "schools.txt")
J <- 8.0
n.groups <- 2
y <- c(28.4,7.9,-2.8,6.8,-0.6,0.6,18.0,12.2)
sigma.y <- c(14.9,10.2,16.3,11.0,9.4,11.4,10.4,17.6)
jags.inits <- function(){ 
    list("mu.theta"=rnorm(1)) 
}

system.time(m1 <- jags.model(model.file,data=Hmisc::llist(J,n.groups,y,sigma.y),
                             inits=c(jags.inits(),.RNG.name="base::Marsaglia-Multicarry"),
                             n.adapt=1000*100,n.chains=n.chains))

system.time(m2 <- cluster.jags.model(cl,file=model.file,data=Hmisc::llist(J,n.groups,y,sigma.y),
                                     inits=jags.inits,
                                     n.adapt=1000*100))

system.time(mc1 <- (coda.samples(m1,c("mu.theta","theta"),n.iter=50000,thin=10)))

system.time(mc2 <- (cluster.coda.samples(cl,m2,c("mu.theta","theta"),n.iter=50000,thin=10)))

## Example 2 - MCMCpack functions using snow package
library(MCMCpack)

##MCMCregress example from MCCMpack help files
line   <- list(X = c(-2,-1,0,1,2), Y = c(1,3,3,3,5))
system.time(posterior1  <- lapply(1:2,function(i) MCMCregress(Y~X, data=line,burnin = 10000, mcmc = 100000,thin=100,seed=list(NA,i))))

##ARGUMENTS MUST BE NAMED!
system.time(posterior2 <- clusterMCMC.snow(cl,f=MCMCregress,formula=Y~X,data=line,burnin = 100000, mcmc = 100000,thin=100,n.chains=4))



library(multicore)


summary(posterior1)
summary(posterior2)
