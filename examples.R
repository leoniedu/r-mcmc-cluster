library(R2jags)
source("~/projects/cluelessR/rjagsCluster.R")
## An example model file is given in: 

library(snow)
n.chains <- 2
if (!exists("cl")) cl <- makeCluster(n.chains, type = "SOCK")

model.file <- system.file(package="R2jags", "model", "schools.txt") 


J <- 8.0
n.groups <- 2
y <- c(28.4,7.9,-2.8,6.8,-0.6,0.6,18.0,12.2)
sd <- c(14.9,10.2,16.3,11.0,9.4,11.4,10.4,17.6) 
jags.data <- list("y","sd","J","n.groups") 
jags.params <- c("mu","sigma","theta") 
jags.inits <- function(){ 
    list("mu"=rnorm(1)) 
}
system.time(m1 <- jags.model(model.file,data=Hmisc::llist(J,n.groups,y,sd),
                             inits=c(jags.inits(),.RNG.name="base::Marsaglia-Multicarry"),
                             n.adapt=10000,nchain=n.chains))


system.time(m1 <- cluster.jags.model(cl,file=model.file,data=Hmisc::llist(J,n.groups,y,sd),
                                     inits=jags.inits,
                                     n.adapt=1000,model.name="test1"))


system.time(m2 <- cluster.jags.model(cl,file=model.file,data=Hmisc::llist(J,n.groups,y,sd),
                                     inits=jags.inits,
                                     n.adapt=1000,model.name="test2"))



system.time(tmp1 <- (cluster.coda.samples(cl,m1,c("mu","theta"),n.iter=50,thin=1)))

system.time(tmp1 <- (cluster.coda.samples(cl,m1,c("mu","theta"),n.iter=50,thin=1)))




## Let's take a look: 
# data 






J <- 8.0
n.groups <- 2
y <- c(28.4,7.9,-2.8,6.8,-0.6,0.6,18.0,12.2)
y <- matrix(rnorm(J*n.groups,rep(y,n.groups),1),ncol=n.groups)
sd <- c(14.9,10.2,16.3,11.0,9.4,11.4,10.4,17.6) 
jags.data <- list("y","sd","J","n.groups") 
jags.params <- c("mu","sigma","theta") 
jags.inits <- function(){ 
    list("mu"=rnorm(1)) 
} 
##library(snow)
## if (!exists("cl")) cl <- makeCluster(2, type = "SOCK")
## clusterExport(cl,c("lp","llp","pprior","inits","pm","theta","theta.sd","n.states","n.parties","n.groups"))
jagsfit <- jags3(data=jags.data, inits=jags.inits, jags.params, 
                 n.iter=5000, model.file="model2.bug",cluster=cl
                 )


library(rjags)
m1 <- jags.model("model.bug",##data=listunlist(jags.data),inits=jags.inits,
                 nchain=4,n.adapt=1000)



library(snow)
if (!exists("cl")) cl <- makeCluster(2, type = "SOCK")
data1 <- list(y=y,sd=sd,J=J,n.groups=n.groups)
if (!exists("cl")) cl <- makeCluster(2, type = "SOCK")
cluster.jags.model(cl,"model2.bug",data=data1,inits=jags.inits,
                   n.adapt=100)

m.s <- cluster.coda.samples(cl,variable.names=c("theta","mu","sigma"),n.iter=1000)


m1 <- jags.model("model2.bug",data=data1,inits=jags.inits,
                 nchain=4,n.adapt=100)
m2 <- jags.samples(m1,c("theta","mu"),n.iter=100)

clusterExport(cl,c("data1","jags.inits"))

clusterEvalQ(cl,library(rjags))

m2 <- clusterEvalQ(cl,
                   m <- jags.model("model2.bug",data=data1,
                                   inits=jags.inits,
                                   nchain=2,n.adapt=100))
m2 <- clusterEvalQ(cl,
                   m.s <- coda.samples(m,variable.names=c("theta"),n.iter=100))




                  

                   

(jagsfit$mean$theta)[2,1]
jagsfit$summary[3:4,1:2]



print(jagsfit) 
plot(jagsfit) 




                                        #=============# 
# using jags # 
#=============# 
jagsfit <- jags(data=jags.data, inits=jags.inits, jags.params, 
n.iter=5000, model.file=model.file) 
print(jagsfit) 
plot(jagsfit) 
# if the model does not converge, update it! 
jagsfit.upd <- update(jagsfit, n.iter=1000) 
print(jagsfit.upd) 
plot(jagsfit.upd)
jags2bugs 5 
# or auto update it until it converges! 
jagsfit.upd <- autojags(jagsfit) 
# to get DIC 
dic.samples(jagsfit.upd$model, n.iter=1000, type="pD") 
# to pick up the last save session 
# for example, load("RWorkspace.Rdata") 
recompile(jagsfit) 
jagsfit.upd <- update(jagsfit) 

jagsfit <- jags2(data=jags.data, inits=jags.inits, jags.params, 
n.iter=5000, model.file=model.file) 
print(jagsfit) 
plot(jagsfit) 


