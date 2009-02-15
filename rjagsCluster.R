library(snow)
## evalute expression
.eval <- function(evaltext,envir=sys.frame()) {
  eval(parse(text=evaltext), envir=envir)
}
## from the snow package (it's not exported)
gets2 <- function (n, v) {
    assign(n, v, env = .GlobalEnv)
    NULL
}
## from the snow package (looks only in the global environment but we need this from the local)
clusterExport2 <- function (cl, name,x) {
    clusterCall(cl, gets2, name, x)
}
cluster.jags.model <- function(cl,file, data = sys.frame(sys.parent()), inits,n.adapt = 1000,model.name="mnow") {
    ##max of 4 chains because jags currently uses only 4 Pseudo-RNGs
    if (length(cl)>4) stop("max number of chains is 4")
    rngnames <- rep(c("base::Wichmann-Hill","base::Marsaglia-Multicarry","base::Super-Duper","base::Mersenne-Twister"),length.out=length(cl))
    ## loads rjags
    clusterEvalQ(cl,library(rjags))
    ## export RNGs for each chain
    for (i in 1:length(cl)) {
        inits.i <- c(inits(),.RNG.name=rngnames[i])
        clusterExport2(cl[i],"inits",get("inits.i"))        
    }
    ## expression to evaluate at each slave
    tev <- paste(model.name," <- jags.model(file, data, inits, nchain=1,n.adapt)")
    ## exports arguments
    for ( i in c("file","data","nchain","n.adapt","model.name","tev")) {
        clusterExport2(cl,i,get(i))        
    }
    ## evalute expression
    res <- clusterEvalQ(cl,.eval(tev))
    ## saves the model name
    attr(res,"model.name") <- model.name
    ## answer
    res
}
##get coda (mcmc) samples from a model
cluster.coda.samples <- function (cl,model,variable.names = NULL, n.iter, thin = 1) {
    ## search for model name in the model attributes
    tev <- paste("coda.samples(",attr(model,"model.name"),",variable.names,n.iter,thin)")
    for ( i in c("variable.names","n.iter","thin","tev")) {
        clusterExport2(cl,i,get(i))        
    }
    res <- clusterEvalQ(cl,.eval(tev))
    mcmc.list(sapply(res,c))
}
## get jags formated samples from the model
cluster.jags.samples <- function (cl,model,variable.names = NULL, n.iter, thin = 1,type="trace",model.name="mnow") {
    tev <- paste("jags.samples(",attr(model,"model.name"),",variable.names,n.iter,thin,type")
    print(tev)
    for ( i in c("variable.names","n.iter","thin","type","tev")) {
        clusterExport2(cl,i,get(i))        
    }
    res <- clusterEvalQ(cl,.eval(tev))
    mcmc.list(sapply(res,c))
}

