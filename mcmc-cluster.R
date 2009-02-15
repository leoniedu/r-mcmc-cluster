## recommended/required packages: snow,multicore,jags,MCMCpack,R2WinBUGS

## evalute expression
.eval <- function(evaltext,envir=sys.frame()) {
  eval(parse(text=evaltext), envir=envir)
}
## from the snow package (it's not exported)
gets2 <- function (n, v) {
    assign(n, v, env = .GlobalEnv)
    NULL
}
## from the snow package (original looks only in the global environment but we need also from the local)
clusterExport2 <- function (cl, name,x) {
    clusterCall(cl, gets2, name, x)
}
cluster.jags.model <- function(cl,file, data = sys.frame(sys.parent()), inits,n.adapt = 1000,model.name="mnow") {
    ##max of 4 chains because jags currently uses only 4 Pseudo-RNGs
    if (length(cl)>4) warning("max number of chains advised to be 4\n starting chains with repeated RNGs and different seeds, but might endanger independence.")
    rngnames <- rep(c("base::Wichmann-Hill","base::Marsaglia-Multicarry","base::Super-Duper","base::Mersenne-Twister"),length.out=length(cl))
    seeds <- sample(1:10000,length(cl))
    ## loads rjags
    clusterEvalQ(cl,library(rjags))
    ## export RNGs for each chain
    for (i in 1:length(cl)) {
        inits.i <- c(inits(),.RNG.name=rngnames[i],.RNG.seed=seeds[i])
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

##This is for MCMCpack functions (not jags) using snow package
clusterMCMC.snow <- function (cl,f, ...) {
    require(MCMCpack)
    require(snow)
    ##FIX seed is ignored
    ##FIX put a line removing the new objects from the memory in the slaves?
    n.chains <- length(cl)
    clusterEvalQ(cl,library(MCMCpack))
    fargs <- list(...)
    lapply(names(fargs),function(x) with(fargs,clusterExport3(cl,x,get(x))))
    post <- parLapply(cl,1:n.chains,function(i) f(...,seed=list(NA,i)))
    mcmc.list(post)
}


##This is for MCMCpack functions (not jags) using multicore package
clusterMCMC <- function (f, ...,n.chains=1) {
    require(MCMCpack)
    require(multicore)
    ##FIX seed is ignored
    post <- mclapply(1:n.chains,function(i) {
        f(...,seed=list(NA,i))
    })
    mcmc.list(post)
}

"coda2bugs" <- function(object.list,bugs=TRUE,parms.todrop=NA,parms.tokeep=NA,engine=NULL) {
  ## try to get engine from object.list attribute
  if (is.null(engine)) engine <- attr(object.list,"engine")
  ## if empty assume it is jags
  if (is.null(engine)) engine <- "jags"
  ## One more piece borrowed from bugs.r
  ## This time from the mcsamp function
  ## this function gets a list of mcmc objects
  ## and transforms it to a bugs object (so we can get
  ## the nice bugs.r summary and plots)
  n.chains <- length(object.list)
  parnames <- colnames(object.list[[1]])
  if (!is.na(parms.tokeep[1])) {
    c.tokeep <- NULL
    for (i in 1:length(parms.tokeep)) {
      c.tokeep <- c(c.tokeep,grep(parms.tokeep[i],parnames))
    }
    c.tokeep <- unique(c.tokeep)
  } else {
    c.tokeep <- 1:length(parnames)
  }
  if (n.chains<2) stop ("n.chains must be at least 2")
  first.chain <- object.list[[1]][,c.tokeep]
  n.parameters <- ncol(first.chain)
  n.iter <- nrow(first.chain)
  sims <- array (NA, c(n.iter, n.chains, n.parameters))
  par.names <- dimnames(first.chain)[[2]]
  sims[,1,] <- first.chain
  for (k in 2:n.chains){
    sims[,k,] <- object.list[[k]][,c.tokeep]
  }
  dimnames(sims) <- list (NULL, NULL, par.names)
  if (!is.na(parms.todrop[1])) {
    c.todrop <- NULL
    for (i in 1:length(parms.todrop)) {
      c.todrop <- c(c.todrop,grep(parms.todrop[i],par.names))
    }
    c.todrop <- unique(c.todrop)
    sims <- sims[,,-c.todrop]
    
}
  return (sims)
}



"as.bugs.array2" <- function (sims.array, pts, DIC=FALSE,engine="bugs")
{
  ## 'sims.array' is supposed to be a 3-way array with
  # n.sims*n.chains*n.parameters simulations, and
  # the 3rd component of dimnames(x) should have the parameter names.

  ## From Andrew Gelman's bugs.r function
  ## a couple of lines commented out by Eduardo Leoni (see comment below)
  require("R2WinBUGS")
  d <- dim(sims.array)
  n.burnin     <- 0
  n.keep       <- d[1]
  n.chains     <- d[2]
  n.parameters <- d[3]
  n.sims       <- n.keep*n.chains
  n.iter       <- n.keep
  n.thin       <- 1
  #
  parameter.names <- dimnames(sims.array)[[3]]
  if (is.null(parameter.names)) {
    parameter.names <- paste("P", 1:n.parameters, sep="")
    dimnames(sims.array)[[3]] <- parameter.names
  }
  parameters.to.save <- unique(sapply(strsplit(parameter.names, "\\["), "[", 1))
  #
  sims <- matrix(NA, n.sims, n.parameters)
  root.long <- character(n.parameters)
  indexes.long <- vector(n.parameters, mode = "list")
  for (i in 1:n.parameters) {
    temp <- R2WinBUGS:::decode.parameter.name(parameter.names[i])
    root.long[i] <- temp$root
    indexes.long[[i]] <- temp$indexes
  }
  n.roots <- length(parameters.to.save)
  left.bracket.short <- as.vector(regexpr("[[]", parameters.to.save))
  right.bracket.short <- as.vector(regexpr("[]]", parameters.to.save))
  root.short <- ifelse(left.bracket.short == -1, parameters.to.save,
      substring(parameters.to.save, 1, left.bracket.short -
          1))
  dimension.short <- rep(0, n.roots)
  indexes.short <- vector(n.roots, mode = "list")
  n.indexes.short <- vector(n.roots, mode = "list")
  long.short <- vector(n.roots, mode = "list")
  length.short <- numeric(n.roots)
  for (j in 1:n.roots) {
      long.short[[j]] <- (1:n.parameters)[root.long == root.short[j]]
      length.short[j] <- length(long.short[[j]])
      if (length.short[j] == 0)
          stop(paste("parameter", root.short[[j]], "is not in the model"))
      else if (length.short[j] > 1) {
          dimension.short[j] <- length(indexes.long[[long.short[[j]][1]]])
          n.indexes.short[[j]] <- numeric(dimension.short[j])
          for (k in 1:dimension.short[j]) n.indexes.short[[j]][k] <- length(unique(unlist(lapply(indexes.long[long.short[[j]]],
              .subset, k))))
          length.short[j] <- prod(n.indexes.short[[j]])
          ## Modified by Eduardo Leoni
          ## this check fails if you take out a part of the simulations
          ## (for example, you don't want the array to have some of the
          ## parameters) so I took them out.

          ## if (length(long.short[[j]]) != length.short[j])
          ##   stop(paste("error in parameter", root.short[[j]],
          ##   "in parameters.to.save"))
          indexes.short[[j]] <- as.list(numeric(length.short[j]))
          for (k in 1:length.short[j]) indexes.short[[j]][[k]] <- indexes.long[[long.short[[j]][k]]]
      }
  }
  rank.long <- unlist(long.short)
  # -----
  # yes, it's inefficient to do this, but for now I'm just letting this be as it is:
  for (k in 1:n.parameters) {
    sims[,k] <- as.vector(sims.array[,,k])
  }
  # ----
  dimnames(sims) <- list(NULL, parameter.names)
  summary <- R2WinBUGS:::monitor(sims.array, n.chains, keep.all = TRUE)
  last.values <- as.list(numeric(n.chains))
  for (i in 1:n.chains) {
    n.roots.0 <- if (DIC)
      n.roots - 1
    else n.roots
    last.values[[i]] <- as.list(numeric(n.roots.0))
    names(last.values[[i]]) <- root.short[1:n.roots.0]
    for (j in 1:n.roots.0) {
      if (dimension.short[j] <= 1) {
        last.values[[i]][[j]] <- sims.array[n.keep, i,
                                            long.short[[j]]]
        names(last.values[[i]][[j]]) <- NULL
      }
      else if (engine=="jags") last.values[[i]][[j]] <- array(sims.array[n.keep,
                 i, long.short[[j]]], n.indexes.short[[j]])
      ## only winbugs have to permute the array.
      else last.values[[i]][[j]] <- aperm(array(sims.array[n.keep,
                                                           i, long.short[[j]]], rev(n.indexes.short[[j]])),
                                          dimension.short[j]:1)
    }
  }
  sims <- sims[sample(n.sims), ]
  sims.list <- summary.mean <- summary.sd <- summary.median <- summary.025 <-  summary.975 <- vector(n.roots,
                                                                      mode = "list")
  names(sims.list) <- names(summary.mean) <- names(summary.sd) <- names(summary.median) <- names(summary.025) <- names(summary.975) <- root.short
  for (j in 1:n.roots) {
    if (length.short[j] == 1) {
        sims.list[[j]] <- sims[, long.short[[j]]]
        summary.mean[[j]] <- summary[long.short[[j]], "mean"]
        summary.sd[[j]] <- summary[long.short[[j]], "sd"]
        summary.median[[j]] <- summary[long.short[[j]], "50%"]
        ##ell: added 025 and 975
        summary.025[[j]] <- summary[long.short[[j]], "2.5%"]
        summary.975[[j]] <- summary[long.short[[j]], "97.5%"]
    }
    else if (engine=="bugs") {
      temp2 <- dimension.short[j]:1
      sims.list[[j]] <- aperm(array(sims[, long.short[[j]]],
                                    c(n.sims, rev(n.indexes.short[[j]]))), c(1, (dimension.short[j] +
                                                                                 1):2))
      summary.mean[[j]] <- aperm(array(summary[long.short[[j]],
                                               "mean"], rev(n.indexes.short[[j]])), temp2)
      summary.sd[[j]] <- aperm(array(summary[long.short[[j]],
                                             "sd"], rev(n.indexes.short[[j]])), temp2)
      summary.median[[j]] <- aperm(array(summary[long.short[[j]],
                                                 "50%"], rev(n.indexes.short[[j]])), temp2)
      ##ell: added 025 and 975
      summary.025[[j]] <- aperm(array(summary[long.short[[j]],
                                              "2.5%"], rev(n.indexes.short[[j]])), temp2)
      summary.975[[j]] <- aperm(array(summary[long.short[[j]],
                                              "97.5%"], rev(n.indexes.short[[j]])), temp2)
    } else if (engine=="jags") {
      ##fix this list
      sims.list[[j]] <- sims[, long.short[[j]]]
      summary.mean[[j]] <- array(summary[long.short[[j]],"mean"],n.indexes.short[[j]])
      summary.sd[[j]] <- array(summary[long.short[[j]],"sd"],n.indexes.short[[j]])
      summary.median[[j]] <- array(summary[long.short[[j]],"50%"],n.indexes.short[[j]])
      ##ell: added 025 and 975
      summary.025[[j]] <- array(summary[long.short[[j]],"2.5%"],n.indexes.short[[j]])
      summary.975[[j]] <- array(summary[long.short[[j]],"97.5%"],n.indexes.short[[j]])
    }
  }
  summary <- summary[rank.long, ]
  all <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
        n.thin = n.thin, n.keep = n.keep, n.sims = n.sims, sims.array = sims.array[,
            , rank.long, drop = FALSE], sims.list = sims.list,
        sims.matrix = sims[, rank.long], summary = summary, mean = summary.mean,
        sd = summary.sd, median = summary.median, root.short = root.short,
        long.short = long.short, dimension.short = dimension.short,
        indexes.short = indexes.short, last.values = last.values, is.DIC=DIC,p02.5=summary.025,p97.5=summary.975)
    if (DIC) {
        deviance <- all$sims.array[, , dim(sims.array)[3], drop = FALSE]
        dim(deviance) <- dim(deviance)[1:2]
        pD <- numeric(n.chains)
        DIC <- numeric(n.chains)
        for (i in 1:n.chains) {
            pD[i] <- var(deviance[, i])/2
            DIC[i] <- mean(deviance[, i]) + pD[i]
        }
        all <- c(all, list(pD = mean(pD), DIC = mean(DIC)))
    }
  class(all) <- "bugs"
  return(all)
}
