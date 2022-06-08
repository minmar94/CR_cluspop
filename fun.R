#This script will contain the main functions for simulating CR histories and fitting the candidate models

## SIMULATION ----

#1) function for simulating from the Pledger et al. (2003)'s model
#   with G class-specific detection probabilities and constant survival;
#   p and w are vectors of G elements; rho is a vector of T.simul elements.

sim.JSpmix = function(M.simul,T.simul,J.simul,G,rho,p,phi,w,nzeros){
  
  if(is.unsorted(p)){  #detection probabilities sorted by increasing order
    p=sort(p)
  }
  
  ##latent variables
  
  clust=sample(1:G,size = M.simul,replace = T,prob = w) #class labels
  
  r=matrix(NA,nrow = M.simul,ncol = T.simul)   
  z=r                                   
  
  r[,1]=rep(1,M.simul)  #individual recruitability (binary variable)
  z[,1]=rbinom(M.simul,1,rho[1])  #individual belonging to the population (binary variable)
  
  for(i in 1:M.simul){
    for(t in 2:T.simul){
      r[i,t]=min(r[i,t-1],1-z[i,t-1])
      z[i,t]=rbinom(1,1,phi*z[i,t-1]+rho[t]*r[i,t])
    }
  }
  
  ##derived parameters
  
  N=colSums(z) #population size at time t
  
  ind=ifelse(rowSums(z)>0,T,F)  #individual belonging to the superpopulation (binary variable)
  
  Nsuper=sum(ind)  #superpopulation size
  
  
  ##augmented data
  
  y.sim=matrix(NA,nrow = M.simul,ncol = T.simul)  #matrix of the simulated data (frequency of detection of individual i at period t)
  
  for (i in 1:M.simul) {
    for (t in 1:T.simul) {
      y.sim[i,t]=rbinom(1,J.simul[t],p[clust[i]]*z[i,t])
    }
  }
  
  ind.observed=rowSums(y.sim)>0  #observed individual (binary variable)
  
  y.obs=y.sim[ind.observed,] #matrix of the observed data 
  D=nrow(y.obs)        #number of distinct observed individuals
  
  y=rbind(y.obs,matrix(rep(0,nzeros*T.simul),ncol=T.simul)) #augmented data matrix
  M=nrow(y) #number of detection histories in the augmented data matrix
  
  return(list(M.simul=M.simul,M=M,T.simul=T.simul,J.simul=J.simul,G=G,
              rho=rho,p=p,phi=phi,w=w,
              z=z,r=r,y=y,y.obs=y.obs,y.sim=y.sim,clust=clust,
              N=N,Nsuper=Nsuper,
              D=D,ind=ind,ind.observed=ind.observed))
  
}

#2) function for simulating from the Pledger et al. (2003)'s model
#   with G class-specific couple of parameters (p1,phi1),...,(pG,phiG) 
#   and weights, w1,...,wG, for the G components of the mixture.
#   Detection and survival probabilities are sorted in increasing order (i.e. p1<...<pG and phi1<...<phiG).

sim.JSphipmix = function(M.simul,T.simul,J.simul,G,rho,p,phi,w,nzeros){
  
  if(is.unsorted(p)){  #detection probabilities sorted by increasing order
    p=sort(p)
  }
  
  if(is.unsorted(phi)){  #survival probabilities sorted by increasing order
    phi=sort(phi)
  }
  
  ##latent variables
  
  clust=sample(1:G,size = M.simul,replace = T,prob = w) #class labels
  
  r=matrix(NA,nrow = M.simul,ncol = T.simul)   
  z=r                                   
  
  r[,1]=rep(1,M.simul)  #individual recruitability (binary variable)
  z[,1]=rbinom(M.simul,1,rho[1])  #individual belonging to the population (binary variable)
  
  for(i in 1:M.simul){
    for(t in 2:T.simul){
      r[i,t]=min(r[i,t-1],1-z[i,t-1])
      z[i,t]=rbinom(1,1,phi[clust[i]]*z[i,t-1]+rho[t]*r[i,t])
    }
  }
  
  ##derived parameters
  
  N=colSums(z) #population size at time t
  
  ind=ifelse(rowSums(z)>0,T,F)  #individual belonging to the superpopulation (binary variable)
  
  Nsuper=sum(ind)  #superpopulation size
  
  
  ##augmented data
  
  y.sim=matrix(NA,nrow = M.simul,ncol = T.simul)  #matrix of the simulated data (frequency of detection of individual i at period t)
  
  for (i in 1:M.simul) {
    for (t in 1:T.simul) {
      y.sim[i,t]=rbinom(1,J.simul[t],p[clust[i]]*z[i,t])
    }
  }
  
  ind.observed=rowSums(y.sim)>0  #observed individual (binary variable)
  
  y.obs=y.sim[ind.observed,] #matrix of the observed data 
  D=nrow(y.obs)        #number of distinct observed individuals
  
  y=rbind(y.obs,matrix(rep(0,nzeros*T.simul),ncol=T.simul)) #augmented data matrix
  M=nrow(y) #number of detection histories in the augmented data matrix
  
  return(list(M.simul=M.simul,M=M,T.simul=T.simul,J.simul=J.simul,G=G,
              rho=rho,p=p,phi=phi,w=w,
              z=z,r=r,y=y,y.obs=y.obs,y.sim=y.sim,clust=clust,
              N=N,Nsuper=Nsuper,
              D=D,ind=ind,ind.observed=ind.observed))
  
}


#3) function for simulating from RPT model:
#   residents have higher detection and survival probabilities;
#   transients have lower detection and survial probabilities;
#   part-times have the same detection as transients and the same survival as residents
#
#   p and phi are vectors of two elements; w is a vector of three elements; rho is a vector of T.simul elements.

sim.RPT = function(M.simul,T.simul,J.simul,rho,p,phi,w,nzeros){
  
  if(is.unsorted(p)){  #detection probabilities sorted by increasing order
    p=sort(p)
  }
  
  if(is.unsorted(phi)){  #survival probabilities sorted by increasing order
    phi=sort(phi)
  }
  
  ##latent variables
  
  clust=sample(1:3,size = M.simul,replace = T,prob = w) #class labels
  
  
  r=matrix(NA,nrow = M.simul,ncol = T.simul)   
  z=r                                   
  
  r[,1]=rep(1,M.simul)  #individual recruitability (binary variable)
  z[,1]=rbinom(M.simul,1,rho[1])  #individual belonging to the population (binary variable)
  
  for(i in 1:M.simul){
    phi.temp = ifelse(clust[i]!=3,phi[2],phi[1])
    for(t in 2:T.simul){
      r[i,t]=min(r[i,t-1],1-z[i,t-1])
      z[i,t]=rbinom(1,1,phi.temp*z[i,t-1]+rho[t]*r[i,t])
    }
  }
  
  ##derived parameters
  
  N=colSums(z) #population size at time t
  
  ind=ifelse(rowSums(z)>0,T,F)  #individual belonging to the superpopulation (binary variable)
  
  Nsuper=sum(ind)  #superpopulation size
  
  
  ##augmented data
  
  y.sim=matrix(NA,nrow = M.simul,ncol = T.simul)  #matrix of the simulated data (frequency of detection of individual i at period t)
  
  for (i in 1:M.simul) {
    p.temp = ifelse(clust[i]==1,p[2],p[1])
    for (t in 1:T.simul) {
      y.sim[i,t]=rbinom(1,J.simul[t],p.temp*z[i,t])
    }
  }
  
  ind.observed=rowSums(y.sim)>0  #individual observed (binary variable)
  
  y.obs=y.sim[ind.observed,] #matrix of the observed data 
  D=nrow(y.obs)        #number of distinct individuals observed
  
  y=rbind(y.obs,matrix(rep(0,nzeros*T.simul),ncol=T.simul)) #augmented data matrix
  M=nrow(y)
  
  return(list(M.simul=M.simul,M=M,T.simul=T.simul,J.simul=J.simul,
              rho=rho,p=p,phi=phi,w=w,
              z=z,r=r,y=y,y.obs=y.obs,y.sim=y.sim,clust=clust,
              N=N,Nsuper=Nsuper,
              D=D,ind=ind,ind.observed=ind.observed))
  
}


## MODEL FITTING VIA JAGS ----

#The following functions return a matrix containing binded parameter chains on each column,
#along with some details about the MCMC run

#1) function to fit the classical Royle & Dorazio (2008)'s JS model.
#   The corresponding BUGS model is "JSvanilla.txt".

JS.fit.jags=function(CR.data.matrix, #matrix with elements y_it (capture frequency of individual i at period t)
                     J, #vector of the number of secondary sampling occasions for each primary sampling occasion
                     nc,  #number of chains
                     sample, #samples per chain
                     burnin, #burn in per chain
                     thin, #thinning rate   
                     pars.to.save = c("p","phi","rho","Nsuper","N","z","loglik_i"),
                     seed=123 #the seed will be set once, before generating the starting values
){
  
  M=nrow(CR.data.matrix)
  T.occ=ncol(CR.data.matrix)
  
  require(R2jags)
  
  set.seed(seed)
  
  start.values = replicate(nc,
                           simplify=FALSE,
                           list(z=matrix(rep(1,M*T.occ),ncol=T.occ),
                                rho=runif(T.occ,0,1),
                                phi=runif(1,0,1),
                                p=runif(1,0,1)))
  
  input.data = list("y"=CR.data.matrix,
                    "T"=T.occ,
                    "M"=M,
                    "J"=J) 
  
  
  start.time=Sys.time()
  
  out = jags(data = input.data,
             inits = start.values,
             parameters.to.save = pars.to.save,
             model.file = "BUGS models/JSvanilla.txt",
             n.chains = nc,
             n.iter = sample,
             n.burnin = burnin,
             n.thin = thin)
  
  
  #Let's add the information about the WAIC:
    
  loglik = out$BUGSoutput$sims.list$loglik_i
  #columns are data points, rows are samples
    
  waic = loo::waic(loglik)
  
  
  #let's simplify the jags output in order to reduce the size of the RData file (excluding the chains of 'loglik')
  
  chains_mat = out$BUGSoutput$sims.matrix[,!colnames(out$BUGSoutput$sims.matrix) %in% paste0("loglik_i[",1:M,"]")]
  
  end.time=Sys.time()
  print(end.time-start.time)
  
  return(list(chains_mat=chains_mat,
              start.values=start.values,
              WAIC=waic,
              model=out$model,
              nchains=out$BUGSoutput$n.chains,
              niter=out$BUGSoutput$n.iter,
              nburnin=out$BUGSoutput$n.burnin,
              thin=out$BUGSoutput$n.thin,
              summary=out$BUGSoutput$summary,
              elapsed.time=end.time-start.time))
  
}


#2) function for fitting the Pledger et al. (2003)'s model
#   with G class-specific detection probabilities, constant survival
#   and weights, w1,...,wG, for the G components of the mixture.
#   Detection probabilities are sorted in increasing order (i.e. p1<...<pG).
#   The corresponding BUGS model is "JSpmix.txt". 


JSpmix.fit.jags=function(CR.data.matrix, #matrix with elements y_it (capture frequency of individual i at period t)
                         J, #vector of the number of secondary sampling occasions for each primary sampling occasion
                         G, #number of classes
                         nc,  #number of chains
                         sample, #samples per chain
                         burnin, #burn in per chain
                         thin, #thinning rate
                         pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","z","loglik_i"),
                         seed=123 #the seed will be set once, before generating the starting values
){
  
  M=nrow(CR.data.matrix)
  T.occ=ncol(CR.data.matrix)
  
  require(R2jags)
  
  set.seed(seed)

  start.values = replicate(nc,
                           simplify=FALSE,
                           list(z=matrix(rep(1,M*T.occ),ncol=T.occ),
                                rho=runif(T.occ,0,1),
                                phi=runif(1,0,1),
                                p=sort(runif(G,0,1)),
                                w=as.vector(MCMCpack::rdirichlet(1,rep(1,G)))))
  
  input.data = list("y"=CR.data.matrix,
                    "T"=T.occ,
                    "M"=M,
                    "J"=J,
                    "G"=G) 
  
  
  start.time=Sys.time()
  
  out = jags(data = input.data,
             inits = start.values,
             parameters.to.save = pars.to.save,
             model.file = "BUGS models/JSpmix.txt",
             n.chains = nc,
             n.iter = sample,
             n.burnin = burnin,
             n.thin = thin)
  
  
  #Let's add the information about the WAIC:
  
  loglik = out$BUGSoutput$sims.list$loglik_i
  #columns are data points, rows are samples
  
  waic = loo::waic(loglik)
  
  
  #let's simplify the jags output in order to reduce the size of the RData file (excluding the chains of 'loglik')
  
  chains_mat = out$BUGSoutput$sims.matrix[,!colnames(out$BUGSoutput$sims.matrix) %in% paste0("loglik_i[",1:M,"]")]
  
  end.time=Sys.time()
  print(end.time-start.time)
  
  return(list(chains_mat=chains_mat,
              start.values=start.values,
              WAIC=waic,
              model=out$model,
              nchains=out$BUGSoutput$n.chains,
              niter=out$BUGSoutput$n.iter,
              nburnin=out$BUGSoutput$n.burnin,
              thin=out$BUGSoutput$n.thin,
              summary=out$BUGSoutput$summary,
              elapsed.time=end.time-start.time))
  
}


#3) The following function fits the Pledger et al. (2003)'s model
#   with G class-specific couple of parameters (p1,phi1),...,(pG,phiG) 
#   and weights, w1,...,wG, for the G components of the mixture.
#   Detection and survival probabilities are sorted in increasing order (i.e. p1<...<pG and phi1<...<phiG).
#   The corresponding BUGS model is "JSphipmix.txt".

JSphipmix.fit.jags=function(CR.data.matrix, #matrix with elements y_it (capture frequency of individual i at period t)
                            J, #vector of the number of secondary sampling occasions for each primary sampling occasion
                            G, #number of classes
                            nc,  #number of chains
                            sample, #samples per chain
                            burnin, #burn in per chain
                            thin, #thinning rate
                            pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","z","loglik_i"),
                            seed=123 #the seed will be set once, before generating the starting values
){
  
  M=nrow(CR.data.matrix)
  T.occ=ncol(CR.data.matrix)
  
  require(R2jags)
  
  set.seed(seed)
  
  start.values = replicate(nc,
                           simplify=FALSE,
                           list(z=matrix(rep(1,M*T.occ),ncol=T.occ),
                                rho=runif(T.occ,0,1),
                                phi=sort(runif(G,0,1)),
                                p=sort(runif(G,0,1)),
                                w=as.vector(MCMCpack::rdirichlet(1,rep(1,G)))))
  
  input.data = list("y"=CR.data.matrix,
                    "T"=T.occ,
                    "M"=M,
                    "J"=J,
                    "G"=G) 
  
  
  start.time=Sys.time()
  
  out = jags(data = input.data,
             inits = start.values,
             parameters.to.save = pars.to.save,
             model.file = "BUGS models/JSphipmix.txt",
             n.chains = nc,
             n.iter = sample,
             n.burnin = burnin,
             n.thin = thin)
  
  
  #Let's add the information about the WAIC:
  
  loglik = out$BUGSoutput$sims.list$loglik_i
  #columns are data points, rows are samples
  
  waic = loo::waic(loglik)
  
  
  #let's simplify the jags output in order to reduce the size of the RData file (excluding the chains of 'loglik')
  
  chains_mat = out$BUGSoutput$sims.matrix[,!colnames(out$BUGSoutput$sims.matrix) %in% paste0("loglik_i[",1:M,"]")]
  
  end.time=Sys.time()
  print(end.time-start.time)
  
  return(list(chains_mat=chains_mat,
              start.values=start.values,
              WAIC=waic,
              model=out$model,
              nchains=out$BUGSoutput$n.chains,
              niter=out$BUGSoutput$n.iter,
              nburnin=out$BUGSoutput$n.burnin,
              thin=out$BUGSoutput$n.thin,
              summary=out$BUGSoutput$summary,
              elapsed.time=end.time-start.time))
  
}


#4) The following function fits the RPT model:
#   residents are the first class, with mixture weight w1, detection p2 and survival phi2;
#   part-times are the second class, with mixture weight w2, detection p1 and survival phi2;
#   transients are the third class, with mixture weight w3, detection p1 and survival phi1.
#   Detection and survival probabilities are sorted in increasing order (i.e. p1<p2 and phi1<phi2).
#   The corresponding BUGS model is "RPT.txt".

RPT.fit.jags=function(CR.data.matrix, #matrix with elements y_it (capture frequency of individual i at period t)
                      J, #vector of the number of secondary sampling occasions for each primary sampling occasion
                      nc,  #number of chains
                      sample, #samples per chain
                      burnin, #burn in per chain
                      thin, #thinning rate
                      pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","z","loglik_i"),
                      seed=123 #the seed will be set once, before generating the starting values
){
  
  M=nrow(CR.data.matrix)
  T.occ=ncol(CR.data.matrix)
  
  require(R2jags)
  
  set.seed(seed)
  
  start.values = replicate(nc,
                           simplify=FALSE,
                           list(z=matrix(rep(1,M*T.occ),ncol=T.occ),
                                rho=runif(T.occ,0,1),
                                p=sort(runif(2,0,1)),
                                phi=sort(runif(2,0,1)),
                                w=as.vector(MCMCpack::rdirichlet(1,rep(1,3)))))
  
  input.data = list("y"=CR.data.matrix,
                    "T"=T.occ,
                    "M"=M,
                    "J"=J) 
  
  start.time=Sys.time()
  
  out = jags(data = input.data,
             inits = start.values,
             parameters.to.save = pars.to.save,
             model.file = "BUGS models/RPT.txt",
             n.chains = nc,
             n.iter = sample,
             n.burnin = burnin,
             n.thin = thin)
  
  
  #Let's add the information about the WAIC:
  
  loglik = out$BUGSoutput$sims.list$loglik_i
  #columns are data points, rows are samples
  
  waic = loo::waic(loglik)
  
  
  #let's simplify the jags output in order to reduce the size of the RData file (excluding the chains of 'loglik')
  
  chains_mat = out$BUGSoutput$sims.matrix[,!colnames(out$BUGSoutput$sims.matrix) %in% paste0("loglik_i[",1:M,"]")]
  
  end.time=Sys.time()
  print(end.time-start.time)
  
  return(list(chains_mat=chains_mat,
              start.values=start.values,
              WAIC=waic,
              model=out$model,
              nchains=out$BUGSoutput$n.chains,
              niter=out$BUGSoutput$n.iter,
              nburnin=out$BUGSoutput$n.burnin,
              thin=out$BUGSoutput$n.thin,
              summary=out$BUGSoutput$summary,
              elapsed.time=end.time-start.time))
  
}
