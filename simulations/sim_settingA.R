### Simulating from RPT model

source("fun.R") #load the functions for simulating and fitting

## 1) scenarios description ----

load("simulations/scenariosA.RData") 
scenarios_table_A #it contains details of the 16 scenarios of setting A

#T is the number of primary occasions
#J is the number of secondary occasions within each primary one
#M.sim size of the augmented populations (individuals not recruited in the population or not detected will be discarded)
#p1 and p2 are the detection probabilities (p1<p2)
#phi1 and phi2 are the survival probabilities (phi1<phi2)

T.period = unique(scenarios_table_A[,"T"])

#mixture weights
w = c(0.15, 0.25, 0.6)

#time-varying recruitment parameters (s1: when T=4; s2: when T=8)
rho = list(s1=c(0.3, rep(0.05,T.period[1]-1)),
           s2=c(0.3, rep(0.05,T.period[2]-1))) 

#inflation parameter (for T=4,8)
psi = c(1-prod(1-rho$s1),
        1-prod(1-rho$s2))

#The expected superpopulation abundance, conditional on psi and M, is equal to 300 


## 2) simulation ----

scenario = 1 #insert the scenario you want to simulate from

set.seed(123); seeds = sample(x=1000:9999,size=9000,replace=F)

K = 50 #number of replicas for each scenario
sim = vector("list", K) #it will contain a list, where each element is a different dataset from the selected scenario

for(j in 1:K){
  
  set.seed(seeds[j])
  sim[[j]] = sim.RPT(M.simul = scenarios_table_A[scenario,"M.sim"],
                     T.simul = scenarios_table_A[scenario,"T"], 
                     J.simul = rep(scenarios_table_A[scenario,"J"],scenarios_table_A[scenario,"T"]),
                     rho = rho[[ifelse(scenarios_table_A[scenario,"T"]==4, 1, 2)]],
                     p = scenarios_table_A[scenario,c("p1","p2")],
                     phi = scenarios_table_A[scenario,c("phi1","phi2")],
                     w = w,
                     nzeros = 500) #matrix of observed histories augmented with 500 rows of zeros
  
}


## 3) model fitting ----

#some details of the MCMC fitting via JAGS:

n.chains = 2  #number of parallel chains
n.sample = 1e5  #number of total iterations per chain (including burn in)
n.burnin = 5e4  #length of burn in (i.e. number of iterations to discard at the beginning)
thin.rate = 10 #thinning rate (may be useful to save memory and storage space)
  

#I) Pledger et al. (2003)'s model with 2 class-specific detection probabilities (i.e. p1,p2),
#   constant survival probabilities and mixture weights w=(w1,w2)

fit_JSpmix_G2 = vector("list", K)
set.seed(123)

for(k in 1:K){
  
  fit_JSpmix_G2[[k]] = JSpmix.fit.jags(CR.data.matrix = sim[[k]]$y,
                                       J = sim[[k]]$J,
                                       G = 2,
                                       nc = n.chains,
                                       sample = n.sample,
                                       burnin = n.burnin,
                                       thin = thin.rate,
                                       pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","loglik_i"))
  
}


#II) Pledger et al. (2003)'s model with 2 class-specific couple of parameters, i.e. (p1,phi1) and (p2,phi2),
#    and mixture weights w=(w1,w2)

fit_JSphipmix_G2 = vector("list", K)
set.seed(123)

for(k in 1:K){
  
  fit_JSphipmix_G2[[k]] = JSphipmix.fit.jags(CR.data.matrix = sim[[k]]$y,
                                             J = sim[[k]]$J,
                                             G = 2,
                                             nc = n.chains,
                                             sample = n.sample,
                                             burnin = n.burnin,
                                             thin = thin.rate,
                                             pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","loglik_i"))
  
}


#III) Pledger et al. (2003)'s model with 3 class-specific couple of parameters,  
#     i.e. (p1,phi1), (p2,phi2) and (p3,phi3) and mixture weights w=(w1,w2,w3)

fit_JSphipmix_G3 = vector("list", K)
set.seed(123)

for(k in 1:K){
  
  fit_JSphipmix_G3[[k]] = JSphipmix.fit.jags(CR.data.matrix = sim[[k]]$y,
                                             J = sim[[k]]$J,
                                             G = 3,
                                             nc = n.chains,
                                             sample = n.sample,
                                             burnin = n.burnin,
                                             thin = thin.rate,
                                             pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","loglik_i"))
  
}


#IV) RPT model with mxiture weights w=(w1,w2,w3)

fit_RPT = vector("list", K)

set.seed(123)

for(k in 1:K){
  
  fit_RPT[[k]] = RPT.fit.jags(CR.data.matrix = sim[[k]]$y,
                              J = sim[[k]]$J,
                              nc = n.chains,
                              sample = n.sample,
                              burnin = n.burnin,
                              thin = thin.rate,
                              pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","loglik_i"))
  
}
