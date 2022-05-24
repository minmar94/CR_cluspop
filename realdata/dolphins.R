### REAL DATA ANALYSIS ###

## Preparing the dataset -----

datJS.y = read.csv("realdata/dolphins_4years") #load the data (.csv file)

#primary (on columns) and secondary occasions (as counts for each primary occasion)
#primary occasion: years (T=4 periods); secondary occasion: days of capture in each year

J.y = c(6,15,35,37) #number of secondary occasion (i.e. days of capture) per each primary occasion/period (i.e. season)
T.y = ncol(datJS.y) #number of primary occasions

nzeros = 500 #number of rows of zeros to add to the data matrix (for data augmentation)
datJS.y.aug = rbind(as.matrix(datJS.y),matrix(rep(0,nzeros*T.y),ncol=T.y)) #augmented data matrix


## Fitting the models ----

source("fun.R") #load the script which contains the functions for CR model fitting via JAGS

n.chains = 2  #number of parallel chains
n.sample = 1e5  #number of total iterations per chain (including burn in)
n.burnin = 5e4  #length of burn in (i.e. number of iterations to discard at the beginning)
thin.rate = 2 #thinning rate 

#I) classical JS model (Royle and Dorazio, 2008) with homogeneous detection and survival

set.seed(123)
fit_m1 = JS.fit.jags(CR.data.matrix = datJS.y.aug,
                     J = J.y,
                     nc = n.chains,
                     sample = n.sample,
                     burnin = n.burnin,
                     thin = thin.rate,
                     pars.to.save = c("p","phi","rho","Nsuper","N","clust","z","loglik_i"))

#II) Pledger et al. (2003)'s model with 2 class-specific detection probabilities (i.e. p1,p2),
#   constant survival probabilities and mixture weights w=(w1,w2)

set.seed(123)
fit_m2 = JSpmix.fit.jags(CR.data.matrix = datJS.y.aug,
                         J = J.y,
                         G = 2,
                         nc = n.chains,
                         sample = n.sample,
                         burnin = n.burnin,
                         thin = thin.rate,
                         pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","z","loglik_i"))

#III) Pledger et al. (2003)'s model with 3 class-specific detection probabilities (i.e. p1,p2,p3),
#   constant survival probabilities and mixture weights w=(w1,w2,w3)

set.seed(123)
fit_m3 = JSpmix.fit.jags(CR.data.matrix = datJS.y.aug,
                         J = J.y,
                         G = 3,
                         nc = n.chains,
                         sample = n.sample,
                         burnin = n.burnin,
                         thin = thin.rate,
                         pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","z","loglik_i"))

#IV) Pledger et al. (2003)'s model with 2 class-specific couple of parameters, i.e. (p1,phi1) and (p2,phi2),
#    and mixture weights w=(w1,w2)

set.seed(123)
fit_m4 = JSphipmix.fit.jags(CR.data.matrix = datJS.y.aug,
                            J = J.y,
                            G = 2,
                            nc = n.chains,
                            sample = n.sample,
                            burnin = n.burnin,
                            thin = thin.rate,
                            pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","z","loglik_i"))

#V) Pledger et al. (2003)'s model with 3 class-specific couple of parameters,  
#     i.e. (p1,phi1), (p2,phi2) and (p3,phi3) and mixture weights w=(w1,w2,w3)

set.seed(123)
fit_m5 = JSphipmix.fit.jags(CR.data.matrix = datJS.y.aug,
                            J = J.y,
                            G = 3,
                            nc = n.chains,
                            sample = n.sample,
                            burnin = n.burnin,
                            thin = thin.rate,
                            pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","z","loglik_i"))


#VI) RPT model with mixture weights w=(w1,w2,w3)

set.seed(123)
fit_m6 = RPT.fit.jags(CR.data.matrix = datJS.y.aug,
                      J = J.y,
                      nc = n.chains,
                      sample = n.sample,
                      burnin = n.burnin,
                      thin = thin.rate,
                      pars.to.save = c("w","p","phi","rho","Nsuper","N","clust","z","loglik_i"))
