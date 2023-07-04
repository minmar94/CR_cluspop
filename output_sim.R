# Packages
packs <- c("tidyverse", "magrittr", "pROC", "overlapping", "HDInterval", "Rfast2")
load_packs <- sapply(packs, require, character.only = T)
if(!any(load_packs)) install.packages(packs[!load_packs])
load_packs <- sapply(packs, require, character.only = T)

# Auxiliary functions -----------------------------------------------------
# Function to compute the multiclass AUC
# Inputs: k-th simulated dataset, k-th fitted model
clust.analysis = function(data_sim, fitted_model){
  
  D = data_sim$D
  true_labels_obs = data_sim$clust[data_sim$ind.observed] #true labels for the observed individuals
  
  #estimated labels for the observed individuals
  clust_chains = fitted_model$chains_mat[,paste0("clust[",1:D,"]")]
  
  G = max(clust_chains) #number of mixture components
  
  #relative frequency distribution of the group labels for each individual
  rel_freq_estimated_labels = matrix(unlist(apply(clust_chains,2,function(x) table(factor(x,levels=1:G))/length(x))), ncol=G,byrow=T)
  colnames(rel_freq_estimated_labels) = 1:G
  
  #estimated labels for the observed individuals
  estimated_labels_obs = apply(rel_freq_estimated_labels,1,which.max) 
  
  #confusion matrix for the observed
  conf_matrix_obs = table(true=true_labels_obs, est=estimated_labels_obs) 
  
  #misclassification error
  miscl_error = 1-sum(diag(conf_matrix_obs))/sum(conf_matrix_obs)
  
  #misclassified (observed) individuals 
  miscl_ind = (true_labels_obs!= estimated_labels_obs) #is the observed individual misclassified?
  
  miscl_ind_histories = NA #when no individuals are misclassified
  
  if(sum(miscl_ind)>1){
    miscl_ind_histories = cbind(data_sim$y.obs,true_lab=true_labels_obs,rel_freq_estimated_labels)[miscl_ind,]
    rownames(miscl_ind_histories) = which(miscl_ind)
  }
  
  if(sum(miscl_ind)==1){
    miscl_ind_histories = cbind(data_sim$y.obs,true_lab=true_labels_obs,rel_freq_estimated_labels)[miscl_ind,]
    miscl_ind_histories = t(as.matrix(miscl_ind_histories))
    rownames(miscl_ind_histories) = which(miscl_ind)
  }
  
  #probability to belong to the estimated group 
  post_prob_est_lab = apply(rel_freq_estimated_labels,1,max)
  
  #probability to belong to the estimated group (only for misclassified individuals)
  post_prob_est_lab_obs_miscl = NA
  if(sum(miscl_ind)>0) post_prob_est_lab_obs_miscl = post_prob_est_lab[miscl_ind]
  
  #probability to belong to the estimated group (for observed individuals & correctly classified)
  post_prob_est_lab_obs_correct = NA
  if(sum(!miscl_ind)>0) post_prob_est_lab_obs_correct = post_prob_est_lab[!miscl_ind] 
  
  #multiclass AUC
  multi.AUC = pROC::multiclass.roc(response = true_labels_obs, predictor = rel_freq_estimated_labels)
  
  # output
  return(list(clust_chains=clust_chains,
              rel_freq_estimated_labels=rel_freq_estimated_labels,
              estimated_labels_obs=estimated_labels_obs,
              true_labels_obs=true_labels_obs,
              conf_matrix_obs=conf_matrix_obs,
              miscl_error=miscl_error,
              miscl_ind_histories=miscl_ind_histories,
              post_prob_est_lab=post_prob_est_lab,
              post_prob_est_lab_obs_correct=post_prob_est_lab_obs_correct,
              post_prob_est_lab_obs_miscl=post_prob_est_lab_obs_miscl,
              multi.AUC=multi.AUC))
  
}


# Output analysis ---------------------------------------------------------
# Create object with names of all output files from the simulation study
rdata_list <- list.files("Simulations/", full.names = T)
rdata_list <- rdata_list[c(1:8, 13:44, 9:12)] # reorder: RPT first, than from M1 to M10

# The output files are named following this scheme: rpt_ + fittedmodel_ + T
# where fittedmodel = c("rpt", "mod1", ..., "mod10")
# and T = c(1="T=10", 2="T=20", 3="T=30", 4="T=40") 

# Since it is impossible to load all files without going out-of-memory, the output is analyzed one at a time 
load(rdata_list[1]) # change the index to load a different output

# Name of the model object
mod_name <- ls()[grepl("fit\\_", ls())]
est_mod <- get(mod_name)
rm(list = ls()[grepl("fit\\_", ls())], pos = ".GlobalEnv")

# Name of the simulation object
sim_name <- ls()[grepl("sim\\_", ls())]
sims <- get(sim_name)
rm(list = ls()[grepl("sim\\_", ls())], pos = ".GlobalEnv")

########### RUN ONLY FOR RPT MODELS ###########
# Compute overlapping index (Pastore et al. 2018) for the posterior estimates of the survival probabilities
map_dbl(
  est_mod,
  function(x){
    phi1 = x$chains_mat[,"phi[1]"]
    phi2 = x$chains_mat[,"phi[2]"]
    overlapping::overlap(list(phi1, phi2))$OV
  }
) %>% median # median over the replicas

# Compute mAUCs
mAUCs <- map2(sims, est_mod, function(x, y) clust.analysis(x,y)$multi.AUC)
map_dbl(mAUCs, function(x) as.numeric(x$auc)) %>% median
map_dbl(mAUCs, function(x) as.numeric(x$auc)) %>% range

############  Table 2  ############  
### PHI ###
phi_chains <- map(est_mod, function(x){ x$chains_mat[, grepl("phi", colnames(x$chains_mat))] })

# MAE
map(phi_chains, function(x) t(abs(t(x) - sims[[1]]$phi))) %>% reduce(rbind) %>% Rfast2::colQuantile(probs = .5)

# Coverages
HDIs_phi <- map(phi_chains, function(x) apply(x, 2, HDInterval::hdi))
map_dfr(HDIs_phi, function(x){
  cov1 <- sims[[1]]$phi[1] >= x[1,1] & sims[[1]]$phi[1] <= x[2,1]
  cov2 <- sims[[1]]$phi[2] >= x[1,2] & sims[[1]]$phi[2] <= x[2,2]
  tibble(cov1, cov2)
}) %>% colMeans()

# CIW
map_dfr(HDIs_phi, function(x){
  ciw1 <- x[2,1]-x[1,1]
  ciw2 <- x[2,2]-x[1,2]
  tibble(ciw1, ciw2)
}) %>% colMeans()

### DELTA ###
delta_chains <- map(est_mod, function(x) x$chains_mat[,"delta"])

# MAE
map_dfc(delta_chains, function(x) (abs(x - sims[[1]]$delta))) %>% as.matrix() %>% median

# Coverages
HDIs_delta <- map(delta_chains, HDInterval::hdi)

(Covs <- map_dbl(HDIs_delta, function(x){ x[1] <= sims[[1]]$delta & sims[[1]]$delta <= x[2] }) %>% mean)

# CIW rel
(CIW <- map_dbl(HDIs_delta, function(x){ (x[2]-x[1]) }) %>% median)

### MU ###
mu_chains <- map(est_mod, function(x) x$chains_mat[,"mu"])

# MAE
map_dfc(mu_chains, function(x) (abs(x - sims[[1]]$mu))) %>% as.matrix() %>% median

# Coverages
HDIs_mu <- map(mu_chains, HDInterval::hdi)

(Covs <- map_dbl(HDIs_mu, function(x){ x[1] <= sims[[1]]$mu & sims[[1]]$mu <= x[2] }) %>% mean)

# CIW rel
(CIW <- map_dbl(HDIs_mu, function(x){ (x[2]-x[1]) }) %>% median)

#################### END RUN FOR RPT MODELS ONLY###################################

# Extract Ntrue for each replica
Ntrue <- map_dbl(sims, `[[`, "Nsuper")
# Nsuper estimates: posterior chains
Nsup_chains <- map(est_mod, function(x){ x$chains_mat[, grepl("Nsup", colnames(x$chains_mat))] })

# Relative errors
relerrors_mat <- map2_dfc(Nsup_chains, Ntrue, function(x,y){(x-y)/y})
median(abs(as.matrix(relerrors_mat))) # MAE rel

# Coverage
HDIs <- map(Nsup_chains, HDInterval::hdi)
(Covs <- map2_dbl(HDIs, Ntrue, function(x, y){ x[1] <= y & y <= x[2] }) %>% mean)

# CIW rel
(CIW <- map2_dbl(HDIs, Ntrue, function(x, y){ (x[2]-x[1])/y }) %>% median)

# WAIC
median(map_dbl(est_mod, function(x) x$WAIC$waic))

rm(list = ls()[!grepl("rdata\\_list|clust.analysis", ls())], pos = ".GlobalEnv") # clean and start over from 


# Plots simulations -------------------------------------------------------
# Here, each output is loaded, but only the WAIC and the MAE are kept in two separate matrices (K x n. models)

waic_mat <- MAEs <- list()

for(marco in 1:length(rdata_list)){
  load(rdata_list[marco])
  # Name of the model object
  mod_name <- ls()[grepl("fit\\_", ls())]
  est_mod <- get(mod_name)
  rm(list = ls()[grepl("fit\\_", ls())], pos = ".GlobalEnv")
  
  # Name of the simulation object
  sim_name <- ls()[grepl("sim\\_", ls())]
  sims <- get(sim_name)
  rm(list = ls()[grepl("sim\\_", ls())], pos = ".GlobalEnv")
  
  waic_mat[[marco]] <- map_dbl(est_mod, function(x) x$WAIC$waic)
  Ntrue <- map_dbl(sims, `[[`, "Nsuper")
  # Nsuper estimates
  Nsup_est <- map_dbl(est_mod, function(x){ median(x$chains_mat[, grepl("Nsup", colnames(x$chains_mat))]) })
  MAEs[[marco]] <- -(Ntrue - Nsup_est) # MAE
  rm(est_mod, sims)
  print(marco)
}

######## Figure 2 #######
jpeg("MAE_new.jpg", width = 1200, height = 700, res = 150)
tibble(Model = rep(c("RPT", paste0("M", 1:10)), each = 50*4), Maes = unlist(MAEs), TT = rep(rep(c("T=10","T=20","T=30","T=40"), each = 50), 11)) %>% 
  mutate(Model = factor(Model, levels = c("RPT", paste0("M", 1:10))),
         Ntrue = case_when(TT == "T=10"~170, TT == "T=20"~209, TT == "T=30"~243, TT == "T=40"~271), # expected Nsuper varying T
         Maes = Maes/Ntrue) %>% # compute MAE rel
  ggplot(aes(Model, Maes, fill = factor(TT))) +
  geom_boxplot(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = I("grey")) +
  labs(x = "", y = "Relative absolute error", fill = "") +
  scale_fill_manual(values = c("lightblue", "skyblue3", "blue", "darkblue")) +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 18))
dev.off()

# Model choice
TT10 <- seq(1, 44, 4) # get index of models with same T
app <- reduce(waic_mat[TT10], cbind) %>% as_tibble() %>% set_colnames(value = c("RPT", paste0("M", 1:10)))
(apply(app, 1, function(x) colnames(app)[which.min(x)]) %>% table)/50 # How many times is each model chosen?

TT20 <- seq(2, 44, 4)
app <- reduce(waic_mat[TT20], cbind) %>% as_tibble() %>% set_colnames(value = c("RPT", paste0("M", 1:10)))
(apply(app, 1, function(x) colnames(app)[which.min(x)]) %>% table)/50

TT30 <- seq(3, 44, 4)
app <- reduce(waic_mat[TT30], cbind) %>% as_tibble() %>% set_colnames(value = c("RPT", paste0("M", 1:10)))
(apply(app, 1, function(x) colnames(app)[which.min(x)]) %>% table)/50

TT40 <- seq(4, 44, 4)
app <- reduce(waic_mat[TT40], cbind) %>% as_tibble() %>% set_colnames(value = c("RPT", paste0("M", 1:10)))
(apply(app, 1, function(x) colnames(app)[which.min(x)]) %>% table)/50

