## SIMULATION SET-UP 5:
## Both replication and discovery GWAS available
## Quantitative trait
## Normal effect size distribution
## Significance threshold of alpha=5e-8

################################################################################

##library(devtools)
##devtools::install_github("amandaforde/winnerscurse")
library(winnerscurse)
library(tidyverse)
library(parallel)

set.seed(1998)

## Total number of simulations: 10
tot_sim <- 100
## Fixed total number of SNPs:
n_snps <- 10^6

## Set of scenarios to be tested
sim_params <- expand.grid(
  sim = c(1:tot_sim),
  n_samples = c(30000,300000),
  h2 = c(0.3,0.8),
  prop_effect = c(0.01,0.001),
  S = c(-1, 0, 1)
)

## Run 'useful_funs.R' here in order to define functions required for the
## simulations below.
## NOTE: Simulations currently being run on Windows, hence mc.cores=1.

## Bias Evaluation Metrics:
## 1. Evaluating the fraction of significant SNPs that have been improved due to
## method implementation - effect size estimates have been adjusted so that they
## are closer to the true effect size
## 2. Evaluating the change in average MSE of significant SNPs due to method
## implementation
## 3. Evaluating the relative change in average MSE of significant SNPs due to
## method implementation

################################################################################

## 4a. Replication and discovery GWAS are of equal size

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S,rep_nid=1)
  ss2 <- data.frame(true_beta = ss$true_beta, se = ss$rep_se)
  disc_stats <- simulate_est(ss)
  rep_stats <- simulate_est(ss2)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,alpha=5e-8)

  ## replication - EB:
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,alpha=5e-8)

  ## UMVCUE
  out_UMVCUE <- UMVCUE(disc_stats,rep_stats,alpha=5e-8)
  names(out_UMVCUE)[names(out_UMVCUE) == "beta_disc"] <- "beta"
  names(out_UMVCUE)[names(out_UMVCUE) == "se_disc"] <- "se"
  flb_UMVCUE <- frac_sig_less_bias(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)
  mse_UMVCUE <- mse_sig_improve(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_UMVCUE <- mse_sig_improve_per(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)

  ## Conditional Likelihood 1 - com
  out_cl <- condlike_rep(disc_stats,rep_stats,alpha=5e-8)
  names(out_cl)[names(out_cl) == "beta_disc"] <- "beta"
  names(out_cl)[names(out_cl) == "se_disc"] <- "se"
  flb_com <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=6)
  mse_com <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_com <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=6)

  ## Conditional Likelihood 2 - MLE
  flb_MLE <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=7)
  mse_MLE <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=7)
  rel_mse_MLE <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=7)

  ## Conditional Likelihood 3 - MSE
  flb_MSE <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=8)
  mse_MSE <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=8)
  rel_mse_MSE <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=8)

  ## MSE minimizer - joint
  out_joint <- MSE_minimizer(disc_stats,rep_stats,alpha=5e-8,spline=FALSE)
  names(out_joint)[names(out_joint) == "beta_disc"] <- "beta"
  names(out_joint)[names(out_joint) == "se_disc"] <- "se"
  flb_joint <- frac_sig_less_bias(out_joint,ss$true_beta,alpha=5e-8,i=6)
  mse_joint <- mse_sig_improve(out_joint,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_joint <- mse_sig_improve_per(out_joint,ss$true_beta,alpha=5e-8,i=6)

  ## MSE minimizer - joint spline version 2
  out_joint_sp <- MSE_minimizer(disc_stats,rep_stats,alpha=5e-8,spline=TRUE)
  names(out_joint_sp)[names(out_joint_sp) == "beta_disc"] <- "beta"
  names(out_joint_sp)[names(out_joint_sp) == "se_disc"] <- "se"
  flb_joint_sp <- frac_sig_less_bias(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)
  mse_joint_sp <- mse_sig_improve(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_joint_sp <- mse_sig_improve_per(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)

  ## Empirical Bayes on combined estimator:
  com_stats <- disc_stats
  com_stats$beta <-  ((((rep_stats$se)^2)*(disc_stats$beta))+(((disc_stats$se)^2)*(rep_stats$beta)))/(((disc_stats$se)^2) + ((rep_stats$se)^2))
  com_stats$se <- sqrt((((disc_stats$se)^2)*((rep_stats$se)^2))/(((disc_stats$se)^2) + ((rep_stats$se)^2)))
  out_EB_com <- empirical_bayes(com_stats)
  out_EB_com <- dplyr::arrange(out_EB_com,out_EB_com$rsid)
  out_EB_com <- cbind(disc_stats,out_EB_com$beta_EB)
  flb_EB_com <- frac_sig_less_bias(out_EB_com,ss$true_beta,alpha=5e-8)
  mse_EB_com <- mse_sig_improve(out_EB_com,ss$true_beta,alpha=5e-8)
  rel_mse_EB_com <- mse_sig_improve_per(out_EB_com,ss$true_beta,alpha=5e-8)

  return(c(flb_EB,mse_EB,rel_mse_EB,flb_rep,mse_rep,rel_mse_rep,flb_UMVCUE,mse_UMVCUE,rel_mse_UMVCUE,flb_com,mse_com,rel_mse_com,flb_MLE,mse_MLE,rel_mse_MLE,flb_MSE,mse_MSE,rel_mse_MSE,flb_joint,mse_joint,rel_mse_joint,flb_joint_sp,mse_joint_sp,rel_mse_joint_sp,flb_EB_com,mse_EB_com,rel_mse_EB_com))
}

res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)


################################################################################

## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rel_mse[i] <- res[[i]][3]
}
res_EB <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)


## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][4]
  mse[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_rep <- cbind(sim_params,flb,mse,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## UMVCUE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rel_mse[i] <- res[[i]][9]
}
res_UMVCUE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_UMVCUE <- ave_results(res_UMVCUE,tot_sim)

## Conditional Likelihood 1 - com
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][10]
  mse[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_com <- cbind(sim_params,flb,mse,rel_mse)
ave_res_com <- ave_results(res_com,tot_sim)

## Conditional Likelihood 2 - MLE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rel_mse[i] <- res[[i]][15]
}
res_MLE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_MLE <- ave_results(res_MLE,tot_sim)

## Conditional Likelihood 3 - MSE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][16]
  mse[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_MSE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_MSE <- ave_results(res_MSE,tot_sim)

## MSE minimizer - joint
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rel_mse[i] <- res[[i]][21]
}
res_joint <- cbind(sim_params,flb,mse,rel_mse)
ave_res_joint <- ave_results(res_joint,tot_sim)

## MSE minimizer - joint spline
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][22]
  mse[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_joint_sp <- cbind(sim_params,flb,mse,rel_mse)
ave_res_joint_sp <- ave_results(res_joint_sp,tot_sim)

## Empirical Bayes using combined estimator
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rel_mse[i] <- res[[i]][27]
}
res_EB_com <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB_com <- ave_results(res_EB_com,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_rep,ave_res_UMVCUE,ave_res_com,ave_res_MLE,ave_res_MSE,ave_res_joint,ave_res_joint_sp,ave_res_EB_com)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("rep",(nrow(sim_params)/tot_sim)),rep("UMVCUE",(nrow(sim_params)/tot_sim)),rep("cl1_com",(nrow(sim_params)/tot_sim)),rep("cl2_MLE",(nrow(sim_params)/tot_sim)),rep("cl3_MSE",(nrow(sim_params)/tot_sim)),rep("MSE_min",(nrow(sim_params)/tot_sim)), rep("MSE_min_sp",(nrow(sim_params)/tot_sim)), rep("EB_com",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/replicate_norm_5e-8_100sim.csv")

################################################################################

## 4b. Replication GWAS is half the size of the discovery GWAS

set.seed(1998)
run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S,rep_nid=0.5)
  ss2 <- data.frame(true_beta = ss$true_beta, se = ss$rep_se)
  disc_stats <- simulate_est(ss)
  rep_stats <- simulate_est(ss2)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,alpha=5e-8)

  ## replication - EB:
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,alpha=5e-8)

  ## UMVCUE
  out_UMVCUE <- UMVCUE(disc_stats,rep_stats,alpha=5e-8)
  names(out_UMVCUE)[names(out_UMVCUE) == "beta_disc"] <- "beta"
  names(out_UMVCUE)[names(out_UMVCUE) == "se_disc"] <- "se"
  flb_UMVCUE <- frac_sig_less_bias(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)
  mse_UMVCUE <- mse_sig_improve(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_UMVCUE <- mse_sig_improve_per(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)

  ## Conditional Likelihood 1 - com
  out_cl <- condlike_rep(disc_stats,rep_stats,alpha=5e-8)
  names(out_cl)[names(out_cl) == "beta_disc"] <- "beta"
  names(out_cl)[names(out_cl) == "se_disc"] <- "se"
  flb_com <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=6)
  mse_com <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_com <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=6)

  ## Conditional Likelihood 2 - MLE
  flb_MLE <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=7)
  mse_MLE <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=7)
  rel_mse_MLE <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=7)

  ## Conditional Likelihood 3 - MSE
  flb_MSE <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=8)
  mse_MSE <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=8)
  rel_mse_MSE <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=8)

  ## MSE minimizer - joint
  out_joint <- MSE_minimizer(disc_stats,rep_stats,alpha=5e-8, spline=FALSE)
  names(out_joint)[names(out_joint) == "beta_disc"] <- "beta"
  names(out_joint)[names(out_joint) == "se_disc"] <- "se"
  flb_joint <- frac_sig_less_bias(out_joint,ss$true_beta,alpha=5e-8,i=6)
  mse_joint <- mse_sig_improve(out_joint,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_joint <- mse_sig_improve_per(out_joint,ss$true_beta,alpha=5e-8,i=6)

  ## MSE minimizer - joint spline version 2
  out_joint_sp <- MSE_minimizer(disc_stats,rep_stats,alpha=5e-8,spline=TRUE)
  names(out_joint_sp)[names(out_joint_sp) == "beta_disc"] <- "beta"
  names(out_joint_sp)[names(out_joint_sp) == "se_disc"] <- "se"
  flb_joint_sp <- frac_sig_less_bias(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)
  mse_joint_sp <- mse_sig_improve(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_joint_sp <- mse_sig_improve_per(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)

  ## Empirical Bayes on combined estimator:
  com_stats <- disc_stats
  com_stats$beta <-  ((((rep_stats$se)^2)*(disc_stats$beta))+(((disc_stats$se)^2)*(rep_stats$beta)))/(((disc_stats$se)^2) + ((rep_stats$se)^2))
  com_stats$se <- sqrt((((disc_stats$se)^2)*((rep_stats$se)^2))/(((disc_stats$se)^2) + ((rep_stats$se)^2)))
  out_EB_com <- empirical_bayes(com_stats)
  out_EB_com <- dplyr::arrange(out_EB_com,out_EB_com$rsid)
  out_EB_com <- cbind(disc_stats,out_EB_com$beta_EB)
  flb_EB_com <- frac_sig_less_bias(out_EB_com,ss$true_beta,alpha=5e-8)
  mse_EB_com <- mse_sig_improve(out_EB_com,ss$true_beta,alpha=5e-8)
  rel_mse_EB_com <- mse_sig_improve_per(out_EB_com,ss$true_beta,alpha=5e-8)

  return(c(flb_EB,mse_EB,rel_mse_EB,flb_rep,mse_rep,rel_mse_rep,flb_UMVCUE,mse_UMVCUE,rel_mse_UMVCUE,flb_com,mse_com,rel_mse_com,flb_MLE,mse_MLE,rel_mse_MLE,flb_MSE,mse_MSE,rel_mse_MSE,flb_joint,mse_joint,rel_mse_joint,flb_joint_sp,mse_joint_sp,rel_mse_joint_sp,flb_EB_com,mse_EB_com,rel_mse_EB_com))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)


################################################################################

## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rel_mse[i] <- res[[i]][3]
}
res_EB <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)


## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][4]
  mse[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_rep <- cbind(sim_params,flb,mse,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## UMVCUE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rel_mse[i] <- res[[i]][9]
}
res_UMVCUE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_UMVCUE <- ave_results(res_UMVCUE,tot_sim)

## Conditional Likelihood 1 - com
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][10]
  mse[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_com <- cbind(sim_params,flb,mse,rel_mse)
ave_res_com <- ave_results(res_com,tot_sim)

## Conditional Likelihood 2 - MLE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rel_mse[i] <- res[[i]][15]
}
res_MLE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_MLE <- ave_results(res_MLE,tot_sim)

## Conditional Likelihood 3 - MSE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][16]
  mse[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_MSE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_MSE <- ave_results(res_MSE,tot_sim)

## MSE minimizer - joint
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rel_mse[i] <- res[[i]][21]
}
res_joint <- cbind(sim_params,flb,mse,rel_mse)
ave_res_joint <- ave_results(res_joint,tot_sim)

## MSE minimizer - joint spline
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][22]
  mse[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_joint_sp <- cbind(sim_params,flb,mse,rel_mse)
ave_res_joint_sp <- ave_results(res_joint_sp,tot_sim)

## Empirical Bayes using combined estimator
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rel_mse[i] <- res[[i]][27]
}
res_EB_com <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB_com <- ave_results(res_EB_com,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_rep,ave_res_UMVCUE,ave_res_com,ave_res_MLE,ave_res_MSE,ave_res_joint,ave_res_joint_sp,ave_res_EB_com)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("rep",(nrow(sim_params)/tot_sim)),rep("UMVCUE",(nrow(sim_params)/tot_sim)),rep("cl1_com",(nrow(sim_params)/tot_sim)),rep("cl2_MLE",(nrow(sim_params)/tot_sim)),rep("cl3_MSE",(nrow(sim_params)/tot_sim)),rep("MSE_min",(nrow(sim_params)/tot_sim)), rep("MSE_min_sp",(nrow(sim_params)/tot_sim)), rep("EB_com", (nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/replicate_norm_5e-8_100sim_halfrep.csv")

################################################################################

## 4c. Replication GWAS is 10% the size of the discovery GWAS

set.seed(1998)
run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S,rep_nid=0.1)
  ss2 <- data.frame(true_beta = ss$true_beta, se = ss$rep_se)
  disc_stats <- simulate_est(ss)
  rep_stats <- simulate_est(ss2)


  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,alpha=5e-8)

  ## replication - EB:
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,alpha=5e-8)

  ## UMVCUE
  out_UMVCUE <- UMVCUE(disc_stats,rep_stats,alpha=5e-8)
  names(out_UMVCUE)[names(out_UMVCUE) == "beta_disc"] <- "beta"
  names(out_UMVCUE)[names(out_UMVCUE) == "se_disc"] <- "se"
  flb_UMVCUE <- frac_sig_less_bias(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)
  mse_UMVCUE <- mse_sig_improve(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_UMVCUE <- mse_sig_improve_per(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)

  ## Conditional Likelihood 1 - com
  out_cl <- condlike_rep(disc_stats,rep_stats,alpha=5e-8)
  names(out_cl)[names(out_cl) == "beta_disc"] <- "beta"
  names(out_cl)[names(out_cl) == "se_disc"] <- "se"
  flb_com <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=6)
  mse_com <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_com <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=6)

  ## Conditional Likelihood 2 - MLE
  flb_MLE <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=7)
  mse_MLE <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=7)
  rel_mse_MLE <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=7)

  ## Conditional Likelihood 3 - MSE
  flb_MSE <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=8)
  mse_MSE <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=8)
  rel_mse_MSE <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=8)

  ## MSE minimizer - joint
  out_joint <- MSE_minimizer(disc_stats,rep_stats,alpha=5e-8, spline=FALSE)
  names(out_joint)[names(out_joint) == "beta_disc"] <- "beta"
  names(out_joint)[names(out_joint) == "se_disc"] <- "se"
  flb_joint <- frac_sig_less_bias(out_joint,ss$true_beta,alpha=5e-8,i=6)
  mse_joint <- mse_sig_improve(out_joint,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_joint <- mse_sig_improve_per(out_joint,ss$true_beta,alpha=5e-8,i=6)

  ## MSE minimizer - joint spline version 2
  out_joint_sp <- MSE_minimizer(disc_stats,rep_stats,alpha=5e-8,spline=TRUE)
  names(out_joint_sp)[names(out_joint_sp) == "beta_disc"] <- "beta"
  names(out_joint_sp)[names(out_joint_sp) == "se_disc"] <- "se"
  flb_joint_sp <- frac_sig_less_bias(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)
  mse_joint_sp <- mse_sig_improve(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_joint_sp <- mse_sig_improve_per(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)

  ## Empirical Bayes on combined estimator:
  com_stats <- disc_stats
  com_stats$beta <-  ((((rep_stats$se)^2)*(disc_stats$beta))+(((disc_stats$se)^2)*(rep_stats$beta)))/(((disc_stats$se)^2) + ((rep_stats$se)^2))
  com_stats$se <- sqrt((((disc_stats$se)^2)*((rep_stats$se)^2))/(((disc_stats$se)^2) + ((rep_stats$se)^2)))
  out_EB_com <- empirical_bayes(com_stats)
  out_EB_com <- dplyr::arrange(out_EB_com,out_EB_com$rsid)
  out_EB_com <- cbind(disc_stats,out_EB_com$beta_EB)
  flb_EB_com <- frac_sig_less_bias(out_EB_com,ss$true_beta,alpha=5e-8)
  mse_EB_com <- mse_sig_improve(out_EB_com,ss$true_beta,alpha=5e-8)
  rel_mse_EB_com <- mse_sig_improve_per(out_EB_com,ss$true_beta,alpha=5e-8)

  return(c(flb_EB,mse_EB,rel_mse_EB,flb_rep,mse_rep,rel_mse_rep,flb_UMVCUE,mse_UMVCUE,rel_mse_UMVCUE,flb_com,mse_com,rel_mse_com,flb_MLE,mse_MLE,rel_mse_MLE,flb_MSE,mse_MSE,rel_mse_MSE,flb_joint,mse_joint,rel_mse_joint,flb_joint_sp,mse_joint_sp,rel_mse_joint_sp,flb_EB_com,mse_EB_com,rel_mse_EB_com))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)


################################################################################

## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rel_mse[i] <- res[[i]][3]
}
res_EB <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)


## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][4]
  mse[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_rep <- cbind(sim_params,flb,mse,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## UMVCUE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rel_mse[i] <- res[[i]][9]
}
res_UMVCUE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_UMVCUE <- ave_results(res_UMVCUE,tot_sim)

## Conditional Likelihood 1 - com
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][10]
  mse[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_com <- cbind(sim_params,flb,mse,rel_mse)
ave_res_com <- ave_results(res_com,tot_sim)

## Conditional Likelihood 2 - MLE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rel_mse[i] <- res[[i]][15]
}
res_MLE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_MLE <- ave_results(res_MLE,tot_sim)

## Conditional Likelihood 3 - MSE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][16]
  mse[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_MSE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_MSE <- ave_results(res_MSE,tot_sim)

## MSE minimizer - joint
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rel_mse[i] <- res[[i]][21]
}
res_joint <- cbind(sim_params,flb,mse,rel_mse)
ave_res_joint <- ave_results(res_joint,tot_sim)

## MSE minimizer - joint spline
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][22]
  mse[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_joint_sp <- cbind(sim_params,flb,mse,rel_mse)
ave_res_joint_sp <- ave_results(res_joint_sp,tot_sim)

## Empirical Bayes using combined estimator
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rel_mse[i] <- res[[i]][27]
}
res_EB_com <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB_com <- ave_results(res_EB_com,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_rep,ave_res_UMVCUE,ave_res_com,ave_res_MLE,ave_res_MSE,ave_res_joint,ave_res_joint_sp,ave_res_EB_com)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("rep",(nrow(sim_params)/tot_sim)),rep("UMVCUE",(nrow(sim_params)/tot_sim)),rep("cl1_com",(nrow(sim_params)/tot_sim)),rep("cl2_MLE",(nrow(sim_params)/tot_sim)),rep("cl3_MSE",(nrow(sim_params)/tot_sim)),rep("MSE_min",(nrow(sim_params)/tot_sim)), rep("MSE_min_sp",(nrow(sim_params)/tot_sim)), rep("EB_com", (nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/replicate_norm_5e-8_100sim_10pc.csv")

################################################################################

## 4d. Replication GWAS is twice the size of the discovery GWAS

set.seed(1998)
run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S,rep_nid=2)
  ss2 <- data.frame(true_beta = ss$true_beta, se = ss$rep_se)
  disc_stats <- simulate_est(ss)
  rep_stats <- simulate_est(ss2)


  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,alpha=5e-8)

  ## replication - EB:
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,alpha=5e-8)

  ## UMVCUE
  out_UMVCUE <- UMVCUE(disc_stats,rep_stats,alpha=5e-8)
  names(out_UMVCUE)[names(out_UMVCUE) == "beta_disc"] <- "beta"
  names(out_UMVCUE)[names(out_UMVCUE) == "se_disc"] <- "se"
  flb_UMVCUE <- frac_sig_less_bias(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)
  mse_UMVCUE <- mse_sig_improve(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_UMVCUE <- mse_sig_improve_per(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)

  ## Conditional Likelihood 1 - com
  out_cl <- condlike_rep(disc_stats,rep_stats,alpha=5e-8)
  names(out_cl)[names(out_cl) == "beta_disc"] <- "beta"
  names(out_cl)[names(out_cl) == "se_disc"] <- "se"
  flb_com <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=6)
  mse_com <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_com <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=6)

  ## Conditional Likelihood 2 - MLE
  flb_MLE <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=7)
  mse_MLE <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=7)
  rel_mse_MLE <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=7)

  ## Conditional Likelihood 3 - MSE
  flb_MSE <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=8)
  mse_MSE <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=8)
  rel_mse_MSE <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=8)

  ## MSE minimizer - joint
  out_joint <- MSE_minimizer(disc_stats,rep_stats,alpha=5e-8, spline=FALSE)
  names(out_joint)[names(out_joint) == "beta_disc"] <- "beta"
  names(out_joint)[names(out_joint) == "se_disc"] <- "se"
  flb_joint <- frac_sig_less_bias(out_joint,ss$true_beta,alpha=5e-8,i=6)
  mse_joint <- mse_sig_improve(out_joint,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_joint <- mse_sig_improve_per(out_joint,ss$true_beta,alpha=5e-8,i=6)

  ## MSE minimizer - joint spline version 2
  out_joint_sp <- MSE_minimizer(disc_stats,rep_stats,alpha=5e-8,spline=TRUE)
  names(out_joint_sp)[names(out_joint_sp) == "beta_disc"] <- "beta"
  names(out_joint_sp)[names(out_joint_sp) == "se_disc"] <- "se"
  flb_joint_sp <- frac_sig_less_bias(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)
  mse_joint_sp <- mse_sig_improve(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_joint_sp <- mse_sig_improve_per(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)

  ## Empirical Bayes on combined estimator:
  com_stats <- disc_stats
  com_stats$beta <-  ((((rep_stats$se)^2)*(disc_stats$beta))+(((disc_stats$se)^2)*(rep_stats$beta)))/(((disc_stats$se)^2) + ((rep_stats$se)^2))
  com_stats$se <- sqrt((((disc_stats$se)^2)*((rep_stats$se)^2))/(((disc_stats$se)^2) + ((rep_stats$se)^2)))
  out_EB_com <- empirical_bayes(com_stats)
  out_EB_com <- dplyr::arrange(out_EB_com,out_EB_com$rsid)
  out_EB_com <- cbind(disc_stats,out_EB_com$beta_EB)
  flb_EB_com <- frac_sig_less_bias(out_EB_com,ss$true_beta,alpha=5e-8)
  mse_EB_com <- mse_sig_improve(out_EB_com,ss$true_beta,alpha=5e-8)
  rel_mse_EB_com <- mse_sig_improve_per(out_EB_com,ss$true_beta,alpha=5e-8)

  return(c(flb_EB,mse_EB,rel_mse_EB,flb_rep,mse_rep,rel_mse_rep,flb_UMVCUE,mse_UMVCUE,rel_mse_UMVCUE,flb_com,mse_com,rel_mse_com,flb_MLE,mse_MLE,rel_mse_MLE,flb_MSE,mse_MSE,rel_mse_MSE,flb_joint,mse_joint,rel_mse_joint,flb_joint_sp,mse_joint_sp,rel_mse_joint_sp,flb_EB_com,mse_EB_com,rel_mse_EB_com))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)


################################################################################

## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rel_mse[i] <- res[[i]][3]
}
res_EB <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)


## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][4]
  mse[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_rep <- cbind(sim_params,flb,mse,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## UMVCUE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rel_mse[i] <- res[[i]][9]
}
res_UMVCUE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_UMVCUE <- ave_results(res_UMVCUE,tot_sim)

## Conditional Likelihood 1 - com
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][10]
  mse[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_com <- cbind(sim_params,flb,mse,rel_mse)
ave_res_com <- ave_results(res_com,tot_sim)

## Conditional Likelihood 2 - MLE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rel_mse[i] <- res[[i]][15]
}
res_MLE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_MLE <- ave_results(res_MLE,tot_sim)

## Conditional Likelihood 3 - MSE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][16]
  mse[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_MSE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_MSE <- ave_results(res_MSE,tot_sim)

## MSE minimizer - joint
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rel_mse[i] <- res[[i]][21]
}
res_joint <- cbind(sim_params,flb,mse,rel_mse)
ave_res_joint <- ave_results(res_joint,tot_sim)

## MSE minimizer - joint spline
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][22]
  mse[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_joint_sp <- cbind(sim_params,flb,mse,rel_mse)
ave_res_joint_sp <- ave_results(res_joint_sp,tot_sim)

## Empirical Bayes using combined estimator
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rel_mse[i] <- res[[i]][27]
}
res_EB_com <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB_com <- ave_results(res_EB_com,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_rep,ave_res_UMVCUE,ave_res_com,ave_res_MLE,ave_res_MSE,ave_res_joint,ave_res_joint_sp,ave_res_EB_com)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("rep",(nrow(sim_params)/tot_sim)),rep("UMVCUE",(nrow(sim_params)/tot_sim)),rep("cl1_com",(nrow(sim_params)/tot_sim)),rep("cl2_MLE",(nrow(sim_params)/tot_sim)),rep("cl3_MSE",(nrow(sim_params)/tot_sim)),rep("MSE_min",(nrow(sim_params)/tot_sim)), rep("MSE_min_sp",(nrow(sim_params)/tot_sim)), rep("EB_com", (nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/replicate_norm_5e-8_100sim_2.csv")

################################################################################

## 4e. Replication GWAS is same size of the discovery GWAS and effect size
## distribution is bimodal

set.seed(1998)
run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_bim(H=h2,Pi=prop_effect,nid=n_samples,sc=S,rep_nid=1)
  ss2 <- data.frame(true_beta = ss$true_beta, se = ss$rep_se)
  disc_stats <- simulate_est(ss)
  rep_stats <- simulate_est(ss2)


  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,alpha=5e-8)

  ## replication - EB:
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,alpha=5e-8)

  ## UMVCUE
  out_UMVCUE <- UMVCUE(disc_stats,rep_stats,alpha=5e-8)
  names(out_UMVCUE)[names(out_UMVCUE) == "beta_disc"] <- "beta"
  names(out_UMVCUE)[names(out_UMVCUE) == "se_disc"] <- "se"
  flb_UMVCUE <- frac_sig_less_bias(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)
  mse_UMVCUE <- mse_sig_improve(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_UMVCUE <- mse_sig_improve_per(out_UMVCUE,ss$true_beta,alpha=5e-8,i=6)

  ## Conditional Likelihood 1 - com
  out_cl <- condlike_rep(disc_stats,rep_stats,alpha=5e-8)
  names(out_cl)[names(out_cl) == "beta_disc"] <- "beta"
  names(out_cl)[names(out_cl) == "se_disc"] <- "se"
  flb_com <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=6)
  mse_com <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_com <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=6)

  ## Conditional Likelihood 2 - MLE
  flb_MLE <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=7)
  mse_MLE <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=7)
  rel_mse_MLE <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=7)

  ## Conditional Likelihood 3 - MSE
  flb_MSE <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=8)
  mse_MSE <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=8)
  rel_mse_MSE <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=8)

  ## MSE minimizer - joint
  out_joint <- MSE_minimizer(disc_stats,rep_stats,alpha=5e-8, spline=FALSE)
  names(out_joint)[names(out_joint) == "beta_disc"] <- "beta"
  names(out_joint)[names(out_joint) == "se_disc"] <- "se"
  flb_joint <- frac_sig_less_bias(out_joint,ss$true_beta,alpha=5e-8,i=6)
  mse_joint <- mse_sig_improve(out_joint,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_joint <- mse_sig_improve_per(out_joint,ss$true_beta,alpha=5e-8,i=6)

  ## MSE minimizer - joint spline
  out_joint_sp <- MSE_minimizer(disc_stats,rep_stats,alpha=5e-8,spline=TRUE)
  names(out_joint_sp)[names(out_joint_sp) == "beta_disc"] <- "beta"
  names(out_joint_sp)[names(out_joint_sp) == "se_disc"] <- "se"
  flb_joint_sp <- frac_sig_less_bias(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)
  mse_joint_sp <- mse_sig_improve(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_joint_sp <- mse_sig_improve_per(out_joint_sp,ss$true_beta,alpha=5e-8,i=6)

  ## Empirical Bayes on combined estimator:
  com_stats <- disc_stats
  com_stats$beta <-  ((((rep_stats$se)^2)*(disc_stats$beta))+(((disc_stats$se)^2)*(rep_stats$beta)))/(((disc_stats$se)^2) + ((rep_stats$se)^2))
  com_stats$se <- sqrt((((disc_stats$se)^2)*((rep_stats$se)^2))/(((disc_stats$se)^2) + ((rep_stats$se)^2)))
  out_EB_com <- empirical_bayes(com_stats)
  out_EB_com <- dplyr::arrange(out_EB_com,out_EB_com$rsid)
  out_EB_com <- cbind(disc_stats,out_EB_com$beta_EB)
  flb_EB_com <- frac_sig_less_bias(out_EB_com,ss$true_beta,alpha=5e-8)
  mse_EB_com <- mse_sig_improve(out_EB_com,ss$true_beta,alpha=5e-8)
  rel_mse_EB_com <- mse_sig_improve_per(out_EB_com,ss$true_beta,alpha=5e-8)

  return(c(flb_EB,mse_EB,rel_mse_EB,flb_rep,mse_rep,rel_mse_rep,flb_UMVCUE,mse_UMVCUE,rel_mse_UMVCUE,flb_com,mse_com,rel_mse_com,flb_MLE,mse_MLE,rel_mse_MLE,flb_MSE,mse_MSE,rel_mse_MSE,flb_joint,mse_joint,rel_mse_joint,flb_joint_sp,mse_joint_sp,rel_mse_joint_sp,flb_EB_com,mse_EB_com,rel_mse_EB_com))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)


################################################################################

## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rel_mse[i] <- res[[i]][3]
}
res_EB <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)


## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][4]
  mse[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_rep <- cbind(sim_params,flb,mse,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## UMVCUE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rel_mse[i] <- res[[i]][9]
}
res_UMVCUE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_UMVCUE <- ave_results(res_UMVCUE,tot_sim)

## Conditional Likelihood 1 - com
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][10]
  mse[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_com <- cbind(sim_params,flb,mse,rel_mse)
ave_res_com <- ave_results(res_com,tot_sim)

## Conditional Likelihood 2 - MLE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rel_mse[i] <- res[[i]][15]
}
res_MLE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_MLE <- ave_results(res_MLE,tot_sim)

## Conditional Likelihood 3 - MSE
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][16]
  mse[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_MSE <- cbind(sim_params,flb,mse,rel_mse)
ave_res_MSE <- ave_results(res_MSE,tot_sim)

## MSE minimizer - joint
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rel_mse[i] <- res[[i]][21]
}
res_joint <- cbind(sim_params,flb,mse,rel_mse)
ave_res_joint <- ave_results(res_joint,tot_sim)

## MSE minimizer - joint spline
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][22]
  mse[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_joint_sp <- cbind(sim_params,flb,mse,rel_mse)
ave_res_joint_sp <- ave_results(res_joint_sp,tot_sim)

## Empirical Bayes using combined estimator
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rel_mse[i] <- res[[i]][27]
}
res_EB_com <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB_com <- ave_results(res_EB_com,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_rep,ave_res_UMVCUE,ave_res_com,ave_res_MLE,ave_res_MSE,ave_res_joint,ave_res_joint_sp,ave_res_EB_com)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("rep",(nrow(sim_params)/tot_sim)),rep("UMVCUE",(nrow(sim_params)/tot_sim)),rep("cl1_com",(nrow(sim_params)/tot_sim)),rep("cl2_MLE",(nrow(sim_params)/tot_sim)),rep("cl3_MSE",(nrow(sim_params)/tot_sim)),rep("MSE_min",(nrow(sim_params)/tot_sim)), rep("MSE_min_sp",(nrow(sim_params)/tot_sim)), rep("EB_com", (nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/replicate_bim_5e-8_100sim.csv")

################################################################################

