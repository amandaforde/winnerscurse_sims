## SIMULATION SET-UP 4:
## Discovery GWAS only
## Quantitative trait
## Skewed effect size distribution: see simulate_ss_exp() in 'useful_funs.R' for
## more details
## Significance threshold of alpha=5e-8

################################################################################

library(winnerscurse)
library(tidyverse)
library(parallel)
library(scam)

set.seed(1998)

## Total number of simulations: 10
tot_sim <- 20
## Fixed total number of SNPs:
n_snps <- 10^6

## Set of scenarios to be tested
sim_params <- expand.grid(
  sim = c(1:tot_sim),
  n_samples = c(30000,300000),
  h2 = c(0.3,0.8),
  prop_effect = c(0.01, 0.001),
  S = c(-1, 0, 1)
)

## Run 'useful_funs.R' here in order to define functions required for the
## simulations below.
## NOTE: Simulations currently being run on Windows, hence mc.cores=1.


## Bias Evaluation Metrics:
## 1. Evaluating the fraction of significant SNPs that are now less biased due
## to method implementation
## 2. Evaluating the change in average MSE of significant SNPs due to method
## implementation
## 3. Evaluating the relative change in average MSE of significant SNPs due to
## method implementation


run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_exp(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,alpha=5e-8)

  ## Empirical Bayes scam:
  out_EB2 <- empirical_bayes_scam(disc_stats)
  flb_EB2 <- frac_sig_less_bias(out_EB2,ss$true_beta,alpha=5e-8)
  mse_EB2 <- mse_sig_improve(out_EB2,ss$true_beta,alpha=5e-8)
  rel_mse_EB2 <- mse_sig_improve_per(out_EB2,ss$true_beta,alpha=5e-8)

  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,ss$true_beta,alpha=5e-8)
  mse_FIQT <- mse_sig_improve(out_FIQT,ss$true_beta,alpha=5e-8)
  rel_mse_FIQT <- mse_sig_improve_per(out_FIQT,ss$true_beta,alpha=5e-8)

  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,ss$true_beta,alpha=5e-8)
  mse_BR <- mse_sig_improve(out_BR,ss$true_beta,alpha=5e-8)
  rel_mse_BR <- mse_sig_improve_per(out_BR,ss$true_beta,alpha=5e-8)

  ## cl1:
  out_cl <- conditional_likelihood(disc_stats,alpha=5e-8)
  flb_cl1 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8)
  mse_cl1 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8)
  rel_mse_cl1 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8)

  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=5)
  mse_cl2 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=5)
  rel_mse_cl2 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=5)

  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=6)
  mse_cl3 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_cl3 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=6)

  ## replication
  ss2 <- data.frame(true_beta=ss$true_beta,se=ss$rep_se)
  rep_stats <- simulate_est(ss2)
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,alpha=5e-8)

  return(c(flb_EB,mse_EB,rel_mse_EB,flb_EB2,mse_EB2,rel_mse_EB2,flb_FIQT,mse_FIQT,rel_mse_FIQT,flb_BR,mse_BR,rel_mse_BR,flb_cl1,mse_cl1,rel_mse_cl1,flb_cl2,mse_cl2,rel_mse_cl2,flb_cl3,mse_cl3,rel_mse_cl3,flb_rep,mse_rep,rel_mse_rep))
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

## Empirical Bayes scam:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][4]
  mse[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_EB2 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB2 <- ave_results(res_EB2,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rel_mse[i] <- res[[i]][9]
}
res_FIQT <- cbind(sim_params,flb,mse,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## Bootstrap:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][10]
  mse[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_BR <- cbind(sim_params,flb,mse,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)

## Conditional Likelihood 1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rel_mse[i] <- res[[i]][15]
}
res_cl1 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## Conditional Likelihood 2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][16]
  mse[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_cl2 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## Conditional Likelihood 3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rel_mse[i] <- res[[i]][21]
}
res_cl3 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][22]
  mse[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_rep <- cbind(sim_params,flb,mse,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_EB2,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3,ave_res_rep)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("EBscam",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)),rep("rep",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/skew_exp_5e-8_20sim.csv")

################################################################################

