## SIMULATION SET-UP 1:
## Normal effect size distribution
## Significance threshold of alpha=5e-8

library(winnerscurse)
library(tidyverse)
library(parallel)

set.seed(1998)

## Total number of simulations: 10
tot_sim <- 10
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
## 1. Evaluating the fraction of significant SNPs that are now less biased due
## to method implementation
## 2. Evaluating the change in average MSE of significant SNPs due to method
## implementation
## 3. Evaluating the relative change in average MSE of significant SNPs due to
## method implementation


## Empirical Bayes:
run_sim_EB <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  out <- empirical_bayes(disc_stats)
  flb <- frac_sig_less_bias(out,ss$true_beta,alpha=5e-8)
  mse <- mse_sig_improve(out,ss$true_beta,alpha=5e-8)
  rel_mse <- mse_sig_improve_per(out,ss$true_beta,alpha=5e-8)
  return(c(flb,mse,rel_mse))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim_EB, args=as.list(sim_params[i,]))}, mc.cores=1)

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


## FIQT:
run_sim_FIQT <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  out <- FDR_IQT(disc_stats)
  flb <- frac_sig_less_bias(out,ss$true_beta,alpha=5e-8)
  mse <- mse_sig_improve(out,ss$true_beta,alpha=5e-8)
  rel_mse <- mse_sig_improve_per(out,ss$true_beta,alpha=5e-8)
  return(c(flb,mse,rel_mse))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim_FIQT, args=as.list(sim_params[i,]))}, mc.cores=1)

flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rel_mse[i] <- res[[i]][3]
}
res_FIQT <- cbind(sim_params,flb,mse,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)


## Bootstrap:
run_sim_BR <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  out <- BR_ss(disc_stats)
  flb <- frac_sig_less_bias(out,ss$true_beta,alpha=5e-8)
  mse <- mse_sig_improve(out,ss$true_beta,alpha=5e-8)
  rel_mse <- mse_sig_improve_per(out,ss$true_beta,alpha=5e-8)
  return(c(flb,mse,rel_mse))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim_BR, args=as.list(sim_params[i,]))}, mc.cores=1)

flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rel_mse[i] <- res[[i]][3]
}
res_BR <- cbind(sim_params,flb,mse,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)


## Conditional Likelihood 1:
run_sim_cl1 <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  out <- conditional_likelihood(disc_stats,alpha=5e-8)
  flb <- frac_sig_less_bias(out,ss$true_beta,alpha=5e-8)
  mse <- mse_sig_improve(out,ss$true_beta,alpha=5e-8)
  rel_mse <- mse_sig_improve_per(out,ss$true_beta,alpha=5e-8)
  return(c(flb,mse,rel_mse))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim_cl1, args=as.list(sim_params[i,]))}, mc.cores=1)

flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rel_mse[i] <- res[[i]][3]
}
res_cl1 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)


## Conditional Likelihood 2:
run_sim_cl2 <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  out <- conditional_likelihood(disc_stats,alpha=5e-8)
  flb <- frac_sig_less_bias(out,ss$true_beta,alpha=5e-8,i=5)
  mse <- mse_sig_improve(out,ss$true_beta,alpha=5e-8,i=5)
  rel_mse <- mse_sig_improve_per(out,ss$true_beta,alpha=5e-8,i=5)
  return(c(flb,mse,rel_mse))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim_cl2, args=as.list(sim_params[i,]))}, mc.cores=1)

flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rel_mse[i] <- res[[i]][3]
}
res_cl2 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)


## Conditional Likelihood 3:
run_sim_cl3 <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  out <- conditional_likelihood(disc_stats,alpha=5e-8)
  flb <- frac_sig_less_bias(out,ss$true_beta,alpha=5e-8,i=6)
  mse <- mse_sig_improve(out,ss$true_beta,alpha=5e-8,i=6)
  rel_mse <- mse_sig_improve_per(out,ss$true_beta,alpha=5e-8,i=6)
  return(c(flb,mse,rel_mse))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim_cl3, args=as.list(sim_params[i,]))}, mc.cores=1)

flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rel_mse[i] <- res[[i]][3]
}
res_cl3 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)


## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/norm_5e-8_10sim.csv")

################################################################################


