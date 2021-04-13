## PRELIMINARY SIMULATIONS:
## Obtain the average number of significant SNPs and the average proportion of
## these SNPs for which their association estimate is more extreme than their
## true effect size over simulations in which at least one significant SNP has
## been detected, for four different simulation set-ups

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
  prop_effect = c(0.01, 0.001),
  S = c(-1, 0, 1)
)

## Run 'useful_funs.R' here in order to define functions required for the
## simulations below.
## NOTE: Simulations currently being run on Windows, hence mc.cores=1.


##############################################################################
## SIMULATION SET-UP 1:
## Normal effect size distribution
## Significance threshold of alpha=5e-8

run_sim <- function(sim,n_samples, h2, prop_effect, S)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  n_sig <- sum(abs(out$beta/out$se) > qnorm(1-(5e-8)/2))
  return(n_sig)
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)
results <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  results[i] <- res[[i]]
}
res <- cbind(sim_params,results)
ave_res <- ave_results1(res,tot_sim)
ave_res <- ave_res %>%
  rename(
    n_sig = results,
    n_sig_sd = res_error
  )
run_sim2 <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(5e-8)/2),]
  if (nrow(snp_sig) == 0){return(0)}
  prop_bias <- sum(abs(snp_sig$beta) > abs(ss$true_beta[snp_sig$rsid]))/nrow(snp_sig)
  return(prop_bias)
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim2, args=as.list(sim_params[i,]))}, mc.cores=1)
results <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  results[i] <- res[[i]]
}
res <- cbind(sim_params,results)
ave_res2 <- ave_results(res,tot_sim)
ave_res <- cbind(ave_res, ave_res2[,5:6])
ave_res <- ave_res %>%
  rename(
    prop_bias = results,
    prop_bias_sd = res_error
  )
write.csv(ave_res, "norm_nsig_prop_bias_5e-8.csv")


##############################################################################
## SIMULATION SET-UP 2:
## Normal effect size distribution
## Significance threshold of alpha=5e-4

run_sim <- function(sim,n_samples, h2, prop_effect, S)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  n_sig <- sum(abs(out$beta/out$se) > qnorm(1-(5e-4)/2))
  return(n_sig)
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)
results <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  results[i] <- res[[i]]
}
res <- cbind(sim_params,results)
ave_res <- ave_results1(res,tot_sim)
ave_res <- ave_res %>%
  rename(
    n_sig = results,
    n_sig_sd = res_error
  )
run_sim2 <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(5e-4)/2),]
  if (nrow(snp_sig) == 0){return(0)}
  prop_bias <- sum(abs(snp_sig$beta) > abs(ss$true_beta[snp_sig$rsid]))/nrow(snp_sig)
  return(prop_bias)
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim2, args=as.list(sim_params[i,]))}, mc.cores=1)
results <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  results[i] <- res[[i]]
}
res <- cbind(sim_params,results)
ave_res2 <- ave_results(res,tot_sim)
ave_res <- cbind(ave_res, ave_res2[,5:6])
ave_res <- ave_res %>%
  rename(
    prop_bias = results,
    prop_bias_sd = res_error
  )
write.csv(ave_res, "norm_nsig_prop_bias_5e-4.csv")

##############################################################################
## SIMULATION SET-UP 3:
## Skewed effect size distribution: when S=0, 50% of effect sizes are generated
## from a N(0,1) distribution while the other 50% are generated from a N(2.5,1)
## distribution - see simulate_ss_skew() in 'useful_funs.R' for more details
## Significance threshold of alpha=5e-8

run_sim <- function(sim,n_samples, h2, prop_effect, S)
{
  ss <- simulate_ss_skew(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  n_sig <- sum(abs(out$beta/out$se) > qnorm(1-(5e-8)/2))
  return(n_sig)
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)
results <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  results[i] <- res[[i]]
}
res <- cbind(sim_params,results)
ave_res <- ave_results1(res,tot_sim)
ave_res <- ave_res %>%
  rename(
    n_sig = results,
    n_sig_sd = res_error
  )
run_sim2 <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_skew(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(5e-8)/2),]
  if (nrow(snp_sig) == 0){return(0)}
  prop_bias <- sum(abs(snp_sig$beta) > abs(ss$true_beta[snp_sig$rsid]))/nrow(snp_sig)
  return(prop_bias)
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim2, args=as.list(sim_params[i,]))}, mc.cores=1)
results <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  results[i] <- res[[i]]
}
res <- cbind(sim_params,results)
ave_res2 <- ave_results(res,tot_sim)
ave_res <- cbind(ave_res, ave_res2[,5:6])
ave_res <- ave_res %>%
  rename(
    prop_bias = results,
    prop_bias_sd = res_error
  )
write.csv(ave_res, "skew_nsig_prop_bias_5e-8.csv")


##############################################################################
## SIMULATION SET-UP 4:
## Skewed effect size distribution: when S=0, 50% of effect sizes are generated
## from a N(0,1) distribution while the other 50% are generated from a N(2.5,1)
## distribution - see simulate_ss_skew() in 'useful_funs.R' for more details
## Significance threshold of alpha=5e-4

run_sim <- function(sim,n_samples, h2, prop_effect, S)
{
  ss <- simulate_ss_skew(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  n_sig <- sum(abs(out$beta/out$se) > qnorm(1-(5e-4)/2))
  return(n_sig)
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)
results <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  results[i] <- res[[i]]
}
res <- cbind(sim_params,results)
ave_res <- ave_results1(res,tot_sim)
ave_res <- ave_res %>%
  rename(
    n_sig = results,
    n_sig_sd = res_error
  )
run_sim2 <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_skew(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(5e-4)/2),]
  if (nrow(snp_sig) == 0){return(0)}
  prop_bias <- sum(abs(snp_sig$beta) > abs(ss$true_beta[snp_sig$rsid]))/nrow(snp_sig)
  return(prop_bias)
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim2, args=as.list(sim_params[i,]))}, mc.cores=1)
results <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  results[i] <- res[[i]]
}
res <- cbind(sim_params,results)
ave_res2 <- ave_results(res,tot_sim)
ave_res <- cbind(ave_res, ave_res2[,5:6])
ave_res <- ave_res %>%
  rename(
    prop_bias = results,
    prop_bias_sd = res_error
  )
write.csv(ave_res, "skew_nsig_prop_bias_5e-4.csv")
