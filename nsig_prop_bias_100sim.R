## PRELIMINARY SIMULATIONS:
## Obtain the average number of significant SNPs and the average proportion of
## these SNPs for which their association estimate is more extreme than their
## true effect size over simulations in which at least one significant SNP has
## been detected, for four different simulation set-ups

library(tidyverse)
library(parallel)

set.seed(1998)

## Total number of simulations: 100
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

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(5e-8)/2),]
  n_sig <- nrow(snp_sig)
  if (n_sig == 0){
    prop_bias <- -1
    mse <- -1
    prop_x <- -1
    }else{
    prop_bias <- sum(abs(snp_sig$beta) > abs(ss$true_beta[snp_sig$rsid]))/n_sig
    prop_x <- sum(abs(snp_sig$beta) > (abs(ss$true_beta[snp_sig$rsid]) + 1.96*ss$se[snp_sig$rsid]))/n_sig
    mse <- mean((ss$true_beta[snp_sig$rsid]-snp_sig$beta)^2)
  }
  return(c(n_sig,prop_bias,prop_x,mse))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

n_sig <- c(rep(0,nrow(sim_params)))
prop_bias <- c(rep(0,nrow(sim_params)))
prop_x <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  n_sig[i] <- res[[i]][1]
  prop_bias[i] <- res[[i]][2]
  prop_x[i] <- res[[i]][3]
  mse[i] <- res[[i]][4]
}

results <- cbind(sim_params,n_sig,prop_bias,prop_x,mse)
ave_res <- ave_results1(results,tot_sim)
write.csv(ave_res, "results/norm_nsig_prop_bias_5e-8.csv")


##############################################################################
## SIMULATION SET-UP 2:
## Normal effect size distribution
## Significance threshold of alpha=5e-4

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(5e-4)/2),]
  n_sig <- nrow(snp_sig)
  if (n_sig == 0){
    prop_bias <- -1
    mse <- -1
    prop_x <- -1
  }else{
    prop_bias <- sum(abs(snp_sig$beta) > abs(ss$true_beta[snp_sig$rsid]))/n_sig
    prop_x <- sum(abs(snp_sig$beta) > (abs(ss$true_beta[snp_sig$rsid]) + 1.96*ss$se[snp_sig$rsid]))/n_sig
    mse <- mean((ss$true_beta[snp_sig$rsid]-snp_sig$beta)^2)
  }
  return(c(n_sig,prop_bias,prop_x,mse))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

n_sig <- c(rep(0,nrow(sim_params)))
prop_bias <- c(rep(0,nrow(sim_params)))
prop_x <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  n_sig[i] <- res[[i]][1]
  prop_bias[i] <- res[[i]][2]
  prop_x[i] <- res[[i]][3]
  mse[i] <- res[[i]][4]
}

results <- cbind(sim_params,n_sig,prop_bias,prop_x,mse)
ave_res <- ave_results1(results,tot_sim)
write.csv(ave_res, "results/norm_nsig_prop_bias_5e-4.csv")


##############################################################################
## SIMULATION SET-UP 3:
## Bimodal effect size distribution: when S=0, 50% of effect sizes are generated
## from a N(0,1) distribution while the other 50% are generated from a N(2.5,1)
## distribution - see simulate_ss_bim() in 'useful_funs.R' for more details
## Significance threshold of alpha=5e-8

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_bim(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  out <- simulate_est(ss)
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(5e-8)/2),]
  n_sig <- nrow(snp_sig)
  if (n_sig == 0){
    prop_bias <- -1
    mse <- -1
    prop_x <- -1
  }else{
    prop_bias <- sum(abs(snp_sig$beta) > abs(ss$true_beta[snp_sig$rsid]))/n_sig
    prop_x <- sum(abs(snp_sig$beta) > (abs(ss$true_beta[snp_sig$rsid]) + 1.96*ss$se[snp_sig$rsid]))/n_sig
    mse <- mean((ss$true_beta[snp_sig$rsid]-snp_sig$beta)^2)
  }
  return(c(n_sig,prop_bias,prop_x,mse))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

n_sig <- c(rep(0,nrow(sim_params)))
prop_bias <- c(rep(0,nrow(sim_params)))
prop_x <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  n_sig[i] <- res[[i]][1]
  prop_bias[i] <- res[[i]][2]
  prop_x[i] <- res[[i]][3]
  mse[i] <- res[[i]][4]
}

results <- cbind(sim_params,n_sig,prop_bias,prop_x,mse)
ave_res <- ave_results1(results,tot_sim)
write.csv(ave_res, "results/bim_nsig_prop_bias_5e-8.csv")



