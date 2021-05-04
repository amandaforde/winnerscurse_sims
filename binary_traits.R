## SIMULATION SET-UP 5: BINARY TRAIT


##library(devtools)
##devtools::install_github("amandaforde/winnerscurse")
library(winnerscurse)
library(tidyverse)
library(parallel)


## Extra functions required for simulating binary trait:

###  convert logit values to probabilities
expit <- function(x){exp(x)/(1+exp(x))}

## find Gamma0 for logistic model so prevalence is 0.1
myfunc <- function(MAF,Gamma0, Gamma1, prev=0.1){
  MAF^2*expit(Gamma0+Gamma1*2)+2*MAF*(1-MAF)*expit(Gamma0+Gamma1)+(1-MAF)^2*expit(Gamma0)-prev
}

## calculate asymptotic version of (X^t W X)
asymp_var_logistic <- function(n,G_prob,Gamma_0,Gamma_1){
  N <- length(Gamma_0)
  disease_probs <- expit(outer(Gamma_0,c(1,1,1)) +  outer(Gamma_1,c(0,1,2)))
  diag_weights <- disease_probs*(1-disease_probs)
  a <- n*apply(G_prob*diag_weights,1,sum)
  b <- n*apply(G_prob*diag_weights*matrix(rep(c(0,1,2),N),byrow=TRUE,nrow=N),1,sum)
  c <-  b
  d <- n*apply(G_prob*diag_weights*matrix(rep(c(0,1,4),N),byrow=TRUE,nrow=N),1,sum)
  ##  invert matrix and take element for Gamma1
  return(a/(a*d-b*c))

}


## Grid-like computation of beta_0 corresponding to values for maf and true_beta
## with fixed prevalence of 0.1 saves computational time
params <- expand.grid(
  maf = seq(0.01,0.5,by=0.01),
  true_beta = seq(-1.5,1.5,by=0.01)
)
params$maf <- round(params$maf,2)
params$true_beta <- round(params$true_beta,2)
params$Gamma_0 <- c(rep(0,nrow(params)))
for(i in 1:nrow(params)) {params$Gamma_0[i] <- uniroot(myfunc,Gamma1=params$true_beta[i],MAF=params$maf[i],prev=0.1,lower=-10, upper=10)$root}


simulate_ss_bin <- function(H,Pi,nid,sc){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- rnorm(effect_snps,0,sd=sqrt((2*maf*(1-maf))^sc))
  scaling <- (1.6^2*H^2)/((sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2))*(1-H^2))
  true_beta <- true_beta*sqrt(scaling)
  true_beta <- c(true_beta, rep(0,n_snps-effect_snps))

  rounded_maf <- round(maf,2)
  rounded_true_beta <- round(true_beta,2)
  Gamma0 <- c(rep(0,n_snps))

  dat <- data.frame(rounded_maf, rounded_true_beta)
  for(i in 1:n_snps) Gamma0[i] <- params$Gamma_0[7500 + 5000*rounded_true_beta[i] + 100*rounded_maf[i]]

  ##  assuming HWE
  G_prob <- cbind(maf^2,2*maf*(1-maf),(1-maf)^2)
  Gamma1 <- true_beta
  var_Gamma_y <- asymp_var_logistic(nid,G_prob,Gamma0,Gamma1)
  se <- sqrt(var_Gamma_y)
  stats <- data.frame(true_beta,se)
  return(stats)
}


################################################################################

## SIMULATION SET-UP 5a.:
## Obtain the average number of significant SNPs, the average proportion of
## these SNPs for which their association estimate is more extreme than their
## true effect size and the average proportion of these SNPs which are
## significantly overexaggerated over simulations in which at least one
## significant SNP has been detected
## Normal effect size distribution
## Significance threshold of alpha=5e-8


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
  prop_effect = c(0.01, 0.001),
  S = c(-1, 0, 1)
)

## Run 'useful_funs.R' here in order to define functions required for the
## simulations below.
## NOTE: Simulations currently being run on Windows, hence mc.cores=1.


run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_bin(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
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
write.csv(ave_res, "results/binary_norm_nsig_prop_bias_5e-8_10sim.csv")


##############################################################################

## SIMULATION SET-UP 5b.:
## Obtain the average number of significant SNPs, the average proportion of
## these SNPs for which their association estimate is more extreme than their
## true effect size and the average proportion of these SNPs which are
## significantly overexaggerated over simulations in which at least one
## significant SNP has been detected
## Normal effect size distribution
## Significance threshold of alpha=5e-4

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_bin(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
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
write.csv(ave_res, "results/binary_norm_nsig_prop_bias_5e-4_10sim.csv")

################################################################################


## SIMULATION SET-UP 5c.:
## Evaluating Winner's Curse adjustment methods with a binary trait
## Discovery GWAS only
## Normal effect size distribution
## Significance threshold of alpha=5e-8


set.seed(1998)

## Bias Evaluation Metrics:
## 1. Evaluating the fraction of significant SNPs that have been improved due to
## method implementation - effect size estimates have been adjusted so that they
## are closer to the true effect size
## 2. Evaluating the change in average MSE of significant SNPs due to method
## implementation
## 3. Evaluating the relative change in average MSE of significant SNPs due to
## method implementation

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_bin_v2(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  rep_stats <- simulate_est(ss)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,alpha=5e-8)

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
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,alpha=5e-8)

  return(c(flb_EB,mse_EB,rel_mse_EB,flb_FIQT,mse_FIQT,rel_mse_FIQT,flb_BR,mse_BR,rel_mse_BR,flb_cl1,mse_cl1,rel_mse_cl1,flb_cl2,mse_cl2,rel_mse_cl2,flb_cl3,mse_cl3,rel_mse_cl3,flb_rep,mse_rep,rel_mse_rep))
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

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][4]
  mse[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_FIQT <- cbind(sim_params,flb,mse,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## Bootstrap:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rel_mse[i] <- res[[i]][9]
}
res_BR <- cbind(sim_params,flb,mse,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)

## Conditional Likelihood 1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][10]
  mse[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_cl1 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## Conditional Likelihood 2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rel_mse[i] <- res[[i]][15]
}
res_cl2 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## Conditional Likelihood 3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][16]
  mse[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_cl3 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rel_mse[i] <- res[[i]][21]
}
res_rep <- cbind(sim_params,flb,mse,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)


## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3,ave_res_rep)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)),rep("rep",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/binary_norm_5e-8_10sim.csv")

################################################################################

