## COMPARING TRUE BAYES RULE WITH EMPIRICAL BAYES

## SIMULATION SET-UP 1:
## Discovery GWAS only
## Quantitative trait
## Normal effect size distribution
## Significance threshold of alpha=5e-8

################################################################################

#library(devtools)
#devtools::install_github("amandaforde/winnerscurse")
library(winnerscurse)
library(tidyverse)
library(parallel)

set.seed(1998)

## Total number of simulations: 5
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
## 1. Evaluating the fraction of significant SNPs that have been improved due to
## method implementation - effect size estimates have been adjusted so that they
## are closer to the true effect size
## 2. Evaluating the change in average MSE of significant SNPs due to method
## implementation
## 3. Evaluating the relative change in average MSE of significant SNPs due to
## method implementation

## first write a function to apply the true_bayes rule!

true_bayes <- function(ss, summary_disc, alpha=5e-8){
  stats <- data.frame(rsid=summary_disc$rsid, true_beta = ss$true_beta, beta = summary_disc$beta, se = ss$se, mu = ss$true_beta/ss$se, z = summary_disc$beta/ss$se)
  if (sum(abs(stats$z) > qnorm(1-(alpha)/2))== 0){return(stats)}

  mu_mid <- seq(from=(min(stats$mu)-0.01),to=(max(stats$mu)+0.01),by=0.01)
  mu_mid <- round(mu_mid,2)
  prob_mid <- c(rep(0,length(mu_mid)))

  for (i in 1:length(prob_mid)){
    x <- (min(stats$mu)-0.01) + (i-1)*0.01
    prob_mid[i] <- (sum(stats$mu > x-0.005 & stats$mu < x+0.005))/nrow(stats)
  }

  mu_prob <- data.frame(mu_mid,prob_mid)

  stats_sig <- stats[abs(stats$z) > qnorm(1-(alpha)/2),]
  stats_sig$mu_bayes <- c(rep(0,nrow(stats_sig)))

  for(i in 1:nrow(stats_sig)){
    if (sum(mu_prob$prob_mid*dnorm(stats_sig$z[i] - mu_prob$mu_mid)) == 0){stats_sig$mu_bayes[i] <- stats_sig$z[i]}else{
      stats_sig$mu_bayes[i] <- sum(mu_prob$mu_mid*mu_prob$prob_mid*dnorm(stats_sig$z[i] - mu_prob$mu_mid))/sum(mu_prob$prob_mid*dnorm(stats_sig$z[i] - mu_prob$mu_mid))
    }
  }

  stats_sig$beta_bayes <- stats_sig$mu_bayes*stats_sig$se
  return(stats_sig)

}





run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,alpha=5e-8)

  ## True Bayes:
  out_bayes <- true_bayes(ss,disc_stats)
  flb_bayes <- frac_sig_less_bias(out_bayes,ss$true_beta,i=8,alpha=5e-8)
  mse_bayes <- mse_sig_improve(out_bayes,ss$true_beta,i=8,alpha=5e-8)
  rel_mse_bayes <- mse_sig_improve_per(out_bayes,ss$true_beta,i=8,alpha=5e-8)

  ## replication
  ss2 <- data.frame(true_beta=ss$true_beta,se=ss$rep_se)
  rep_stats <- simulate_est(ss2)
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,alpha=5e-8)

  return(c(flb_EB,mse_EB,rel_mse_EB,flb_bayes,mse_bayes,rel_mse_bayes,flb_rep,mse_rep,rel_mse_rep))
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

## TRUE BAYES:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][4]
  mse[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_bayes <- cbind(sim_params,flb,mse,rel_mse)
ave_res_bayes <- ave_results(res_bayes,tot_sim)


## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rel_mse[i] <- res[[i]][9]
}
res_rep <- cbind(sim_params,flb,mse,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)


## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_bayes,ave_res_rep)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("bayes",(nrow(sim_params)/tot_sim)),rep("rep",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/norm_5e-8_10sim_bayes_compare.csv")

################################################################################



## compare with old empirical bayes approach!

empirical_bayes_old <- function(summary_data)

{

  if(nrow(summary_data) == 1){return(summary_data)}

  z <- summary_data$beta/summary_data$se

  bins <- seq(min(z),max(z),length.out=120)
  mids <- (bins[-length(bins)]+bins[-1])/2
  counts <- graphics::hist(z,breaks=bins,plot=F)$counts

  most_extreme <- 10
  boundary_lower <- sort(z)[most_extreme]
  boundary_upper <- sort(z,decreasing=TRUE)[most_extreme]

  df <- 7
  AIC_vector <- c(rep(0,28))
  for (best_df in 3:30){
    model <- stats::glm(counts ~ splines::ns(mids,knots = (seq(from=boundary_lower,to=boundary_upper,length=best_df+1)[2:best_df]), Boundary.knots=c(boundary_lower,boundary_upper)),stats::poisson,weights=rep(10^-50,length(counts)))
    minus2loglike <- 10^50*(model$deviance)
    AIC_vector[best_df-2] <- minus2loglike + 2*(best_df-2)
  }
  df <- 2 + which.min(AIC_vector)

  f <- stats::glm(counts ~ splines::ns(mids,knots = (seq(from=boundary_lower,to=boundary_upper,length=df+1)[2:df]), Boundary.knots=c(boundary_lower,boundary_upper)),stats::poisson,weight=rep(10^-50,length(counts)))$fit
  f[f==0] <- min(f[f>0])

  log_f <- as.vector(log(f))
  diff <- diff(log_f)/diff(mids)
  mids2 <- (mids[-length(mids)]+mids[-1])/2

  diff_interpol <- stats::approx(mids2,diff,mids,rule=2,ties=mean)$y

  mids_est <- c(rep(0,length(mids)))
  mids_est[mids>0] <- pmax(0, mids[mids>0] + diff_interpol[mids>0])
  mids_est[mids<0] <- pmin(0, mids[mids<0] + diff_interpol[mids<0])

  z_hat <- stats::approx(mids,mids_est,z,rule=2,ties=mean)$y
  z_hat <- sign(z)*pmin(abs(z),abs(z_hat))

  beta_EB <- z_hat*summary_data$se
  summary_data <- cbind(summary_data,beta_EB)

  summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))

  return(summary_data)

}


run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,alpha=5e-8)

  ## True Bayes:
  out_EB_2 <- empirical_bayes_old(disc_stats)
  flb_EB_2 <- frac_sig_less_bias(out_EB_2,ss$true_beta,alpha=5e-8)
  mse_EB_2 <- mse_sig_improve(out_EB_2,ss$true_beta,alpha=5e-8)
  rel_mse_EB_2 <- mse_sig_improve_per(out_EB_2,ss$true_beta,alpha=5e-8)

  ## replication
  ss2 <- data.frame(true_beta=ss$true_beta,se=ss$rep_se)
  rep_stats <- simulate_est(ss2)
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,alpha=5e-8)

  return(c(flb_EB,mse_EB,rel_mse_EB,flb_EB_2,mse_EB_2,rel_mse_EB_2,flb_rep,mse_rep,rel_mse_rep))
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

## old empirical BAYES:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][4]
  mse[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_EB_2 <- cbind(sim_params,flb,mse,rel_mse)
ave_res_EB_2 <- ave_results(res_EB_2,tot_sim)


## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rel_mse[i] <- res[[i]][9]
}
res_rep <- cbind(sim_params,flb,mse,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)


## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_EB_2,ave_res_rep)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("old EB",(nrow(sim_params)/tot_sim)),rep("rep",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/norm_5e-8_10sim_EB_compare.csv")

################################################################################
