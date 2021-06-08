## More formal analyses of factors which influence degree of winner's curse and
## optimal performance of EB...

## What we want to do is compute mse_naive/mse_rep for all DISCOVERY only
## scenarios... normal quantitative/binary trait, skewed and bimodal
## quantitative trait .. do this first! compare lets say mse for top 10 SNPs and
## top 100 SNPs - ensure fair comparisons across situations!

## then after.. when does EB behave optimally?? repeat but with
## mse_EB/mse_true_bayes .. compute this again for top 10 SNPs and top 100 SNPs



## normal distribution, 1e+6 SNPs, 5e-8 threshold


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


true_bayes <- function(ss, summary_disc, alpha=5e-8){
  stats <- data.frame(rsid=summary_disc$rsid, true_beta = ss$true_beta, beta = summary_disc$beta, se = ss$se, mu = ss$true_beta/ss$se, z = summary_disc$beta/ss$se)
  #if (sum(abs(stats$z) > qnorm(1-(alpha)/2))== 0){return(stats)}

  mu_mid <- seq(from=(min(stats$mu)-0.01),to=(max(stats$mu)+0.01),by=0.01)
  mu_mid <- round(mu_mid,2)
  prob_mid <- c(rep(0,length(mu_mid)))

  for (i in 1:length(prob_mid)){
    x <- (min(stats$mu)-0.01) + (i-1)*0.01
    prob_mid[i] <- (sum(stats$mu > x-0.005 & stats$mu < x+0.005))/nrow(stats)
  }

  mu_prob <- data.frame(mu_mid,prob_mid)

  stats <- dplyr::arrange(stats,desc(abs(stats$z)))
  stats_sig <- stats[1:100,]
  stats_sig$mu_bayes <- c(rep(0,nrow(stats_sig)))

  for(i in 1:nrow(stats_sig)){
    if (sum(mu_prob$prob_mid*dnorm(stats_sig$z[i] - mu_prob$mu_mid)) == 0){stats_sig$mu_bayes[i] <- stats_sig$z[i]}else{
      stats_sig$mu_bayes[i] <- sum(mu_prob$mu_mid*mu_prob$prob_mid*dnorm(stats_sig$z[i] - mu_prob$mu_mid))/sum(mu_prob$prob_mid*dnorm(stats_sig$z[i] - mu_prob$mu_mid))
    }
  }

  stats_sig$beta_bayes <- stats_sig$mu_bayes*stats_sig$se
  return(stats_sig)

}


mse_sig_evaluate_10 <- function(out,true_beta,i=4,alpha=5e-8){
  out <- dplyr::arrange(out,desc(abs(out$beta/out$se)))
  snp_sig <- out[1:10,]
  mse_sig <- mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2)
  return(mse_sig)
}

mse_sig_evaluate_100 <- function(out,true_beta,i=4,alpha=5e-8){
  out <- dplyr::arrange(out,desc(abs(out$beta/out$se)))
  snp_sig <- out[1:100,]
  mse_sig <- mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2)
  return(mse_sig)
}


## consider different architectures..


run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ## Quantitative trait - normal distribution
  ss_norm <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats_norm <- simulate_est(ss_norm)

  naive_mse_10_norm <- mse_sig_evaluate_10(disc_stats_norm,ss_norm$true_beta,i=2)
  naive_mse_100_norm <- mse_sig_evaluate_100(disc_stats_norm,ss_norm$true_beta,i=2)

  ss2_norm <- data.frame(true_beta=ss_norm$true_beta,se=ss_norm$rep_se)
  rep_stats_norm <- simulate_est(ss2_norm)
  out_rep_norm <- cbind(disc_stats_norm,rep_stats_norm$beta)
  mse_rep_10_norm <- mse_sig_evaluate_10(out_rep_norm,ss_norm$true_beta,alpha=5e-8)
  mse_rep_100_norm <- mse_sig_evaluate_100(out_rep_norm,ss_norm$true_beta,alpha=5e-8)

  wc_10_norm <- 100*((naive_mse_10_norm/mse_rep_10_norm)-1)
  wc_100_norm <- 100*((naive_mse_100_norm/mse_rep_100_norm)-1)

  ## Quantitative trait - bimodal distribution
  ss_bim <- simulate_ss_bim(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats_bim <- simulate_est(ss_bim)

  naive_mse_10_bim <- mse_sig_evaluate_10(disc_stats_bim,ss_bim$true_beta,i=2)
  naive_mse_100_bim <- mse_sig_evaluate_100(disc_stats_bim,ss_bim$true_beta,i=2)

  ss2_bim <- data.frame(true_beta=ss_bim$true_beta,se=ss_bim$se)
  rep_stats_bim <- simulate_est(ss2_bim)
  out_rep_bim <- cbind(disc_stats_bim,rep_stats_bim$beta)
  mse_rep_10_bim <- mse_sig_evaluate_10(out_rep_bim,ss_bim$true_beta,alpha=5e-8)
  mse_rep_100_bim <- mse_sig_evaluate_100(out_rep_bim,ss_bim$true_beta,alpha=5e-8)

  wc_10_bim <- (naive_mse_10_bim/mse_rep_10_bim)
  wc_100_bim <- (naive_mse_100_bim/mse_rep_100_bim)

  ## Quantitative trait - skewed distribution
  ss_exp <- simulate_ss_exp(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats_exp <- simulate_est(ss_exp)

  naive_mse_10_exp <- mse_sig_evaluate_10(disc_stats_exp,ss_exp$true_beta,i=2)
  naive_mse_100_exp <- mse_sig_evaluate_100(disc_stats_exp,ss_exp$true_beta,i=2)

  ss2_exp <- data.frame(true_beta=ss_exp$true_beta,se=ss_exp$se)
  rep_stats_exp <- simulate_est(ss2_exp)
  out_rep_exp <- cbind(disc_stats_exp,rep_stats_exp$beta)
  mse_rep_10_exp <- mse_sig_evaluate_10(out_rep_exp,ss_exp$true_beta,alpha=5e-8)
  mse_rep_100_exp <- mse_sig_evaluate_100(out_rep_exp,ss_exp$true_beta,alpha=5e-8)

  wc_10_exp <- (naive_mse_10_exp/mse_rep_10_exp)
  wc_100_exp <- (naive_mse_100_exp/mse_rep_100_exp)


  return(c(wc_10_norm,wc_100_norm,wc_10_bim,wc_100_bim,wc_10_exp,wc_100_exp))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)


################################################################################

## Organising results:


ave_results_2 <- function(res_vec, n_sim){
  ave_wc_10_norm <- c(rep(0,nrow(res_vec)/n_sim))
  sd_wc_10_norm <- c(rep(0,nrow(res_vec)/n_sim))
  ave_wc_100_norm <- c(rep(0,nrow(res_vec)/n_sim))
  sd_wc_100_norm <- c(rep(0,nrow(res_vec)/n_sim))
  ave_wc_10_bim <- c(rep(0,nrow(res_vec)/n_sim))
  sd_wc_10_bim <- c(rep(0,nrow(res_vec)/n_sim))
  ave_wc_100_bim <- c(rep(0,nrow(res_vec)/n_sim))
  sd_wc_100_bim <- c(rep(0,nrow(res_vec)/n_sim))
  ave_wc_10_exp <- c(rep(0,nrow(res_vec)/n_sim))
  sd_wc_10_exp <- c(rep(0,nrow(res_vec)/n_sim))
  ave_wc_100_exp <- c(rep(0,nrow(res_vec)/n_sim))
  sd_wc_100_exp <- c(rep(0,nrow(res_vec)/n_sim))
  for (i in 1:(nrow(res_vec)/n_sim)){
    ave_wc_10_norm[i] <- average(res_vec$wc_10_norm[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_wc_10_norm[i] <- std(res_vec$wc_10_norm[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_wc_100_norm[i] <- average(res_vec$wc_100_norm[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_wc_100_norm[i] <- std(res_vec$wc_100_norm[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_wc_10_bim[i] <- average(res_vec$wc_10_bim[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_wc_10_bim[i] <- std(res_vec$wc_10_bim[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_wc_100_bim[i] <- average(res_vec$wc_100_bim[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_wc_100_bim[i] <- std(res_vec$wc_100_bim[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_wc_10_exp[i] <- average(res_vec$wc_10_exp[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_wc_10_exp[i] <- std(res_vec$wc_10_exp[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_wc_100_exp[i] <- average(res_vec$wc_100_exp[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_wc_100_exp[i] <- std(res_vec$wc_100_exp[(i*n_sim-(n_sim-1)):(i*n_sim)])
  }
  keep <- seq(from=1,to=nrow(res_vec)-n_sim+1,by=n_sim)
  res_vec_ave <- res_vec[rownames(res_vec) %in% keep,]
  res_vec_ave$wc_10_norm <- ave_wc_10_norm
  res_vec_ave$wc_10_norm_error <- sd_wc_10_norm
  res_vec_ave$wc_100_norm <- ave_wc_100_norm
  res_vec_ave$wc_100_norm_error <- sd_wc_100_norm
  res_vec_ave$wc_10_bim <- ave_wc_10_bim
  res_vec_ave$wc_10_bim_error <- sd_wc_10_bim
  res_vec_ave$wc_100_bim <- ave_wc_100_bim
  res_vec_ave$wc_100_bim_error <- sd_wc_100_bim
  res_vec_ave$wc_10_exp <- ave_wc_10_exp
  res_vec_ave$wc_10_exp_error <- sd_wc_10_exp
  res_vec_ave$wc_100_exp <- ave_wc_100_exp
  res_vec_ave$wc_100_exp_error <- sd_wc_100_exp
  res_vec_ave <- subset(res_vec_ave, select = -sim )
  return(res_vec_ave)
}



wc_10_norm <- c(rep(0,nrow(sim_params)))
wc_100_norm <- c(rep(0,nrow(sim_params)))
wc_10_bim <- c(rep(0,nrow(sim_params)))
wc_100_bim <- c(rep(0,nrow(sim_params)))
wc_10_exp <- c(rep(0,nrow(sim_params)))
wc_100_exp <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  wc_10_norm[i] <- res[[i]][1]
  wc_100_norm[i] <- res[[i]][2]
  wc_10_bim[i] <- res[[i]][3]
  wc_100_bim[i] <- res[[i]][4]
  wc_10_exp[i] <- res[[i]][5]
  wc_100_exp[i] <- res[[i]][6]
}
res_2 <- cbind(sim_params,wc_10_norm,wc_100_norm,wc_10_bim,wc_100_bim,wc_10_exp,wc_100_exp)
ave_res <- ave_results_2(res_2,tot_sim)

write.csv(ave_res,"results/wc_measure_rank_100sim.csv")


################################################################################

## optimality of empirical Bayes

set.seed(1998)

## Total number of simulations: 100
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

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ## normal quantitative
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  mse_EB_10 <- mse_sig_evaluate_10(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB_100 <- mse_sig_evaluate_100(out_EB,ss$true_beta,alpha=5e-8)

  ## True Bayes:
  out_bayes <- true_bayes(ss,disc_stats)
  mse_bayes_10 <- mse_sig_evaluate_10(out_bayes,ss$true_beta,i=8,alpha=5e-8)
  mse_bayes_100 <- mse_sig_evaluate_100(out_bayes,ss$true_beta,i=8,alpha=5e-8)

  perc_opt_EB_10_norm <- mse_EB_10/mse_bayes_10
  perc_opt_EB_100_norm <- mse_EB_100/mse_bayes_100

  ## bimodal quantitative
  ss <- simulate_ss_bim(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  mse_EB_10 <- mse_sig_evaluate_10(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB_100 <- mse_sig_evaluate_100(out_EB,ss$true_beta,alpha=5e-8)

  ## True Bayes:
  out_bayes <- true_bayes(ss,disc_stats)
  mse_bayes_10 <- mse_sig_evaluate_10(out_bayes,ss$true_beta,i=8,alpha=5e-8)
  mse_bayes_100 <- mse_sig_evaluate_100(out_bayes,ss$true_beta,i=8,alpha=5e-8)

  perc_opt_EB_10_bim <- mse_EB_10/mse_bayes_10
  perc_opt_EB_100_bim <- mse_EB_100/mse_bayes_100

  ## skewed quantitative
  ss <- simulate_ss_exp(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  mse_EB_10 <- mse_sig_evaluate_10(out_EB,ss$true_beta,alpha=5e-8)
  mse_EB_100 <- mse_sig_evaluate_100(out_EB,ss$true_beta,alpha=5e-8)

  ## True Bayes:
  out_bayes <- true_bayes(ss,disc_stats)
  mse_bayes_10 <- mse_sig_evaluate_10(out_bayes,ss$true_beta,i=8,alpha=5e-8)
  mse_bayes_100 <- mse_sig_evaluate_100(out_bayes,ss$true_beta,i=8,alpha=5e-8)

  perc_opt_EB_10_skew <- mse_EB_10/mse_bayes_10
  perc_opt_EB_100_skew <- mse_EB_100/mse_bayes_100


  return(c(perc_opt_EB_10_norm,perc_opt_EB_100_norm,perc_opt_EB_10_bim,perc_opt_EB_100_bim,perc_opt_EB_10_skew,perc_opt_EB_100_skew))
}
res <- mclapply(1:nrow(sim_params), function(i){do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

ave_results_EB <- function(res_vec, n_sim){
  ave_perc_opt_EB_10_norm <- c(rep(0,nrow(res_vec)/n_sim))
  sd_perc_opt_EB_10_norm <- c(rep(0,nrow(res_vec)/n_sim))
  ave_perc_opt_EB_100_norm <- c(rep(0,nrow(res_vec)/n_sim))
  sd_perc_opt_EB_100_norm <- c(rep(0,nrow(res_vec)/n_sim))
  ave_perc_opt_EB_10_bim <- c(rep(0,nrow(res_vec)/n_sim))
  sd_perc_opt_EB_10_bim <- c(rep(0,nrow(res_vec)/n_sim))
  ave_perc_opt_EB_100_bim <- c(rep(0,nrow(res_vec)/n_sim))
  sd_perc_opt_EB_100_bim <- c(rep(0,nrow(res_vec)/n_sim))
  ave_perc_opt_EB_10_skew <- c(rep(0,nrow(res_vec)/n_sim))
  sd_perc_opt_EB_10_skew <- c(rep(0,nrow(res_vec)/n_sim))
  ave_perc_opt_EB_100_skew <- c(rep(0,nrow(res_vec)/n_sim))
  sd_perc_opt_EB_100_skew <- c(rep(0,nrow(res_vec)/n_sim))
  for (i in 1:(nrow(res_vec)/n_sim)){
    ave_perc_opt_EB_10_norm[i] <- average(res_vec$perc_opt_EB_10_norm[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_perc_opt_EB_10_norm[i] <- std(res_vec$perc_opt_EB_10_norm[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_perc_opt_EB_100_norm[i] <- average(res_vec$perc_opt_EB_100_norm[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_perc_opt_EB_100_norm[i] <- std(res_vec$perc_opt_EB_100_norm[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_perc_opt_EB_10_bim[i] <- average(res_vec$perc_opt_EB_10_bim[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_perc_opt_EB_10_bim[i] <- std(res_vec$perc_opt_EB_10_bim[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_perc_opt_EB_100_bim[i] <- average(res_vec$perc_opt_EB_100_bim[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_perc_opt_EB_100_bim[i] <- std(res_vec$perc_opt_EB_100_bim[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_perc_opt_EB_10_skew[i] <- average(res_vec$perc_opt_EB_10_skew[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_perc_opt_EB_10_skew[i] <- std(res_vec$perc_opt_EB_10_skew[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_perc_opt_EB_100_skew[i] <- average(res_vec$perc_opt_EB_100_skew[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_perc_opt_EB_100_skew[i] <- std(res_vec$perc_opt_EB_100_skew[(i*n_sim-(n_sim-1)):(i*n_sim)])
  }
  keep <- seq(from=1,to=nrow(res_vec)-n_sim+1,by=n_sim)
  res_vec_ave <- res_vec[rownames(res_vec) %in% keep,]
  res_vec_ave$perc_opt_EB_10_norm <- ave_perc_opt_EB_10_norm
  res_vec_ave$perc_opt_EB_10_norm_error <- sd_perc_opt_EB_10_norm
  res_vec_ave$perc_opt_EB_100_norm <- ave_perc_opt_EB_100_norm
  res_vec_ave$perc_opt_EB_100_norm_error <- sd_perc_opt_EB_100_norm
  res_vec_ave$perc_opt_EB_10_bim <- ave_perc_opt_EB_10_bim
  res_vec_ave$perc_opt_EB_10_bim_error <- sd_perc_opt_EB_10_bim
  res_vec_ave$perc_opt_EB_100_bim <- ave_perc_opt_EB_100_bim
  res_vec_ave$perc_opt_EB_100_bim_error <- sd_perc_opt_EB_100_bim
  res_vec_ave$perc_opt_EB_10_skew <- ave_perc_opt_EB_10_skew
  res_vec_ave$perc_opt_EB_10_skew_error <- sd_perc_opt_EB_10_skew
  res_vec_ave$perc_opt_EB_100_skew <- ave_perc_opt_EB_100_skew
  res_vec_ave$perc_opt_EB_100_skew_error <- sd_perc_opt_EB_100_skew
  res_vec_ave <- subset(res_vec_ave, select = -sim )
  return(res_vec_ave)
}


perc_opt_EB_10_norm <- c(rep(0,nrow(sim_params)))
perc_opt_EB_100_norm <- c(rep(0,nrow(sim_params)))
perc_opt_EB_10_bim <- c(rep(0,nrow(sim_params)))
perc_opt_EB_100_bim <- c(rep(0,nrow(sim_params)))
perc_opt_EB_10_skew <- c(rep(0,nrow(sim_params)))
perc_opt_EB_100_skew <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  perc_opt_EB_10_norm[i] <- res[[i]][1]
  perc_opt_EB_100_norm[i] <- res[[i]][2]
  perc_opt_EB_10_bim[i] <- res[[i]][3]
  perc_opt_EB_100_bim[i] <- res[[i]][4]
  perc_opt_EB_10_skew[i] <- res[[i]][5]
  perc_opt_EB_100_skew[i] <- res[[i]][6]
}
res_2 <- cbind(sim_params,perc_opt_EB_10_norm,perc_opt_EB_100_norm,perc_opt_EB_10_bim,perc_opt_EB_100_bim,perc_opt_EB_10_skew,perc_opt_EB_100_skew)
ave_res <- ave_results_EB(res_2,tot_sim)

write.csv(ave_res,"results/EB_perform_measure_rank_10sim.csv")


################################################################################


