## WINNER'S CURSE SIMULATION STUDY SCRIPT 4:
## Method evaluation and comparison with independence and normal effect size
## distribution

## 1) Quantitative trait with normal effect size distribution - simulate_ss()

## Number of repetitions: 100
## Significance thresholds: alpha=5e-8 and alpha=5e-4
## Methods evaluated: Empirical Bayes, Empirical Bayes with fixed df-7, Empirical
## Bayes using scam, Empirical Bayes using gam and poisson, Empirical Bayes using
## gam and negative binomial, FIQT, Bootstrap, Conditional Likelihood methods
## Evaluation metrics: mse, rmse, rel_mse, bias, flb
## Assumption: SNPs are independent

## Evaluating relative improvement of average RMSE for significant SNPs due
## to method implementation
rmse_sig_improve_per <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm((alpha)/2, lower.tail=FALSE),]
  if (nrow(snp_sig) == 0){return(100)}
  mse_sig_improve <- (sqrt(mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2))-sqrt(mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2)))/(sqrt(mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2)))
  return(mse_sig_improve)
}

################################################################################

## Total number of simulations:
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

################################################################################
################################################################################

## 1A) Quantitative trait with normal effect size distribution - 5e-8

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,i=4,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,i=4,alpha=5e-8)
  rmse_EB <- mse_sig_improve_root(out_EB,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_up <- bias_sig_up(out_EB,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_down <- bias_sig_down(out_EB,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_EB <- rmse_sig_improve_per(out_EB,ss$true_beta,i=4,alpha=5e-8)

  ## Empirical Bayes df=7:
  out_EB_df <- empirical_bayes(disc_stats, method="fix_df")
  flb_EB_df <- frac_sig_less_bias(out_EB_df,ss$true_beta,i=4,alpha=5e-8)
  mse_EB_df <- mse_sig_improve(out_EB_df,ss$true_beta,i=4,alpha=5e-8)
  rmse_EB_df <- mse_sig_improve_root(out_EB_df,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_df_up <- bias_sig_up(out_EB_df,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_df_down <- bias_sig_down(out_EB_df,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_EB_df <- rmse_sig_improve_per(out_EB_df,ss$true_beta,i=4,alpha=5e-8)

  ## Empirical Bayes scam:
  out_EB_scam <- empirical_bayes(disc_stats, method="scam")
  flb_EB_scam <- frac_sig_less_bias(out_EB_scam,ss$true_beta,i=4,alpha=5e-8)
  mse_EB_scam <- mse_sig_improve(out_EB_scam,ss$true_beta,i=4,alpha=5e-8)
  rmse_EB_scam <- mse_sig_improve_root(out_EB_scam,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_scam_up <- bias_sig_up(out_EB_scam,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_scam_down <- bias_sig_down(out_EB_scam,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_EB_scam <- rmse_sig_improve_per(out_EB_scam,ss$true_beta,i=4,alpha=5e-8)

  ## Empirical Bayes gam_po:
  out_EB_gam_po <- empirical_bayes(disc_stats, method="gam_po")
  flb_EB_gam_po <- frac_sig_less_bias(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-8)
  mse_EB_gam_po <- mse_sig_improve(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-8)
  rmse_EB_gam_po <- mse_sig_improve_root(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_gam_po_up <- bias_sig_up(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_gam_po_down <- bias_sig_down(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_EB_gam_po <- rmse_sig_improve_per(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-8)

  ## Empirical Bayes gam_nb:
  out_EB_gam_nb <- empirical_bayes(disc_stats, method="gam_nb")
  flb_EB_gam_nb <- frac_sig_less_bias(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-8)
  mse_EB_gam_nb <- mse_sig_improve(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-8)
  rmse_EB_gam_nb <- mse_sig_improve_root(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_gam_nb_up <- bias_sig_up(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_gam_nb_down <- bias_sig_down(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_EB_gam_nb <- rmse_sig_improve_per(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-8)

  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  mse_FIQT <- mse_sig_improve(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  rmse_FIQT <- mse_sig_improve_root(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  bias_FIQT_up <- bias_sig_up(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  bias_FIQT_down <- bias_sig_down(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_FIQT <- rmse_sig_improve_per(out_FIQT,ss$true_beta,i=4,alpha=5e-8)

  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,ss$true_beta,i=4,alpha=5e-8)
  mse_BR <- mse_sig_improve(out_BR,ss$true_beta,i=4,alpha=5e-8)
  rmse_BR <- mse_sig_improve_root(out_BR,ss$true_beta,i=4,alpha=5e-8)
  bias_BR_up <- bias_sig_up(out_BR,ss$true_beta,i=4,alpha=5e-8)
  bias_BR_down <- bias_sig_down(out_BR,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_BR <- rmse_sig_improve_per(out_BR,ss$true_beta,i=4,alpha=5e-8)

  ## cl1:
  out_cl <- conditional_likelihood_ind(disc_stats,alpha=5e-8)
  flb_cl1 <- frac_sig_less_bias(out_cl,ss$true_beta,i=4,alpha=5e-8)
  mse_cl1 <- mse_sig_improve(out_cl,ss$true_beta,i=4,alpha=5e-8)
  rmse_cl1 <- mse_sig_improve_root(out_cl,ss$true_beta,i=4,alpha=5e-8)
  bias_cl1_up <- bias_sig_up(out_cl,ss$true_beta,i=4,alpha=5e-8)
  bias_cl1_down <- bias_sig_down(out_cl,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_cl1 <- rmse_sig_improve_per(out_cl,ss$true_beta,i=4,alpha=5e-8)

  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=5)
  mse_cl2 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=5)
  rmse_cl2 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=5)
  bias_cl2_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=5)
  bias_cl2_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=5)
  rel_mse_cl2 <- rmse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=5)

  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=6)
  mse_cl3 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rmse_cl3 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=6)
  bias_cl3_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=6)
  bias_cl3_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_cl3 <- rmse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=6)

  ## replication
  ss2 <- data.frame(true_beta=ss$true_beta,se=ss$rep_se)
  rep_stats <- simulate_est(ss2)
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,i=4,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,i=4,alpha=5e-8)
  rmse_rep <- mse_sig_improve_root(out_rep,ss$true_beta,i=4,alpha=5e-8)
  bias_rep_up <- bias_sig_up(out_rep,ss$true_beta,i=4,alpha=5e-8)
  bias_rep_down <- bias_sig_down(out_rep,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_rep <- rmse_sig_improve_per(out_rep,ss$true_beta,i=4,alpha=5e-8)

  ## naive - only interested in bias (improvement should be 0):
  flb_naive <- frac_sig_less_bias(out_cl,disc_stats$true_beta,alpha=5e-8,i=2)
  mse_naive <- mse_sig_improve(out_cl,disc_stats$true_beta,alpha=5e-8,i=2)
  rmse_naive <- mse_sig_improve_root(out_cl,disc_stats$true_beta,alpha=5e-8,i=2)
  bias_naive_up <- bias_sig_up(out_cl,disc_stats$true_beta,alpha=5e-8,i=2)
  bias_naive_down <- bias_sig_down(out_cl,disc_stats$true_beta,alpha=5e-8,i=2)
  rel_mse_naive <- rmse_sig_improve_per(out_cl,disc_stats$true_beta,alpha=5e-8,i=2)

  return(c(flb_EB,mse_EB,rmse_EB,bias_EB_up,bias_EB_down,rel_mse_EB,flb_EB_df,mse_EB_df,rmse_EB_df,bias_EB_df_up,bias_EB_df_down,rel_mse_EB_df,flb_EB_scam,mse_EB_scam,rmse_EB_scam,bias_EB_scam_up,bias_EB_scam_down,rel_mse_EB_scam,flb_EB_gam_po,mse_EB_gam_po,rmse_EB_gam_po,bias_EB_gam_po_up,bias_EB_gam_po_down,rel_mse_EB_gam_po,flb_EB_gam_nb,mse_EB_gam_nb,rmse_EB_gam_nb,bias_EB_gam_nb_up,bias_EB_gam_nb_down,rel_mse_EB_gam_nb,flb_FIQT,mse_FIQT,rmse_FIQT,bias_FIQT_up,bias_FIQT_down,rel_mse_FIQT,flb_BR,mse_BR,rmse_BR,bias_BR_up,bias_BR_down,rel_mse_BR,flb_cl1,mse_cl1,rmse_cl1,bias_cl1_up,bias_cl1_down,rel_mse_cl1,flb_cl2,mse_cl2,rmse_cl2,bias_cl2_up,bias_cl2_down,rel_mse_cl2,flb_cl3,mse_cl3,rmse_cl3,bias_cl3_up,bias_cl3_down,rel_mse_cl3,flb_rep,mse_rep,rmse_rep,bias_rep_up,bias_rep_down,rel_mse_rep,flb_naive,mse_naive,rmse_naive,bias_naive_up,bias_naive_down,rel_mse_naive))
}
res <- mclapply(1:nrow(sim_params), function(i){
  print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################

## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  bias_up[i] <- res[[i]][4]
  bias_down[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_EB <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)

## Empirical Bayes df=7:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rmse[i] <- res[[i]][9]
  bias_up[i] <- res[[i]][10]
  bias_down[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_EB_df <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB_df <- ave_results(res_EB_df,tot_sim)

## Empirical Bayes scam:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rmse[i] <- res[[i]][15]
  bias_up[i] <- res[[i]][16]
  bias_down[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_EB_scam <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB_scam <- ave_results(res_EB_scam,tot_sim)

## Empirical Bayes gam_po:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rmse[i] <- res[[i]][21]
  bias_up[i] <- res[[i]][22]
  bias_down[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_EB_gam_po <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB_gam_po <- ave_results(res_EB_gam_po,tot_sim)

## Empirical Bayes gam_nb:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rmse[i] <- res[[i]][27]
  bias_up[i] <- res[[i]][28]
  bias_down[i] <- res[[i]][29]
  rel_mse[i] <- res[[i]][30]
}
res_EB_gam_nb <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB_gam_nb <- ave_results(res_EB_gam_nb,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][31]
  mse[i] <- res[[i]][32]
  rmse[i] <- res[[i]][33]
  bias_up[i] <- res[[i]][34]
  bias_down[i] <- res[[i]][35]
  rel_mse[i] <- res[[i]][36]
}
res_FIQT <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## Bootstrap:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][37]
  mse[i] <- res[[i]][38]
  rmse[i] <- res[[i]][39]
  bias_up[i] <- res[[i]][40]
  bias_down[i] <- res[[i]][41]
  rel_mse[i] <- res[[i]][42]
}
res_BR <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)

## Conditional Likelihood 1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][43]
  mse[i] <- res[[i]][44]
  rmse[i] <- res[[i]][45]
  bias_up[i] <- res[[i]][46]
  bias_down[i] <- res[[i]][47]
  rel_mse[i] <- res[[i]][48]
}
res_cl1 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## Conditional Likelihood 2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][49]
  mse[i] <- res[[i]][50]
  rmse[i] <- res[[i]][51]
  bias_up[i] <- res[[i]][52]
  bias_down[i] <- res[[i]][53]
  rel_mse[i] <- res[[i]][54]
}
res_cl2 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## Conditional Likelihood 3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][55]
  mse[i] <- res[[i]][56]
  rmse[i] <- res[[i]][57]
  bias_up[i] <- res[[i]][58]
  bias_down[i] <- res[[i]][59]
  rel_mse[i] <- res[[i]][60]
}
res_cl3 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][61]
  mse[i] <- res[[i]][62]
  rmse[i] <- res[[i]][63]
  bias_up[i] <- res[[i]][64]
  bias_down[i] <- res[[i]][65]
  rel_mse[i] <- res[[i]][66]
}
res_rep <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## Naive:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][67]
  mse[i] <- res[[i]][68]
  rmse[i] <- res[[i]][69]
  bias_up[i] <- res[[i]][70]
  bias_down[i] <- res[[i]][71]
  rel_mse[i] <- res[[i]][72]
}
res_naive <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_naive <- ave_results(res_naive,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_EB_df,ave_res_EB_scam,ave_res_EB_gam_po,ave_res_EB_gam_nb,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3, ave_res_rep,ave_res_naive)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("EB_df",(nrow(sim_params)/tot_sim)),rep("EB_scam",(nrow(sim_params)/tot_sim)),rep("EB_gam_po",(nrow(sim_params)/tot_sim)),rep("EB_gam_nb",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)), rep("rep",(nrow(sim_params)/tot_sim)),rep("naive",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/norm_5e-8_100sim_ave.csv")

## Combine all results:
results_all <- rbind(res_EB,res_EB_df,res_EB_scam,res_EB_gam_po,res_EB_gam_nb,res_FIQT,res_BR,res_cl1,res_cl2,res_cl3,res_rep,res_naive)
results_all$method <- c(rep("EB",nrow(sim_params)),rep("EB_df",nrow(sim_params)),rep("EB_scam",nrow(sim_params)),rep("EB_gam_po",nrow(sim_params)),rep("EB_gam_nb",nrow(sim_params)),rep("FIQT",nrow(sim_params)),rep("BR",nrow(sim_params)),rep("cl1",nrow(sim_params)),rep("cl2",nrow(sim_params)),rep("cl3",nrow(sim_params)),rep("rep",nrow(sim_params)),rep("naive",nrow(sim_params)))
write.csv(results_all,"results/norm_5e-8_100sim_all.csv")

print("PART 1A) complete!")

################################################################################
################################################################################

## 1B) Quantitative trait with normal effect size distribution - 5e-4

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,i=4,alpha=5e-4)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,i=4,alpha=5e-4)
  rmse_EB <- mse_sig_improve_root(out_EB,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_up <- bias_sig_up(out_EB,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_down <- bias_sig_down(out_EB,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,i=4,alpha=5e-4)

  ## Empirical Bayes df=7:
  out_EB_df <- empirical_bayes(disc_stats, method="fix_df")
  flb_EB_df <- frac_sig_less_bias(out_EB_df,ss$true_beta,i=4,alpha=5e-4)
  mse_EB_df <- mse_sig_improve(out_EB_df,ss$true_beta,i=4,alpha=5e-4)
  rmse_EB_df <- mse_sig_improve_root(out_EB_df,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_df_up <- bias_sig_up(out_EB_df,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_df_down <- bias_sig_down(out_EB_df,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_EB_df <- mse_sig_improve_per(out_EB_df,ss$true_beta,i=4,alpha=5e-4)

  ## Empirical Bayes scam:
  out_EB_scam <- empirical_bayes(disc_stats, method="scam")
  flb_EB_scam <- frac_sig_less_bias(out_EB_scam,ss$true_beta,i=4,alpha=5e-4)
  mse_EB_scam <- mse_sig_improve(out_EB_scam,ss$true_beta,i=4,alpha=5e-4)
  rmse_EB_scam <- mse_sig_improve_root(out_EB_scam,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_scam_up <- bias_sig_up(out_EB_scam,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_scam_down <- bias_sig_down(out_EB_scam,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_EB_scam <- mse_sig_improve_per(out_EB_scam,ss$true_beta,i=4,alpha=5e-4)

  ## Empirical Bayes gam_po:
  out_EB_gam_po <- empirical_bayes(disc_stats, method="gam_po")
  flb_EB_gam_po <- frac_sig_less_bias(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-4)
  mse_EB_gam_po <- mse_sig_improve(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-4)
  rmse_EB_gam_po <- mse_sig_improve_root(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_gam_po_up <- bias_sig_up(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_gam_po_down <- bias_sig_down(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_EB_gam_po <- mse_sig_improve_per(out_EB_gam_po,ss$true_beta,i=4,alpha=5e-4)

  ## Empirical Bayes gam_nb:
  out_EB_gam_nb <- empirical_bayes(disc_stats, method="gam_nb")
  flb_EB_gam_nb <- frac_sig_less_bias(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-4)
  mse_EB_gam_nb <- mse_sig_improve(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-4)
  rmse_EB_gam_nb <- mse_sig_improve_root(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_gam_nb_up <- bias_sig_up(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_gam_nb_down <- bias_sig_down(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_EB_gam_nb <- mse_sig_improve_per(out_EB_gam_nb,ss$true_beta,i=4,alpha=5e-4)

  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  mse_FIQT <- mse_sig_improve(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  rmse_FIQT <- mse_sig_improve_root(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  bias_FIQT_up <- bias_sig_up(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  bias_FIQT_down <- bias_sig_down(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_FIQT <- mse_sig_improve_per(out_FIQT,ss$true_beta,i=4,alpha=5e-4)

  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,ss$true_beta,i=4,alpha=5e-4)
  mse_BR <- mse_sig_improve(out_BR,ss$true_beta,i=4,alpha=5e-4)
  rmse_BR <- mse_sig_improve_root(out_BR,ss$true_beta,i=4,alpha=5e-4)
  bias_BR_up <- bias_sig_up(out_BR,ss$true_beta,i=4,alpha=5e-4)
  bias_BR_down <- bias_sig_down(out_BR,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_BR <- mse_sig_improve_per(out_BR,ss$true_beta,i=4,alpha=5e-4)

  ## cl1:
  out_cl <- conditional_likelihood_ind(disc_stats,alpha=5e-4)
  flb_cl1 <- frac_sig_less_bias(out_cl,ss$true_beta,i=4,alpha=5e-4)
  mse_cl1 <- mse_sig_improve(out_cl,ss$true_beta,i=4,alpha=5e-4)
  rmse_cl1 <- mse_sig_improve_root(out_cl,ss$true_beta,i=4,alpha=5e-4)
  bias_cl1_up <- bias_sig_up(out_cl,ss$true_beta,i=4,alpha=5e-4)
  bias_cl1_down <- bias_sig_down(out_cl,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_cl1 <- mse_sig_improve_per(out_cl,ss$true_beta,i=4,alpha=5e-4)

  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=5)
  mse_cl2 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=5)
  rmse_cl2 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=5)
  bias_cl2_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=5)
  bias_cl2_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=5)
  rel_mse_cl2 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=5)

  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=6)
  mse_cl3 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=6)
  rmse_cl3 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=6)
  bias_cl3_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=6)
  bias_cl3_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=6)
  rel_mse_cl3 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=6)

  ## replication
  ss2 <- data.frame(true_beta=ss$true_beta,se=ss$rep_se)
  rep_stats <- simulate_est(ss2)
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,i=4,alpha=5e-4)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,i=4,alpha=5e-4)
  rmse_rep <- mse_sig_improve_root(out_rep,ss$true_beta,i=4,alpha=5e-4)
  bias_rep_up <- bias_sig_up(out_rep,ss$true_beta,i=4,alpha=5e-4)
  bias_rep_down <- bias_sig_down(out_rep,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,i=4,alpha=5e-4)

  ## naive - only interested in bias (improvement should be 0):
  flb_naive <- frac_sig_less_bias(out_cl,disc_stats$true_beta,alpha=5e-4,i=2)
  mse_naive <- mse_sig_improve(out_cl,disc_stats$true_beta,alpha=5e-4,i=2)
  rmse_naive <- mse_sig_improve_root(out_cl,disc_stats$true_beta,alpha=5e-4,i=2)
  bias_naive_up <- bias_sig_up(out_cl,disc_stats$true_beta,alpha=5e-4,i=2)
  bias_naive_down <- bias_sig_down(out_cl,disc_stats$true_beta,alpha=5e-4,i=2)
  rel_mse_naive <- mse_sig_improve_per(out_cl,disc_stats$true_beta,alpha=5e-4,i=2)

  return(c(flb_EB,mse_EB,rmse_EB,bias_EB_up,bias_EB_down,rel_mse_EB,flb_EB_df,mse_EB_df,rmse_EB_df,bias_EB_df_up,bias_EB_df_down,rel_mse_EB_df,flb_EB_scam,mse_EB_scam,rmse_EB_scam,bias_EB_scam_up,bias_EB_scam_down,rel_mse_EB_scam,flb_EB_gam_po,mse_EB_gam_po,rmse_EB_gam_po,bias_EB_gam_po_up,bias_EB_gam_po_down,rel_mse_EB_gam_po,flb_EB_gam_nb,mse_EB_gam_nb,rmse_EB_gam_nb,bias_EB_gam_nb_up,bias_EB_gam_nb_down,rel_mse_EB_gam_nb,flb_FIQT,mse_FIQT,rmse_FIQT,bias_FIQT_up,bias_FIQT_down,rel_mse_FIQT,flb_BR,mse_BR,rmse_BR,bias_BR_up,bias_BR_down,rel_mse_BR,flb_cl1,mse_cl1,rmse_cl1,bias_cl1_up,bias_cl1_down,rel_mse_cl1,flb_cl2,mse_cl2,rmse_cl2,bias_cl2_up,bias_cl2_down,rel_mse_cl2,flb_cl3,mse_cl3,rmse_cl3,bias_cl3_up,bias_cl3_down,rel_mse_cl3,flb_rep,mse_rep,rmse_rep,bias_rep_up,bias_rep_down,rel_mse_rep,flb_naive,mse_naive,rmse_naive,bias_naive_up,bias_naive_down,rel_mse_naive))
}
res <- mclapply(1:nrow(sim_params), function(i){
  print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################


## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  bias_up[i] <- res[[i]][4]
  bias_down[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_EB <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)

## Empirical Bayes df=7:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rmse[i] <- res[[i]][9]
  bias_up[i] <- res[[i]][10]
  bias_down[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_EB_df <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB_df <- ave_results(res_EB_df,tot_sim)

## Empirical Bayes scam:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rmse[i] <- res[[i]][15]
  bias_up[i] <- res[[i]][16]
  bias_down[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_EB_scam <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB_scam <- ave_results(res_EB_scam,tot_sim)

## Empirical Bayes gam_po:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rmse[i] <- res[[i]][21]
  bias_up[i] <- res[[i]][22]
  bias_down[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_EB_gam_po <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB_gam_po <- ave_results(res_EB_gam_po,tot_sim)

## Empirical Bayes gam_nb:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rmse[i] <- res[[i]][27]
  bias_up[i] <- res[[i]][28]
  bias_down[i] <- res[[i]][29]
  rel_mse[i] <- res[[i]][30]
}
res_EB_gam_nb <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB_gam_nb <- ave_results(res_EB_gam_nb,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][31]
  mse[i] <- res[[i]][32]
  rmse[i] <- res[[i]][33]
  bias_up[i] <- res[[i]][34]
  bias_down[i] <- res[[i]][35]
  rel_mse[i] <- res[[i]][36]
}
res_FIQT <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## Bootstrap:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][37]
  mse[i] <- res[[i]][38]
  rmse[i] <- res[[i]][39]
  bias_up[i] <- res[[i]][40]
  bias_down[i] <- res[[i]][41]
  rel_mse[i] <- res[[i]][42]
}
res_BR <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)

## Conditional Likelihood 1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][43]
  mse[i] <- res[[i]][44]
  rmse[i] <- res[[i]][45]
  bias_up[i] <- res[[i]][46]
  bias_down[i] <- res[[i]][47]
  rel_mse[i] <- res[[i]][48]
}
res_cl1 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## Conditional Likelihood 2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][49]
  mse[i] <- res[[i]][50]
  rmse[i] <- res[[i]][51]
  bias_up[i] <- res[[i]][52]
  bias_down[i] <- res[[i]][53]
  rel_mse[i] <- res[[i]][54]
}
res_cl2 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## Conditional Likelihood 3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][55]
  mse[i] <- res[[i]][56]
  rmse[i] <- res[[i]][57]
  bias_up[i] <- res[[i]][58]
  bias_down[i] <- res[[i]][59]
  rel_mse[i] <- res[[i]][60]
}
res_cl3 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## replication
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][61]
  mse[i] <- res[[i]][62]
  rmse[i] <- res[[i]][63]
  bias_up[i] <- res[[i]][64]
  bias_down[i] <- res[[i]][65]
  rel_mse[i] <- res[[i]][66]
}
res_rep <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## Naive:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][67]
  mse[i] <- res[[i]][68]
  rmse[i] <- res[[i]][69]
  bias_up[i] <- res[[i]][70]
  bias_down[i] <- res[[i]][71]
  rel_mse[i] <- res[[i]][72]
}
res_naive <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_naive <- ave_results(res_naive,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_EB_df,ave_res_EB_scam,ave_res_EB_gam_po,ave_res_EB_gam_nb,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3, ave_res_rep, ave_res_naive)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("EB_df",(nrow(sim_params)/tot_sim)),rep("EB_scam",(nrow(sim_params)/tot_sim)),rep("EB_gam_po",(nrow(sim_params)/tot_sim)),rep("EB_gam_nb",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)), rep("rep",(nrow(sim_params)/tot_sim)),rep("naive",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/norm_5e-4_100sim_ave_2.csv")

## Combine all results:
results_all <- rbind(res_EB,res_EB_df,res_EB_scam,res_EB_gam_po,res_EB_gam_nb,res_FIQT,res_BR,res_cl1,res_cl2,res_cl3,res_rep,res_naive)
results_all$method <- c(rep("EB",nrow(sim_params)),rep("EB_df",nrow(sim_params)),rep("EB_scam",nrow(sim_params)),rep("EB_gam_po",nrow(sim_params)),rep("EB_gam_nb",nrow(sim_params)),rep("FIQT",nrow(sim_params)),rep("BR",nrow(sim_params)),rep("cl1",nrow(sim_params)),rep("cl2",nrow(sim_params)),rep("cl3",nrow(sim_params)),rep("rep",nrow(sim_params)),rep("naive",nrow(sim_params)))
write.csv(results_all,"results/norm_5e-4_100sim_all_2.csv")

print("PART 1B) complete!")

################################################################################
################################################################################
