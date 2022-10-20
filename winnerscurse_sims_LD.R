## WINNER'S CURSE SIMULATION STUDY SCRIPT 2:
## Method evaluation and comparison with correlation structure


## 1) Quantitative trait with normal effect size distribution - simulate_ss_ld()

## Number of repetitions: 100
## Significance thresholds: alpha=5e-8 and alpha=5e-4
## Methods evaluated: Empirical Bayes, Empirical Bayes with fixed df-7, Empirical
## Bayes using scam, Empirical Bayes using gam and poisson, Empirical Bayes using
## gam and negative binomial, FIQT, Bootstrap, Conditional Likelihood methods
## Evaluation metrics: mse, rmse, rel_mse, flb
## Assumption: simple correlation structure imposed on independent blocks of 100
## SNPs

################################################################################

## 1) Extra functions required to simulate correlation structure:

x <- 0.9825
s <- function(n,x=0.9825){
  vector <- c()
  for (i in 1:n){
    vector[i] <- x^(i-1)
  }
  return(vector)
}

vec <- function(tot){
  vector <- c()
  for (i in 1:tot){
    vector <- c(vector,s(100-(i-1)))
  }
  return(vector)
}

m1 <- matrix(NA, 100, 100)
m1[lower.tri(m1, diag=TRUE)] <- vec(100)
m2 <- t(m1)
m2[lower.tri(m2, diag=TRUE)] <- vec(100)
R <- m2
R_sqrt <- sqrtm(R) ## remains constant


## Function to simulate LD - takes about 30 seconds to simulate 1 set of
## summary statistics:
simulate_ss_ld <- function(H,Pi,nid,sc,cormat=R,cormatsq=R_sqrt){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- rnorm(effect_snps,0,sd=sqrt((2*maf[1:effect_snps]*(1-maf[1:effect_snps]))^sc))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/H
  true_beta <- true_beta/sqrt(var_y)
  true_beta_all <- c(rep(0,n_snps))
  maf_re <- c(rep(0,n_snps))

  positions <- sample(1:n_snps, effect_snps, replace=F)
  for (i in 1:effect_snps){
    true_beta_all[positions[i]] <- true_beta[i]
    maf_re[positions[i]] <- maf[i]
  }

  vec <- seq(1,10^6)
  vec_new <- vec[! vec %in% positions]
  for (i in (effect_snps+1):n_snps){
    maf_re[vec_new[i-effect_snps]] <- maf[i]
  }

  mu_all <- numeric(n_snps)
  beta_hat_all <- numeric(n_snps)
  se_all <- numeric(n_snps)

  ## R defined outside of function as it is fixed
  for (i in 1:10000){
    D_1 <- diag(x = 1/(sqrt(nid*2*maf_re[((i-1)*100+1):(i*100)]*(1-maf_re[((i-1)*100+1):(i*100)]))), nrow=100, ncol=100)
    D_2 <- diag(x = (sqrt(nid*2*maf_re[((i-1)*100+1):(i*100)]*(1-maf_re[((i-1)*100+1):(i*100)]))), nrow=100, ncol=100)
    A <- D_1 %*% cormat %*% D_2
    B <- cormatsq %*% D_1
    mu_all[((i-1)*100+1):(i*100)] <- A %*% true_beta_all[((i-1)*100+1):(i*100)]
    beta_hat_all[((i-1)*100+1):(i*100)] <- mu_all[((i-1)*100+1):(i*100)] + ( B %*% rnorm(100))
    se_all[((i-1)*100+1):(i*100)] <- diag(D_1)
  }

  summary_stats <- data.frame(rsid=seq(1,n_snps),beta=beta_hat_all,se=se_all,true_beta=mu_all)
  return(summary_stats)
}

################################################################################
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
  prop_effect = c(0.01, 0.001),
  S = c(0)
)

################################################################################
################################################################################

## 1A) Quantitative trait with normal effect size distribution - 5e-8

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  disc_stats <- simulate_ss_ld(H=h2,Pi=prop_effect,nid=n_samples,sc=S)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,disc_stats$true_beta,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,disc_stats$true_beta,alpha=5e-8)
  rmse_EB <- mse_sig_improve_root(out_EB,disc_stats$true_beta,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,disc_stats$true_beta,alpha=5e-8)

  ## Empirical Bayes df=7:
  out_EB_df <- empirical_bayes(disc_stats, method="fix_df")
  flb_EB_df <- frac_sig_less_bias(out_EB_df,disc_stats$true_beta,alpha=5e-8)
  mse_EB_df <- mse_sig_improve(out_EB_df,disc_stats$true_beta,alpha=5e-8)
  rmse_EB_df <- mse_sig_improve_root(out_EB_df,disc_stats$true_beta,alpha=5e-8)
  rel_mse_EB_df <- mse_sig_improve_per(out_EB_df,disc_stats$true_beta,alpha=5e-8)

  ## Empirical Bayes scam:
  out_EB_scam <- empirical_bayes(disc_stats, method="scam")
  flb_EB_scam <- frac_sig_less_bias(out_EB_scam,disc_stats$true_beta,alpha=5e-8)
  mse_EB_scam <- mse_sig_improve(out_EB_scam,disc_stats$true_beta,alpha=5e-8)
  rmse_EB_scam <- mse_sig_improve_root(out_EB_scam,disc_stats$true_beta,alpha=5e-8)
  rel_mse_EB_scam <- mse_sig_improve_per(out_EB_scam,disc_stats$true_beta,alpha=5e-8)

  ## Empirical Bayes gam_po:
  out_EB_gam_po <- empirical_bayes(disc_stats, method="gam_po")
  flb_EB_gam_po <- frac_sig_less_bias(out_EB_gam_po,disc_stats$true_beta,alpha=5e-8)
  mse_EB_gam_po <- mse_sig_improve(out_EB_gam_po,disc_stats$true_beta,alpha=5e-8)
  rmse_EB_gam_po <- mse_sig_improve_root(out_EB_gam_po,disc_stats$true_beta,alpha=5e-8)
  rel_mse_EB_gam_po <- mse_sig_improve_per(out_EB_gam_po,disc_stats$true_beta,alpha=5e-8)

  ## Empirical Bayes gam_nb:
  out_EB_gam_nb <- empirical_bayes(disc_stats, method="gam_nb")
  flb_EB_gam_nb <- frac_sig_less_bias(out_EB_gam_nb,disc_stats$true_beta,alpha=5e-8)
  mse_EB_gam_nb <- mse_sig_improve(out_EB_gam_nb,disc_stats$true_beta,alpha=5e-8)
  rmse_EB_gam_nb <- mse_sig_improve_root(out_EB_gam_nb,disc_stats$true_beta,alpha=5e-8)
  rel_mse_EB_gam_nb <- mse_sig_improve_per(out_EB_gam_nb,disc_stats$true_beta,alpha=5e-8)

  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,disc_stats$true_beta,alpha=5e-8)
  mse_FIQT <- mse_sig_improve(out_FIQT,disc_stats$true_beta,alpha=5e-8)
  rmse_FIQT <- mse_sig_improve_root(out_FIQT,disc_stats$true_beta,alpha=5e-8)
  rel_mse_FIQT <- mse_sig_improve_per(out_FIQT,disc_stats$true_beta,alpha=5e-8)

  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,disc_stats$true_beta,alpha=5e-8)
  mse_BR <- mse_sig_improve(out_BR,disc_stats$true_beta,alpha=5e-8)
  rmse_BR <- mse_sig_improve_root(out_BR,disc_stats$true_beta,alpha=5e-8)
  rel_mse_BR <- mse_sig_improve_per(out_BR,disc_stats$true_beta,alpha=5e-8)

  ## cl1:
  out_cl <- conditional_likelihood_old(disc_stats,alpha=5e-8)
  flb_cl1 <- frac_sig_less_bias(out_cl,disc_stats$true_beta,alpha=5e-8)
  mse_cl1 <- mse_sig_improve(out_cl,disc_stats$true_beta,alpha=5e-8)
  rmse_cl1 <- mse_sig_improve_root(out_cl,disc_stats$true_beta,alpha=5e-8)
  rel_mse_cl1 <- mse_sig_improve_per(out_cl,disc_stats$true_beta,alpha=5e-8)

  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,disc_stats$true_beta,alpha=5e-8,i=6)
  mse_cl2 <- mse_sig_improve(out_cl,disc_stats$true_beta,alpha=5e-8,i=6)
  rmse_cl2 <- mse_sig_improve_root(out_cl,disc_stats$true_beta,alpha=5e-8,i=6)
  rel_mse_cl2 <- mse_sig_improve_per(out_cl,disc_stats$true_beta,alpha=5e-8,i=6)

  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,disc_stats$true_beta,alpha=5e-8,i=7)
  mse_cl3 <- mse_sig_improve(out_cl,disc_stats$true_beta,alpha=5e-8,i=7)
  rmse_cl3 <- mse_sig_improve_root(out_cl,disc_stats$true_beta,alpha=5e-8,i=7)
  rel_mse_cl3 <- mse_sig_improve_per(out_cl,disc_stats$true_beta,alpha=5e-8,i=7)

  return(c(flb_EB,mse_EB,rmse_EB,rel_mse_EB,flb_EB_df,mse_EB_df,rmse_EB_df,rel_mse_EB_df,flb_EB_scam,mse_EB_scam,rmse_EB_scam,rel_mse_EB_scam,flb_EB_gam_po,mse_EB_gam_po,rmse_EB_gam_po,rel_mse_EB_gam_po,flb_EB_gam_nb,mse_EB_gam_nb,rmse_EB_gam_nb,rel_mse_EB_gam_nb,flb_FIQT,mse_FIQT,rmse_FIQT,rel_mse_FIQT,flb_BR,mse_BR,rmse_BR,rel_mse_BR,flb_cl1,mse_cl1,rmse_cl1,rel_mse_cl1,flb_cl2,mse_cl2,rmse_cl2,rel_mse_cl2,flb_cl3,mse_cl3,rmse_cl3,rel_mse_cl3))
}
res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################

## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  rel_mse[i] <- res[[i]][4]
}
res_EB <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)

## Empirical Bayes df=7:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][5]
  mse[i] <- res[[i]][6]
  rmse[i] <- res[[i]][7]
  rel_mse[i] <- res[[i]][8]
}
res_EB_df <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_EB_df <- ave_results(res_EB_df,tot_sim)

## Empirical Bayes scam:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][9]
  mse[i] <- res[[i]][10]
  rmse[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_EB_scam <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_EB_scam <- ave_results(res_EB_scam,tot_sim)

## Empirical Bayes gam_po:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rmse[i] <- res[[i]][15]
  rel_mse[i] <- res[[i]][16]
}
res_EB_gam_po <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_EB_gam_po <- ave_results(res_EB_gam_po,tot_sim)

## Empirical Bayes gam_nb:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][17]
  mse[i] <- res[[i]][18]
  rmse[i] <- res[[i]][19]
  rel_mse[i] <- res[[i]][20]
}
res_EB_gam_nb <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_EB_gam_nb <- ave_results(res_EB_gam_nb,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][21]
  mse[i] <- res[[i]][22]
  rmse[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_FIQT <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## Bootstrap:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rmse[i] <- res[[i]][27]
  rel_mse[i] <- res[[i]][28]
}
res_BR <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)

## Conditional Likelihood 1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][29]
  mse[i] <- res[[i]][30]
  rmse[i] <- res[[i]][31]
  rel_mse[i] <- res[[i]][32]
}
res_cl1 <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## Conditional Likelihood 2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][33]
  mse[i] <- res[[i]][34]
  rmse[i] <- res[[i]][35]
  rel_mse[i] <- res[[i]][36]
}
res_cl2 <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## Conditional Likelihood 3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][37]
  mse[i] <- res[[i]][38]
  rmse[i] <- res[[i]][39]
  rel_mse[i] <- res[[i]][40]
}
res_cl3 <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_EB_df,ave_res_EB_scam,ave_res_EB_gam_po,ave_res_EB_gam_nb,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("EB_df",(nrow(sim_params)/tot_sim)),rep("EB_scam",(nrow(sim_params)/tot_sim)),rep("EB_gam_po",(nrow(sim_params)/tot_sim)),rep("EB_gam_nb",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/norm_5e-8_100sim_LD_ave.csv")

## Combine all results:
results_all <- rbind(res_EB,res_EB_df,res_EB_scam,res_EB_gam_po,res_EB_gam_nb,res_FIQT,res_BR,res_cl1,res_cl2,res_cl3)
results_all$method <- c(rep("EB",nrow(sim_params)),rep("EB_df",nrow(sim_params)),rep("EB_scam",nrow(sim_params)),rep("EB_gam_po",nrow(sim_params)),rep("EB_gam_nb",nrow(sim_params)),rep("FIQT",nrow(sim_params)),rep("BR",nrow(sim_params)),rep("cl1",nrow(sim_params)),rep("cl2",nrow(sim_params)),rep("cl3",nrow(sim_params)))
write.csv(results_all,"results/norm_5e-8_100sim_LD_all.csv")

print("PART 1A) complete!")

################################################################################
################################################################################

## 1B) Quantitative trait with normal effect size distribution - 5e-4

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  disc_stats <- simulate_ss_ld(H=h2,Pi=prop_effect,nid=n_samples,sc=S)

  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,disc_stats$true_beta,alpha=5e-4)
  mse_EB <- mse_sig_improve(out_EB,disc_stats$true_beta,alpha=5e-4)
  rmse_EB <- mse_sig_improve_root(out_EB,disc_stats$true_beta,alpha=5e-4)
  rel_mse_EB <- mse_sig_improve_per(out_EB,disc_stats$true_beta,alpha=5e-4)

  ## Empirical Bayes df=7:
  out_EB_df <- empirical_bayes(disc_stats, method="fix_df")
  flb_EB_df <- frac_sig_less_bias(out_EB_df,disc_stats$true_beta,alpha=5e-4)
  mse_EB_df <- mse_sig_improve(out_EB_df,disc_stats$true_beta,alpha=5e-4)
  rmse_EB_df <- mse_sig_improve_root(out_EB_df,disc_stats$true_beta,alpha=5e-4)
  rel_mse_EB_df <- mse_sig_improve_per(out_EB_df,disc_stats$true_beta,alpha=5e-4)

  ## Empirical Bayes scam:
  out_EB_scam <- empirical_bayes(disc_stats, method="scam")
  flb_EB_scam <- frac_sig_less_bias(out_EB_scam,disc_stats$true_beta,alpha=5e-4)
  mse_EB_scam <- mse_sig_improve(out_EB_scam,disc_stats$true_beta,alpha=5e-4)
  rmse_EB_scam <- mse_sig_improve_root(out_EB_scam,disc_stats$true_beta,alpha=5e-4)
  rel_mse_EB_scam <- mse_sig_improve_per(out_EB_scam,disc_stats$true_beta,alpha=5e-4)

  ## Empirical Bayes gam_po:
  out_EB_gam_po <- empirical_bayes(disc_stats, method="gam_po")
  flb_EB_gam_po <- frac_sig_less_bias(out_EB_gam_po,disc_stats$true_beta,alpha=5e-4)
  mse_EB_gam_po <- mse_sig_improve(out_EB_gam_po,disc_stats$true_beta,alpha=5e-4)
  rmse_EB_gam_po <- mse_sig_improve_root(out_EB_gam_po,disc_stats$true_beta,alpha=5e-4)
  rel_mse_EB_gam_po <- mse_sig_improve_per(out_EB_gam_po,disc_stats$true_beta,alpha=5e-4)

  ## Empirical Bayes gam_nb:
  out_EB_gam_nb <- empirical_bayes(disc_stats, method="gam_nb")
  flb_EB_gam_nb <- frac_sig_less_bias(out_EB_gam_nb,disc_stats$true_beta,alpha=5e-4)
  mse_EB_gam_nb <- mse_sig_improve(out_EB_gam_nb,disc_stats$true_beta,alpha=5e-4)
  rmse_EB_gam_nb <- mse_sig_improve_root(out_EB_gam_nb,disc_stats$true_beta,alpha=5e-4)
  rel_mse_EB_gam_nb <- mse_sig_improve_per(out_EB_gam_nb,disc_stats$true_beta,alpha=5e-4)

  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,disc_stats$true_beta,alpha=5e-4)
  mse_FIQT <- mse_sig_improve(out_FIQT,disc_stats$true_beta,alpha=5e-4)
  rmse_FIQT <- mse_sig_improve_root(out_FIQT,disc_stats$true_beta,alpha=5e-4)
  rel_mse_FIQT <- mse_sig_improve_per(out_FIQT,disc_stats$true_beta,alpha=5e-4)

  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,disc_stats$true_beta,alpha=5e-4)
  mse_BR <- mse_sig_improve(out_BR,disc_stats$true_beta,alpha=5e-4)
  rmse_BR <- mse_sig_improve_root(out_BR,disc_stats$true_beta,alpha=5e-4)
  rel_mse_BR <- mse_sig_improve_per(out_BR,disc_stats$true_beta,alpha=5e-4)

  ## cl1:
  out_cl <- conditional_likelihood_old(disc_stats,alpha=5e-4)
  flb_cl1 <- frac_sig_less_bias(out_cl,disc_stats$true_beta,alpha=5e-4)
  mse_cl1 <- mse_sig_improve(out_cl,disc_stats$true_beta,alpha=5e-4)
  rmse_cl1 <- mse_sig_improve_root(out_cl,disc_stats$true_beta,alpha=5e-4)
  rel_mse_cl1 <- mse_sig_improve_per(out_cl,disc_stats$true_beta,alpha=5e-4)

  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,disc_stats$true_beta,alpha=5e-4,i=6)
  mse_cl2 <- mse_sig_improve(out_cl,disc_stats$true_beta,alpha=5e-4,i=6)
  rmse_cl2 <- mse_sig_improve_root(out_cl,disc_stats$true_beta,alpha=5e-4,i=6)
  rel_mse_cl2 <- mse_sig_improve_per(out_cl,disc_stats$true_beta,alpha=5e-4,i=6)

  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,disc_stats$true_beta,alpha=5e-4,i=7)
  mse_cl3 <- mse_sig_improve(out_cl,disc_stats$true_beta,alpha=5e-4,i=7)
  rmse_cl3 <- mse_sig_improve_root(out_cl,disc_stats$true_beta,alpha=5e-4,i=7)
  rel_mse_cl3 <- mse_sig_improve_per(out_cl,disc_stats$true_beta,alpha=5e-4,i=7)

  return(c(flb_EB,mse_EB,rmse_EB,rel_mse_EB,flb_EB_df,mse_EB_df,rmse_EB_df,rel_mse_EB_df,flb_EB_scam,mse_EB_scam,rmse_EB_scam,rel_mse_EB_scam,flb_EB_gam_po,mse_EB_gam_po,rmse_EB_gam_po,rel_mse_EB_gam_po,flb_EB_gam_nb,mse_EB_gam_nb,rmse_EB_gam_nb,rel_mse_EB_gam_nb,flb_FIQT,mse_FIQT,rmse_FIQT,rel_mse_FIQT,flb_BR,mse_BR,rmse_BR,rel_mse_BR,flb_cl1,mse_cl1,rmse_cl1,rel_mse_cl1,flb_cl2,mse_cl2,rmse_cl2,rel_mse_cl2,flb_cl3,mse_cl3,rmse_cl3,rel_mse_cl3))
}

res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################

## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  rel_mse[i] <- res[[i]][4]
}
res_EB <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)

## Empirical Bayes df=7:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][5]
  mse[i] <- res[[i]][6]
  rmse[i] <- res[[i]][7]
  rel_mse[i] <- res[[i]][8]
}
res_EB_df <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_EB_df <- ave_results(res_EB_df,tot_sim)

## Empirical Bayes scam:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][9]
  mse[i] <- res[[i]][10]
  rmse[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_EB_scam <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_EB_scam <- ave_results(res_EB_scam,tot_sim)


## Empirical Bayes gam_po:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rmse[i] <- res[[i]][15]
  rel_mse[i] <- res[[i]][16]
}
res_EB_gam_po <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_EB_gam_po <- ave_results(res_EB_gam_po,tot_sim)

## Empirical Bayes gam_nb:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][17]
  mse[i] <- res[[i]][18]
  rmse[i] <- res[[i]][19]
  rel_mse[i] <- res[[i]][20]
}
res_EB_gam_nb <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_EB_gam_nb <- ave_results(res_EB_gam_nb,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][21]
  mse[i] <- res[[i]][22]
  rmse[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_FIQT <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## Bootstrap:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rmse[i] <- res[[i]][27]
  rel_mse[i] <- res[[i]][28]
}
res_BR <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)

## Conditional Likelihood 1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][29]
  mse[i] <- res[[i]][30]
  rmse[i] <- res[[i]][31]
  rel_mse[i] <- res[[i]][32]
}
res_cl1 <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## Conditional Likelihood 2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][33]
  mse[i] <- res[[i]][34]
  rmse[i] <- res[[i]][35]
  rel_mse[i] <- res[[i]][36]
}
res_cl2 <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## Conditional Likelihood 3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][37]
  mse[i] <- res[[i]][38]
  rmse[i] <- res[[i]][39]
  rel_mse[i] <- res[[i]][40]
}
res_cl3 <- cbind(sim_params,flb,mse,rmse,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_EB_df,ave_res_EB_scam,ave_res_EB_gam_po,ave_res_EB_gam_nb,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("EB_df",(nrow(sim_params)/tot_sim)),rep("EB_scam",(nrow(sim_params)/tot_sim)),rep("EB_gam_po",(nrow(sim_params)/tot_sim)),rep("EB_gam_nb",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/norm_5e-4_100sim_LD_ave.csv")

## Combine all results:
results_all <- rbind(res_EB,res_EB_df,res_EB_scam,res_EB_gam_po,res_EB_gam_nb,res_FIQT,res_BR,res_cl1,res_cl2,res_cl3)
results_all$method <- c(rep("EB",nrow(sim_params)),rep("EB_df",nrow(sim_params)),rep("EB_scam",nrow(sim_params)),rep("EB_gam_po",nrow(sim_params)),rep("EB_gam_nb",nrow(sim_params)),rep("FIQT",nrow(sim_params)),rep("BR",nrow(sim_params)),rep("cl1",nrow(sim_params)),rep("cl2",nrow(sim_params)),rep("cl3",nrow(sim_params)))
write.csv(results_all,"results/norm_5e-4_100sim_LD_all.csv")

print("PART 1B) complete!")

################################################################################
################################################################################
