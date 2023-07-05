## WINNER'S CURSE SIMULATION STUDY SCRIPT PLUS:
## Bootstrap method evaluation and comparison with correlation structure


## 1) Quantitative trait with normal effect size distribution - simulate_ss_ld()

## Number of repetitions: 100 Significance thresholds: alpha=5e-8 and alpha=5e-4
## Methods evaluated: Bootstrap with and without adjustment to avoid upward
## correction
## Evaluation metrics: mse, rmse, rel_mse, bias_up, bias_down, flb, prop_nochange
## Assumption: simple correlation structure imposed on independent blocks of 100 SNPs

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


## Other form of bootstrap function to be tested:

BR_ss_2 <- function(summary_data,seed_opt = FALSE,seed=1998){
  stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
  stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
  stopifnot(nrow(summary_data) > 5)
  stopifnot(is.numeric(summary_data$beta) && is.numeric(summary_data$se))
  stopifnot(!any(duplicated(summary_data$rsid)))
  summary_data <- dplyr::arrange(summary_data, dplyr::desc((summary_data$beta/summary_data$se)))
  N <- nrow(summary_data)
  if(seed_opt==TRUE){set.seed(seed)}
  beta_boot <- matrix(stats::rnorm(1*N, mean = rep(summary_data$beta,1), sd = rep(summary_data$se,1)), nrow=N, ncol=1, byrow=FALSE)
  beta_mat <- matrix(rep(summary_data$beta,1), nrow=N, ncol=1, byrow=FALSE)
  se_mat <- matrix(rep(summary_data$se,1), nrow=N, ncol=1, byrow=FALSE)
  beta_oob <- beta_mat
  ordering <- apply(beta_boot/se_mat, 2,order,decreasing=TRUE)
  bias_correct <- matrix(nrow=N, ncol=1)
  bias_correct[,1] <- (beta_boot[ordering[,1],1] - beta_oob[ordering[,1],1])/summary_data$se[ordering[,1]]
  z <- summary_data$beta/summary_data$se
  bias_correct <- stats::predict(stats::smooth.spline(z,bias_correct)$fit, z)$y
  beta_BR_ss <- summary_data$beta - summary_data$se*bias_correct[rank(-1*summary_data$beta/summary_data$se)]
  beta_BR_ss[sign(beta_BR_ss) != sign(summary_data$beta)] <- 0
  summary_data <- cbind(summary_data, beta_BR_ss)
  ## following line removed in this version!
  ## for (i in 1:N){
  ##  if(abs(summary_data$beta[i]) < abs(summary_data$beta_BR_ss[i])){summary_data$beta_BR_ss[i] <- summary_data$beta[i]}
  ## }
  summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))
  return(summary_data)
}

## function to test the percentage of estimates that are unchanged
prop_nochange <- function(out,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm((alpha)/2, lower.tail=FALSE),]
  if (nrow(snp_sig) == 0){return(100)}
  prop_nochange <- sum(snp_sig$beta==snp_sig[,i])/length(snp_sig$rsid)
  return(prop_nochange)
}

## create ave_results2 function which includes prop_nochange
ave_results2 <- function(res_vec, n_sim){
  ave_flb <- c(rep(0,nrow(res_vec)/n_sim))
  sd_flb <- c(rep(0,nrow(res_vec)/n_sim))
  ave_mse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_mse <- c(rep(0,nrow(res_vec)/n_sim))
  ave_rmse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_rmse <- c(rep(0,nrow(res_vec)/n_sim))
  ave_bias_up <- c(rep(0,nrow(res_vec)/n_sim))
  sd_bias_up <- c(rep(0,nrow(res_vec)/n_sim))
  ave_bias_down <- c(rep(0,nrow(res_vec)/n_sim))
  sd_bias_down <- c(rep(0,nrow(res_vec)/n_sim))
  ave_rel_mse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_rel_mse <- c(rep(0,nrow(res_vec)/n_sim))
  ave_prop_nochange <- c(rep(0,nrow(res_vec)/n_sim))
  sd_prop_nochange <- c(rep(0,nrow(res_vec)/n_sim))
  for (i in 1:(nrow(res_vec)/n_sim)){
    ave_flb[i] <- average(res_vec$flb[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_flb[i] <- std(res_vec$flb[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_mse[i] <- average(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_mse[i] <- std(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_rmse[i] <- average(res_vec$rmse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_rmse[i] <- std(res_vec$rmse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_bias_up[i] <- average(res_vec$bias_up[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_bias_up[i] <- std(res_vec$bias_up[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_bias_down[i] <- average(res_vec$bias_down[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_bias_down[i] <- std(res_vec$bias_down[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_rel_mse[i] <- average(res_vec$rel_mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_rel_mse[i] <- std(res_vec$rel_mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_prop_nochange[i] <- average(res_vec$prop_nochange[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_prop_nochange[i] <- std(res_vec$prop_nochange[(i*n_sim-(n_sim-1)):(i*n_sim)])
  }
  keep <- seq(from=1,to=nrow(res_vec)-n_sim+1,by=n_sim)
  res_vec_ave <- res_vec[rownames(res_vec) %in% keep,]
  res_vec_ave$flb <- ave_flb
  res_vec_ave$flb_error <- sd_flb
  res_vec_ave$mse <- ave_mse
  res_vec_ave$mse_error <- sd_mse
  res_vec_ave$rmse <- ave_rmse
  res_vec_ave$rmse_error <- sd_rmse
  res_vec_ave$bias_up <- ave_bias_up
  res_vec_ave$bias_up_error <- sd_bias_up
  res_vec_ave$bias_down <- ave_bias_down
  res_vec_ave$bias_down_error <- sd_bias_down
  res_vec_ave$rel_mse <- ave_rel_mse
  res_vec_ave$rel_mse_error <- sd_rel_mse
  res_vec_ave$prop_nochange <- ave_prop_nochange
  res_vec_ave$prop_nochange_error <- sd_prop_nochange
  res_vec_ave <- subset(res_vec_ave, select = -sim )
  return(res_vec_ave)
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

  ## Bootstrap 1:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,disc_stats$true_beta,alpha=5e-8)
  mse_BR <- mse_sig_improve(out_BR,disc_stats$true_beta,alpha=5e-8)
  rmse_BR <- mse_sig_improve_root(out_BR,disc_stats$true_beta,alpha=5e-8)
  bias_BR_up <- bias_sig_up(out_BR,disc_stats$true_beta,alpha=5e-8)
  bias_BR_down <- bias_sig_down(out_BR,disc_stats$true_beta,alpha=5e-8)
  rel_mse_BR <- mse_sig_improve_per(out_BR,disc_stats$true_beta,alpha=5e-8)
  prop_nochange_BR <- prop_nochange(out_BR, alpha=5e-8)

  ## Bootstrap 2:
  out_BR2 <- BR_ss_2(disc_stats)
  flb_BR2 <- frac_sig_less_bias(out_BR2,disc_stats$true_beta,alpha=5e-8)
  mse_BR2 <- mse_sig_improve(out_BR2,disc_stats$true_beta,alpha=5e-8)
  rmse_BR2 <- mse_sig_improve_root(out_BR2,disc_stats$true_beta,alpha=5e-8)
  bias_BR2_up <- bias_sig_up(out_BR2,disc_stats$true_beta,alpha=5e-8)
  bias_BR2_down <- bias_sig_down(out_BR2,disc_stats$true_beta,alpha=5e-8)
  rel_mse_BR2 <- mse_sig_improve_per(out_BR2,disc_stats$true_beta,alpha=5e-8)
  prop_nochange_BR2 <- prop_nochange(out_BR2, alpha=5e-8)


  ## naive - only interested in bias (improvement should be 0):
  flb_naive <- frac_sig_less_bias(out_BR,disc_stats$true_beta,alpha=5e-8,i=2)
  mse_naive <- mse_sig_improve(out_BR,disc_stats$true_beta,alpha=5e-8,i=2)
  rmse_naive <- mse_sig_improve_root(out_BR,disc_stats$true_beta,alpha=5e-8,i=2)
  bias_naive_up <- bias_sig_up(out_BR,disc_stats$true_beta,alpha=5e-8,i=2)
  bias_naive_down <- bias_sig_down(out_BR,disc_stats$true_beta,alpha=5e-8,i=2)
  rel_mse_naive <- mse_sig_improve_per(out_BR,disc_stats$true_beta,alpha=5e-8,i=2)
  prop_nochange_naive <- prop_nochange(out_BR, alpha=5e-8,i=2)

  return(c(flb_BR,mse_BR,rmse_BR,bias_BR_up,bias_BR_down,rel_mse_BR,prop_nochange_BR,flb_BR2,mse_BR2,rmse_BR2,bias_BR2_up,bias_BR2_down,rel_mse_BR2,prop_nochange_BR2,flb_naive,mse_naive,rmse_naive,bias_naive_up,bias_naive_down,rel_mse_naive,prop_nochange_naive))
}
res <- mclapply(1:nrow(sim_params), function(i){
  print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################

## Organising results:

## Bootstrap:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
prop_nochange <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  bias_up[i] <- res[[i]][4]
  bias_down[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
  prop_nochange[i] <- res[[i]][7]
}
res_BR <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse,prop_nochange)
ave_res_BR <- ave_results2(res_BR,tot_sim)

## Bootstrap 2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
prop_nochange <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][8]
  mse[i] <- res[[i]][9]
  rmse[i] <- res[[i]][10]
  bias_up[i] <- res[[i]][11]
  bias_down[i] <- res[[i]][12]
  rel_mse[i] <- res[[i]][13]
  prop_nochange[i] <- res[[i]][14]
}
res_BR2 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse,prop_nochange)
ave_res_BR2 <- ave_results2(res_BR2,tot_sim)

## Naive:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
prop_nochange <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][15]
  mse[i] <- res[[i]][16]
  rmse[i] <- res[[i]][17]
  bias_up[i] <- res[[i]][18]
  bias_down[i] <- res[[i]][19]
  rel_mse[i] <- res[[i]][20]
  prop_nochange[i] <- res[[i]][21]
}
res_naive <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse,prop_nochange)
ave_res_naive <- ave_results2(res_naive,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_BR,ave_res_BR2,ave_res_naive)
results_all$method <- c(rep("BR",(nrow(sim_params)/tot_sim)),rep("BR2",(nrow(sim_params)/tot_sim)),rep("naive",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/norm_5e-8_100sim_LD_ave_boot.csv")

## Combine all results:
results_all <- rbind(res_BR,res_BR2,res_naive)
results_all$method <- c(rep("BR",nrow(sim_params)),rep("BR2",nrow(sim_params)),rep("naive",nrow(sim_params)))
write.csv(results_all,"results/norm_5e-8_100sim_LD_all_boot.csv")

print("PART BOOT A) complete!")

################################################################################
################################################################################

## 1B) Quantitative trait with normal effect size distribution - 5e-4

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  disc_stats <- simulate_ss_ld(H=h2,Pi=prop_effect,nid=n_samples,sc=S)

  ## Bootstrap 1:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,disc_stats$true_beta,alpha=5e-4)
  mse_BR <- mse_sig_improve(out_BR,disc_stats$true_beta,alpha=5e-4)
  rmse_BR <- mse_sig_improve_root(out_BR,disc_stats$true_beta,alpha=5e-4)
  bias_BR_up <- bias_sig_up(out_BR,disc_stats$true_beta,alpha=5e-4)
  bias_BR_down <- bias_sig_down(out_BR,disc_stats$true_beta,alpha=5e-4)
  rel_mse_BR <- mse_sig_improve_per(out_BR,disc_stats$true_beta,alpha=5e-4)
  prop_nochange_BR <- prop_nochange(out_BR, alpha=5e-4)

  ## Bootstrap 2:
  out_BR2 <- BR_ss_2(disc_stats)
  flb_BR2 <- frac_sig_less_bias(out_BR2,disc_stats$true_beta,alpha=5e-4)
  mse_BR2 <- mse_sig_improve(out_BR2,disc_stats$true_beta,alpha=5e-4)
  rmse_BR2 <- mse_sig_improve_root(out_BR2,disc_stats$true_beta,alpha=5e-4)
  bias_BR2_up <- bias_sig_up(out_BR2,disc_stats$true_beta,alpha=5e-4)
  bias_BR2_down <- bias_sig_down(out_BR2,disc_stats$true_beta,alpha=5e-4)
  rel_mse_BR2 <- mse_sig_improve_per(out_BR2,disc_stats$true_beta,alpha=5e-4)
  prop_nochange_BR2 <- prop_nochange(out_BR2, alpha=5e-4)


  ## naive - only interested in bias (improvement should be 0):
  flb_naive <- frac_sig_less_bias(out_BR,disc_stats$true_beta,alpha=5e-4,i=2)
  mse_naive <- mse_sig_improve(out_BR,disc_stats$true_beta,alpha=5e-4,i=2)
  rmse_naive <- mse_sig_improve_root(out_BR,disc_stats$true_beta,alpha=5e-4,i=2)
  bias_naive_up <- bias_sig_up(out_BR,disc_stats$true_beta,alpha=5e-4,i=2)
  bias_naive_down <- bias_sig_down(out_BR,disc_stats$true_beta,alpha=5e-4,i=2)
  rel_mse_naive <- mse_sig_improve_per(out_BR,disc_stats$true_beta,alpha=5e-4,i=2)
  prop_nochange_naive <- prop_nochange(out_BR, alpha=5e-4,i=2)

  return(c(flb_BR,mse_BR,rmse_BR,bias_BR_up,bias_BR_down,rel_mse_BR,prop_nochange_BR,flb_BR2,mse_BR2,rmse_BR2,bias_BR2_up,bias_BR2_down,rel_mse_BR2,prop_nochange_BR2,flb_naive,mse_naive,rmse_naive,bias_naive_up,bias_naive_down,rel_mse_naive,prop_nochange_naive))
}
res <- mclapply(1:nrow(sim_params), function(i){
  print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################

## Organising results:

## Bootstrap:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
prop_nochange <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  bias_up[i] <- res[[i]][4]
  bias_down[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
  prop_nochange[i] <- res[[i]][7]
}
res_BR <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse,prop_nochange)
ave_res_BR <- ave_results2(res_BR,tot_sim)

## Bootstrap 2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
prop_nochange <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][8]
  mse[i] <- res[[i]][9]
  rmse[i] <- res[[i]][10]
  bias_up[i] <- res[[i]][11]
  bias_down[i] <- res[[i]][12]
  rel_mse[i] <- res[[i]][13]
  prop_nochange[i] <- res[[i]][14]
}
res_BR2 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse,prop_nochange)
ave_res_BR2 <- ave_results2(res_BR2,tot_sim)

## Naive:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
prop_nochange <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][15]
  mse[i] <- res[[i]][16]
  rmse[i] <- res[[i]][17]
  bias_up[i] <- res[[i]][18]
  bias_down[i] <- res[[i]][19]
  rel_mse[i] <- res[[i]][20]
  prop_nochange[i] <- res[[i]][21]
}
res_naive <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse,prop_nochange)
ave_res_naive <- ave_results2(res_naive,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_BR,ave_res_BR2,ave_res_naive)
results_all$method <- c(rep("BR",(nrow(sim_params)/tot_sim)),rep("BR2",(nrow(sim_params)/tot_sim)),rep("naive",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/norm_5e-4_100sim_LD_ave_boot.csv")

## Combine all results:
results_all <- rbind(res_BR,res_BR2,res_naive)
results_all$method <- c(rep("BR",nrow(sim_params)),rep("BR2",nrow(sim_params)),rep("naive",nrow(sim_params)))
write.csv(results_all,"results/norm_5e-4_100sim_LD_all_boot.csv")

print("PART BOOT B) complete!")




