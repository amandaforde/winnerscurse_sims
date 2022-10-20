## WINNER'S CURSE SIMULATION STUDY SCRIPT 1:
## Preliminary investigation with correlation structure

## This script obtains the number of significant SNPs, the proportion of these
## SNPs for which their association estimate is more extreme than their true
## effect size, the proportion of these SNPs which are significantly
## overexaggerated and the mean square error (MSE) of significant SNPs

## 1) Quantitative trait with normal effect size distribution

## Number of repetitions: 10
## Significance thresholds: alpha=5e-8
## Assumption: simple correlation structure imposed on independent blocks of 100
## SNPs

###############################################################################

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
  S = c(0)
)

################################################################################
################################################################################

## 1A) Quantitative trait with normal effect size distribution - 5e-8

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  out <- simulate_ss_ld(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  snp_sig <- out[abs(out$beta/out$se) > qnorm((5e-8)/2, lower.tail=FALSE),]
  n_sig <- nrow(snp_sig)
  if (n_sig == 0){
    prop_bias <- -1
    mse <- -1
    prop_x <- -1
  }else{
    prop_bias <- sum(abs(snp_sig$beta) > abs(out$true_beta[snp_sig$rsid]))/n_sig
    prop_x <- sum(abs(snp_sig$beta) > (abs(out$true_beta[snp_sig$rsid]) + 1.96*out$se[snp_sig$rsid]))/n_sig
    mse <- mean((out$true_beta[snp_sig$rsid]-snp_sig$beta)^2)
  }
  return(c(n_sig,prop_bias,prop_x,mse))
}
res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

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

write.csv(results, "results/nsig_prop_bias_5e_8_LD_all.csv")
write.csv(ave_res, "results/nsig_prop_bias_5e-8_LD.csv")

print("PART 1A) complete!")

################################################################################
################################################################################

## 1B) Quantitative trait with normal effect size distribution - 5e-4

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  out <- simulate_ss_ld(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  snp_sig <- out[abs(out$beta/out$se) > qnorm((5e-4)/2, lower.tail=FALSE),]
  n_sig <- nrow(snp_sig)
  if (n_sig == 0){
    prop_bias <- -1
    mse <- -1
    prop_x <- -1
  }else{
    prop_bias <- sum(abs(snp_sig$beta) > abs(out$true_beta[snp_sig$rsid]))/n_sig
    prop_x <- sum(abs(snp_sig$beta) > (abs(out$true_beta[snp_sig$rsid]) + 1.96*out$se[snp_sig$rsid]))/n_sig
    mse <- mean((out$true_beta[snp_sig$rsid]-snp_sig$beta)^2)
  }
  return(c(n_sig,prop_bias,prop_x,mse))
}
res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

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

write.csv(results, "results/nsig_prop_bias_5e_4_LD_all.csv")
write.csv(ave_res, "results/nsig_prop_bias_5e-4_LD.csv")

print("PART 1B) complete!")

################################################################################
################################################################################
