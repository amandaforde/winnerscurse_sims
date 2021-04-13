## FUNCTIONS REQUIRED FOR SIMULATIONS:

## 1. Given values for heritability (H), polygenicity (Pi), sample size (nid)
## and selection coefficient (S), 'true' values of effect size and corresponding
## standard error are simulated in which effect sizes are assumed to be normally
## distributed
simulate_ss <- function(H,Pi,nid,sc){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- rnorm(effect_snps,0,sd=sqrt((2*maf*(1-maf))^sc))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/H
  true_beta <- true_beta/sqrt(var_y)
  true_beta <- c(true_beta, rep(0,n_snps-effect_snps))
  se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(nid-2)*maf*(1-maf)))
  stats <- data.frame(true_beta,se)
  return(stats)
}


## 2. Given values for heritability (H), polygenicity (Pi), sample size (nid)
## and selection coefficient (S), 'true' values of effect size and corresponding
## standard error are simulated in which effect sizes are assumed to have a
## skewed distribution, i.e. 50% of effect sizes come from a normal distribution
## centered at 0 while the other half are generated from a normal distribution
## with mean 2.5
simulate_ss_skew <- function(H,Pi,nid,sc){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- c(rnorm(0.5*effect_snps,2.5,sd=sqrt((2*maf*(1-maf))^sc)),rnorm(0.5*effect_snps,0,sd=sqrt((2*maf[(0.5*effect_snps+1):effect_snps]*(1-maf[(0.5*effect_snps+1):effect_snps]))^sc)))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/H
  true_beta <- true_beta/sqrt(var_y)
  true_beta <- c(true_beta, rep(0,n_snps-effect_snps))
  se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(nid-2)*maf*(1-maf)))
  stats <- data.frame(true_beta,se)
  return(stats)
}

## 3. Simulates discovery GWAS summary statistics given an output from
## functions 1 or 2
simulate_est <- function(stats){
  est <- data.frame(rsid=seq(1,n_snps),beta=rnorm(n=n_snps,mean=stats$true_beta,sd=stats$se),se=stats$se)
  return(est)
}

## 4. Evaluating MSE for significant SNPs
mse_sig_evaluate <- function(out,true_beta,i=4,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(alpha)/2),]
  if (nrow(snp_sig) == 0){return(0)}
  mse_sig <- mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2)
  return(mse_sig)
}

## 5. Evaluating improvement of average MSE for significant SNPs due to method
## implementation
mse_sig_improve <- function(out,true_beta,i=4,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(alpha)/2),]
  if (nrow(snp_sig) == 0){return(0)}
  mse_sig_improve <- mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2)-mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2)
  return(mse_sig_improve)
}

## 6. Evaluating relative improvement of average MSE for significant SNPs due
## to method implementation
mse_sig_improve_per <- function(out,true_beta,i=4,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(alpha)/2),]
  if (nrow(snp_sig) == 0){return(0)}
  mse_sig_improve <- (mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2)-mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2))/(mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2))
  return(mse_sig_improve)
}


## 7. Evaluating fraction of significant SNPs that are less biased due to method
## implementation
frac_sig_less_bias <- function(out,true_beta,i=4,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(alpha)/2),]
  if (nrow(snp_sig) == 0){return(0)}
  flb <- sum(abs(true_beta[snp_sig$rsid] - snp_sig$beta)>abs(true_beta[snp_sig$rsid] - snp_sig[,i]))/length(snp_sig$rsid)
  return(flb)
}


## 8. Obtains mean of bias evaluation measure over only those simulations in
## which at least one significant SNP has been found
average <- function(ave_vector){
  ave_vector <- ave_vector[which(ave_vector !=0)]
  ave <- mean(ave_vector)
  return(ave)
}


## 9. Obtains standard deviation of bias evaluation measure over only those
## simulations in which at least one significant SNP has been found
std <- function(ave_vector){
  ave_vector <- ave_vector[which(ave_vector !=0)]
  error <- sd(ave_vector)
  return(error)
}


## 10. Takes the 'long' results data frame which contains bias evaluation measure
## for each simulation and outputs an average and standard deviation of this
## bias evaluation measure over simulations of the same parameter values, in
## which at least one significant SNP has been found
ave_results <- function(res_vec, n_sim){
  ave_mse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_mse <- c(rep(0,nrow(res_vec)/n_sim))
  for (i in 1:(nrow(res_vec)/n_sim)){
    ave_mse[i] <- average(res_vec$results[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_mse[i] <- std(res_vec$results[(i*n_sim-(n_sim-1)):(i*n_sim)])
  }
  keep <- seq(from=1,to=nrow(res_vec)-n_sim+1,by=n_sim)
  res_vec_ave <- res_vec[rownames(res_vec) %in% keep,]
  res_vec_ave$results <- ave_mse
  res_vec_ave$res_error <- sd_mse
  res_vec_ave <- subset(res_vec_ave, select = -sim )
  return(res_vec_ave)
}


## 11. Takes the 'long' results data frame which contains bias evaluation
## measure for each simulation and outputs an average and standard deviation of
## this bias evaluation measure over all simulations of the same parameter
## values - currently only used in finding average number of significant SNPs in
## nsig_prop_bias_100sim.R
ave_results1 <- ave_results <- function(res_vec, n_sim){
  ave_mse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_mse <- c(rep(0,nrow(res_vec)/n_sim))
  for (i in 1:(nrow(res_vec)/n_sim)){
    ave_mse[i] <- mean(res_vec$results[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_mse[i] <- sd(res_vec$results[(i*n_sim-(n_sim-1)):(i*n_sim)])
  }
  keep <- seq(from=1,to=nrow(res_vec)-n_sim+1,by=n_sim)
  res_vec_ave <- res_vec[rownames(res_vec) %in% keep,]
  res_vec_ave$results <- ave_mse
  res_vec_ave$res_error <- sd_mse
  res_vec_ave <- subset(res_vec_ave, select = -sim )
  return(res_vec_ave)
}
