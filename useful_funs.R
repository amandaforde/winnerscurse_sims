## FUNCTIONS REQUIRED FOR SIMULATIONS:

## 1. Given values for heritability (H), polygenicity (Pi), sample size (nid)
## and selection coefficient (S), 'true' values of effect size and corresponding
## standard error are simulated in which effect sizes are assumed to be normally
## distributed
simulate_ss <- function(H,Pi,nid,sc,rep_nid=1){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- rnorm(effect_snps,0,sd=sqrt((2*maf*(1-maf))^sc))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/H
  true_beta <- true_beta/sqrt(var_y)
  true_beta <- c(true_beta, rep(0,n_snps-effect_snps))
  se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(nid-2)*maf*(1-maf)))
  rep_se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(rep_nid*nid-2)*maf*(1-maf)))
  stats <- data.frame(true_beta,se,rep_se)
  return(stats)
}


## 2. Given values for heritability (H), polygenicity (Pi), sample size (nid)
## and selection coefficient (S), 'true' values of effect size and corresponding
## standard error are simulated in which effect sizes are assumed to have a
## bimodal distribution, i.e. 50% of effect sizes come from a normal distribution
## centered at 0 while the other half are generated from a normal distribution
## with mean 2.5
simulate_ss_bim <- function(H,Pi,nid,sc,rep_nid=1){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- c(rnorm(0.5*effect_snps,2.5,sd=sqrt((2*maf*(1-maf))^sc)),rnorm(0.5*effect_snps,0,sd=sqrt((2*maf[(0.5*effect_snps+1):effect_snps]*(1-maf[(0.5*effect_snps+1):effect_snps]))^sc)))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/H
  true_beta <- true_beta/sqrt(var_y)
  true_beta <- c(true_beta, rep(0,n_snps-effect_snps))
  se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(nid-2)*maf*(1-maf)))
  rep_se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(rep_nid*nid-2)*maf*(1-maf)))
  stats <- data.frame(true_beta,se,rep_se)
  return(stats)
}

## 3. Simulates discovery GWAS summary statistics given an output from
## functions 1 or 2
simulate_est <- function(stats){
  est <- data.frame(rsid=seq(1,n_snps),beta=rnorm(n=n_snps,mean=stats$true_beta,sd=stats$se),se=stats$se)
  return(est)
}

## 4. Evaluating MSE for significant SNPs
mse_sig_evaluate <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(alpha)/2),]
  if (nrow(snp_sig) == 0){return(100)}
  mse_sig <- mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2)
  return(mse_sig)
}

## 5. Evaluating improvement of average MSE for significant SNPs due to method
## implementation
mse_sig_improve <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(alpha)/2),]
  if (nrow(snp_sig) == 0){return(100)}
  mse_sig_improve <- mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2)-mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2)
  return(mse_sig_improve)
}

## 5. Evaluating improvement of average RMSE for significant SNPs due to method
## implementation
mse_sig_improve_root <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(alpha)/2),]
  if (nrow(snp_sig) == 0){return(100)}
  mse_sig_improve <- sqrt(mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2))-sqrt(mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2))
  return(mse_sig_improve)
}

## 6. Evaluating relative improvement of average MSE for significant SNPs due
## to method implementation
mse_sig_improve_per <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(alpha)/2),]
  if (nrow(snp_sig) == 0){return(100)}
  mse_sig_improve <- (mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2)-mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2))/(mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2))
  return(mse_sig_improve)
}


## 7. Evaluating fraction of significant SNPs that are less biased due to method
## implementation
frac_sig_less_bias <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm(1-(alpha)/2),]
  if (nrow(snp_sig) == 0){return(100)}
  flb <- sum(abs(true_beta[snp_sig$rsid] - snp_sig$beta)>abs(true_beta[snp_sig$rsid] - snp_sig[,i]))/length(snp_sig$rsid)
  return(flb)
}


## 8. Obtains mean of bias evaluation measure over only those simulations in
## which at least one significant SNP has been found
average <- function(ave_vector){
  ave_vector <- ave_vector[which(ave_vector != 100)]
  ave <- mean(ave_vector)
  return(ave)
}


## 9. Obtains standard deviation of bias evaluation measure over only those
## simulations in which at least one significant SNP has been found
std <- function(ave_vector){
  ave_vector <- ave_vector[which(ave_vector != 100)]
  error <- sd(ave_vector)
  return(error)
}


## 10. Takes the 'long' results data frame which contains bias evaluation
## measures for each simulation and outputs averages and standard deviations of
## these bias evaluation measures over simulations of the same parameter values,
## in which at least one significant SNP has been found
ave_results <- function(res_vec, n_sim){
  ave_flb <- c(rep(0,nrow(res_vec)/n_sim))
  sd_flb <- c(rep(0,nrow(res_vec)/n_sim))
  ave_mse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_mse <- c(rep(0,nrow(res_vec)/n_sim))
  ave_rmse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_rmse <- c(rep(0,nrow(res_vec)/n_sim))
  ave_rel_mse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_rel_mse <- c(rep(0,nrow(res_vec)/n_sim))
  for (i in 1:(nrow(res_vec)/n_sim)){
    ave_flb[i] <- average(res_vec$flb[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_flb[i] <- std(res_vec$flb[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_mse[i] <- average(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_mse[i] <- std(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_rmse[i] <- average(res_vec$rmse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_rmse[i] <- std(res_vec$rmse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_rel_mse[i] <- average(res_vec$rel_mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_rel_mse[i] <- std(res_vec$rel_mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
  }
  keep <- seq(from=1,to=nrow(res_vec)-n_sim+1,by=n_sim)
  res_vec_ave <- res_vec[rownames(res_vec) %in% keep,]
  res_vec_ave$flb <- ave_flb
  res_vec_ave$flb_error <- sd_flb
  res_vec_ave$mse <- ave_mse
  res_vec_ave$mse_error <- sd_mse
  res_vec_ave$rmse <- ave_rmse
  res_vec_ave$rmse_error <- sd_rmse
  res_vec_ave$rel_mse <- ave_rel_mse
  res_vec_ave$rel_mse_error <- sd_rel_mse
  res_vec_ave <- subset(res_vec_ave, select = -sim )
  return(res_vec_ave)
}


## 11. Takes the 'long' results data frame which contains number of significant
## SNPs and proportion of these SNP that are 'biased' for each simulation and
## outputs averages and standard deviations for both quantities over all
## simulations of the same parameter values
## NOTE: statistics for proportion of biased significant SNPs are only obtained
## over simulations in which at least one significant SNP was detected
ave_results1 <- function(res_vec, n_sim){
  ave_nsig <- c(rep(0,nrow(res_vec)/n_sim))
  sd_nsig <- c(rep(0,nrow(res_vec)/n_sim))
  ave_pb <- c(rep(0,nrow(res_vec)/n_sim))
  sd_pb <- c(rep(0,nrow(res_vec)/n_sim))
  ave_px <- c(rep(0,nrow(res_vec)/n_sim))
  sd_px <- c(rep(0,nrow(res_vec)/n_sim))
  ave_mse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_mse <- c(rep(0,nrow(res_vec)/n_sim))
  for (i in 1:(nrow(res_vec)/n_sim)){
    ave_nsig[i] <- mean(res_vec$n_sig[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_nsig[i] <- sd(res_vec$n_sig[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_pb[i] <- mean(res_vec$prop_bias[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$prop_bias[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
    sd_pb[i] <- sd(res_vec$prop_bias[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$prop_bias[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
    ave_px[i] <- mean(res_vec$prop_x[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$prop_x[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
    sd_px[i] <- sd(res_vec$prop_x[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$prop_x[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
    ave_mse[i] <- mean(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
    sd_mse[i] <- sd(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
  }
  keep <- seq(from=1,to=nrow(res_vec)-n_sim+1,by=n_sim)
  res_vec_ave <- res_vec[rownames(res_vec) %in% keep,]
  res_vec_ave$n_sig <- ave_nsig
  res_vec_ave$n_sig_sd <- sd_nsig
  res_vec_ave$prop_bias <- ave_pb
  res_vec_ave$prop_bias_sd <- sd_pb
  res_vec_ave$prop_x <- ave_px
  res_vec_ave$prop_x_sd <- sd_px
  res_vec_ave$mse <- ave_mse
  res_vec_ave$mse_sd <- sd_mse
  res_vec_ave <- subset(res_vec_ave, select = -sim )
  return(res_vec_ave)
}


## 12. Given values for heritability (H), polygenicity (Pi), sample size (nid)
## and selection coefficient (S), 'true' values of effect size and corresponding
## standard error are simulated in which effect sizes are assumed to have a
## skewed distribution
simulate_ss_exp <- function(H,Pi,nid,sc,rep_nid=1){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- c(-rexp(n=.1*effect_snps,rate = 1/sqrt((2*maf*(1-maf))^sc)), rexp(n=.9*effect_snps,rate = 1/sqrt((2*maf[(0.1*effect_snps+1):effect_snps]*(1-maf[(0.1*effect_snps+1):effect_snps]))^sc)))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/H
  true_beta <- true_beta/sqrt(var_y)
  true_beta <- c(true_beta, rep(0,n_snps-effect_snps))
  se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(nid-2)*maf*(1-maf)))
  rep_se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(rep_nid*nid-2)*maf*(1-maf)))
  stats <- data.frame(true_beta,se,rep_se)
  return(stats)
}


## 13. Empirical Bayes using scam function, very similar to that in the
## winnerscurse package
empirical_bayes_scam <- function(summary_data){
  z <- summary_data$beta/summary_data$se
  bins <- seq(min(z),max(z),length.out=120)
  mids <- (bins[-length(bins)]+bins[-1])/2
  counts <- graphics::hist(z,breaks=bins,plot=F)$counts
  data <- data.frame(counts,mids)

  f1 <- scam(counts[mids >= 0] ~ s(mids[mids >= 0], bs="mpd"), family =poisson(link="log"), data = data)$fit
  f2 <- scam(counts[mids < 0] ~ s(mids[mids < 0], bs="mpi"), family = poisson(link="log"), data = data)$fit

  f <- c(f2,f1)
  f[f==0] <- min(f[f>0])
  log_f <- as.vector(log(f))
  diff <- diff(log_f)/diff(mids)
  mids2 <- (mids[-length(mids)]+mids[-1])/2

  diff_interpol <- stats::approx(mids2,diff,mids,rule=2,ties=mean)$y
  mids_est <- c(rep(0,length(mids)))
  mids_est[mids>0] <- pmax(0, mids[mids>0] + diff_interpol[mids>0])
  mids_est[mids<0] <- pmin(0, mids[mids<0] + diff_interpol[mids<0])
  z_hat <- stats::approx(mids,mids_est,z,rule=1,ties=mean)$y
  z_hat[is.na(z_hat) & z > 0] <-  z[is.na(z_hat) & z > 0] +   diff_interpol[length(diff_interpol)]
  z_hat[is.na(z_hat) & z < 0] <-  z[is.na(z_hat) & z < 0] + diff_interpol[1]
  z_hat <- sign(z)*pmin(abs(z),abs(z_hat))
  beta_EB <- z_hat*summary_data$se
  summary_data <- cbind(summary_data,beta_EB)
  summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))
  return(summary_data)
}


## 14. Empirical Bayes using gam function with poisson
empirical_bayes_gam_po <- function(summary_data){
  z <- summary_data$beta/summary_data$se
  bins <- seq(min(z),max(z),length.out=120)
  mids <- (bins[-length(bins)]+bins[-1])/2
  counts <- graphics::hist(z,breaks=bins,plot=F)$counts
  data <- data.frame(counts,mids)

  f <- mgcv::gam(counts ~ s(mids), family=poisson(), data=data)$fit
  f[f==0] <- min(f[f>0])
  log_f <- as.vector(log(f))
  diff <- diff(log_f)/diff(mids)
  mids2 <- (mids[-length(mids)]+mids[-1])/2

  diff_interpol <- stats::approx(mids2,diff,mids,rule=2,ties=mean)$y
  mids_est <- c(rep(0,length(mids)))
  mids_est[mids>0] <- pmax(0, mids[mids>0] + diff_interpol[mids>0])
  mids_est[mids<0] <- pmin(0, mids[mids<0] + diff_interpol[mids<0])
  z_hat <- stats::approx(mids,mids_est,z,rule=1,ties=mean)$y
  z_hat[is.na(z_hat) & z > 0] <-  z[is.na(z_hat) & z > 0] +   diff_interpol[length(diff_interpol)]
  z_hat[is.na(z_hat) & z < 0] <-  z[is.na(z_hat) & z < 0] + diff_interpol[1]
  z_hat <- sign(z)*pmin(abs(z),abs(z_hat))
  beta_EB <- z_hat*summary_data$se
  summary_data <- cbind(summary_data,beta_EB)
  summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))
  return(summary_data)
}



## 15. Empirical Bayes using gam function with negative binomial
empirical_bayes_gam_nb <- function(summary_data){
  z <- summary_data$beta/summary_data$se
  bins <- seq(min(z),max(z),length.out=120)
  mids <- (bins[-length(bins)]+bins[-1])/2
  counts <- graphics::hist(z,breaks=bins,plot=F)$counts
  data <- data.frame(counts,mids)

  f <- mgcv::gam(counts ~ s(mids), family=nb(), data=data)$fit
  f[f==0] <- min(f[f>0])
  log_f <- as.vector(log(f))
  diff <- diff(log_f)/diff(mids)
  mids2 <- (mids[-length(mids)]+mids[-1])/2

  diff_interpol <- stats::approx(mids2,diff,mids,rule=2,ties=mean)$y
  mids_est <- c(rep(0,length(mids)))
  mids_est[mids>0] <- pmax(0, mids[mids>0] + diff_interpol[mids>0])
  mids_est[mids<0] <- pmin(0, mids[mids<0] + diff_interpol[mids<0])
  z_hat <- stats::approx(mids,mids_est,z,rule=1,ties=mean)$y
  z_hat[is.na(z_hat) & z > 0] <-  z[is.na(z_hat) & z > 0] +   diff_interpol[length(diff_interpol)]
  z_hat[is.na(z_hat) & z < 0] <-  z[is.na(z_hat) & z < 0] + diff_interpol[1]
  z_hat <- sign(z)*pmin(abs(z),abs(z_hat))
  beta_EB <- z_hat*summary_data$se
  summary_data <- cbind(summary_data,beta_EB)
  summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))
  return(summary_data)
}


## 16. Original conditional likelihood function, doesn't include warning -
## better for simulations possibly also introduced conditions to help avoid
## unusual occurrences in beta.cl2, one function for LD and one for ind
conditional_likelihood_old <- function(summary_data, alpha=5e-8){

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
  stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
  stopifnot(is.numeric(summary_data$beta) && is.numeric(summary_data$se))
  stopifnot(!any(duplicated(summary_data$rsid)))


  z <- summary_data$beta/summary_data$se
  p_val <- 2*(1-stats::pnorm(abs(z)))
  summary_data <- cbind(summary_data, z, p_val)


  if(sum(summary_data$p_val<alpha) == 0){
    summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$z)))
    return(summary_data[,1:4])
  }

  summary_data_sig <- summary_data[summary_data$p_val<alpha,]

  c <- stats::qnorm(1-(alpha)/2)

  beta.cl1 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl2 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl3 <- c(rep(0,nrow(summary_data_sig)))

  for (i in 1:nrow(summary_data_sig)){
    cond.like <- function (mu) { return ((stats::dnorm(summary_data_sig$z[i]-mu))/(stats::pnorm(-c+mu)+stats::pnorm(-c-mu)))}
    mean.cond.like <- function (mu) {return(mu*cond.like(mu))}

    beta.cl1[i] <- (stats::optimize(cond.like, c(0,summary_data_sig$z[i]), maximum=TRUE)$maximum)*summary_data_sig$se[i]

    if(abs(summary_data_sig$z[i]) < 37){
      beta.cl2[i] <- ((stats::integrate(mean.cond.like,-37,37)$value)/(stats::integrate(cond.like,-37,37)$value))*(summary_data_sig$se[i])
    }else{
      if(abs(summary_data_sig$z[i]) > 100){
        beta.cl2[i] <- summary_data_sig$beta[i]
      }else{
        beta.cl2[i] <- ((stats::integrate(mean.cond.like,-100,100)$value)/(stats::integrate(cond.like,-100,100)$value))*(summary_data_sig$se[i])
      }
    }

    beta.cl3[i] <- (beta.cl1[i]+beta.cl2[i])/2
  }

  summary_data_sig <- cbind(summary_data_sig,beta.cl1,beta.cl2,beta.cl3)


  if(sum(abs(beta.cl2) > abs(summary_data_sig$beta + sign(summary_data_sig$beta))) >= 1){
    bad.cl2 <- which((abs(beta.cl2) > abs(summary_data_sig$beta + sign(summary_data_sig$beta))) == TRUE)
    bad.cl2.df <- data.frame(rsid = summary_data_sig$rsid[bad.cl2], z = summary_data_sig$z[bad.cl2], z.cl2 = beta.cl2[bad.cl2]/summary_data_sig$se[bad.cl2], beta = summary_data_sig$beta[bad.cl2], beta.cl2 = beta.cl2[bad.cl2], se = summary_data_sig$se[bad.cl2])
    good.cl2.df <- data.frame(rsid = summary_data_sig$rsid[-(bad.cl2)], z = summary_data_sig$z[-(bad.cl2)], z.cl2 = beta.cl2[-(bad.cl2)]/summary_data_sig$se[-(bad.cl2)], beta = summary_data_sig$beta[-bad.cl2], beta.cl2 = beta.cl2[-(bad.cl2)], se = summary_data_sig$se[-(bad.cl2)])
    for(i in 1:nrow(bad.cl2.df)){
      if(length(good.cl2.df$z[good.cl2.df$z >= bad.cl2.df$z[i]]) == 0 || length(good.cl2.df$z[good.cl2.df$z < bad.cl2.df$z[i]]) == 0){bad.cl2.df$z.cl2[i] <- good.cl2.df$z.cl2[which.min(abs(good.cl2.df$z-bad.cl2.df$z[i]))]
      bad.cl2.df$beta.cl2[i] <- bad.cl2.df$z.cl2[i]*bad.cl2.df$se[i]}else{
        min_greater <- min(good.cl2.df$z[good.cl2.df$z >= bad.cl2.df$z[i]])
        min_smaller <- min(good.cl2.df$z[good.cl2.df$z < bad.cl2.df$z[i]])
        bad.cl2.df$z.cl2[i] <- ((abs(min_greater - bad.cl2.df$z[i]))*(good.cl2.df$z.cl2[good.cl2.df$z == min_smaller]) + (abs(bad.cl2.df$z[i] -     min_smaller))*(good.cl2.df$z.cl2[good.cl2.df$z == min_greater]))/(abs(min_greater - min_smaller))
        bad.cl2.df$beta.cl2[i] <- bad.cl2.df$z.cl2[i]*bad.cl2.df$se[i]
      }

      summary_data_sig$beta.cl2[summary_data_sig$rsid == bad.cl2.df$rsid[i]] <- bad.cl2.df$beta.cl2[i]
      summary_data_sig$beta.cl3[summary_data_sig$rsid == bad.cl2.df$rsid[i]] <- (summary_data_sig$beta.cl1[summary_data_sig$rsid == bad.cl2.df$rsid[i]] +bad.cl2.df$beta.cl2[i])/2
    }
  }


  summary_data_sig <- dplyr::arrange(summary_data_sig,dplyr::desc(abs(summary_data_sig$z)))
  summary_data_sig <- cbind(summary_data_sig[,1:4],summary_data_sig[,7:9])
  return(summary_data_sig)

}



conditional_likelihood_ind <- function(summary_data, alpha=5e-8){

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
  stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
  stopifnot(is.numeric(summary_data$beta) && is.numeric(summary_data$se))
  stopifnot(!any(duplicated(summary_data$rsid)))


  z <- summary_data$beta/summary_data$se
  p_val <- 2*(1-stats::pnorm(abs(z)))
  summary_data <- cbind(summary_data, z, p_val)


  if(sum(summary_data$p_val<alpha) == 0){
    summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$z)))
    return(summary_data[,1:4])
  }

  summary_data_sig <- summary_data[summary_data$p_val<alpha,]

  c <- stats::qnorm(1-(alpha)/2)

  beta.cl1 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl2 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl3 <- c(rep(0,nrow(summary_data_sig)))

  for (i in 1:nrow(summary_data_sig)){
    cond.like <- function (mu) { return ((stats::dnorm(summary_data_sig$z[i]-mu))/(stats::pnorm(-c+mu)+stats::pnorm(-c-mu)))}
    mean.cond.like <- function (mu) {return(mu*cond.like(mu))}

    beta.cl1[i] <- (stats::optimize(cond.like, c(0,summary_data_sig$z[i]), maximum=TRUE)$maximum)*summary_data_sig$se[i]

    if(abs(summary_data_sig$z[i]) < 37){
      beta.cl2[i] <- ((stats::integrate(mean.cond.like,-37,37)$value)/(stats::integrate(cond.like,-37,37)$value))*(summary_data_sig$se[i])
    }else{
      if(abs(summary_data_sig$z[i]) > 100){
        beta.cl2[i] <- summary_data_sig$beta[i]
      }else{
        beta.cl2[i] <- ((stats::integrate(mean.cond.like,-100,100)$value)/(stats::integrate(cond.like,-100,100)$value))*(summary_data_sig$se[i])
      }
    }

    beta.cl3[i] <- (beta.cl1[i]+beta.cl2[i])/2
  }

  summary_data_sig <- cbind(summary_data_sig,beta.cl1,beta.cl2,beta.cl3)


  if(sum(abs(beta.cl2) > abs(summary_data_sig$beta + sign(summary_data_sig$beta))) >= 1){
    bad.cl2 <- which((abs(beta.cl2) > abs(summary_data_sig$beta + sign(summary_data_sig$beta))) == TRUE)
    bad.cl2.df <- data.frame(rsid = summary_data_sig$rsid[bad.cl2], z = summary_data_sig$z[bad.cl2], z.cl2 = beta.cl2[bad.cl2]/summary_data_sig$se[bad.cl2], beta = summary_data_sig$beta[bad.cl2], beta.cl2 = beta.cl2[bad.cl2], se = summary_data_sig$se[bad.cl2])
    good.cl2.df <- data.frame(rsid = summary_data_sig$rsid[-(bad.cl2)], z = summary_data_sig$z[-(bad.cl2)], z.cl2 = beta.cl2[-(bad.cl2)]/summary_data_sig$se[-(bad.cl2)], beta = summary_data_sig$beta[-bad.cl2], beta.cl2 = beta.cl2[-(bad.cl2)], se = summary_data_sig$se[-(bad.cl2)])
    for(i in 1:nrow(bad.cl2.df)){
      if(length(good.cl2.df$z[good.cl2.df$z >= bad.cl2.df$z[i]]) == 0 || length(good.cl2.df$z[good.cl2.df$z < bad.cl2.df$z[i]]) == 0){bad.cl2.df$z.cl2[i] <- good.cl2.df$z.cl2[which.min(abs(good.cl2.df$z-bad.cl2.df$z[i]))]
      bad.cl2.df$beta.cl2[i] <- bad.cl2.df$z.cl2[i]*bad.cl2.df$se[i]}else{
        min_greater <- min(good.cl2.df$z[good.cl2.df$z >= bad.cl2.df$z[i]])
        min_smaller <- min(good.cl2.df$z[good.cl2.df$z < bad.cl2.df$z[i]])
        bad.cl2.df$z.cl2[i] <- ((abs(min_greater - bad.cl2.df$z[i]))*(good.cl2.df$z.cl2[good.cl2.df$z == min_smaller]) + (abs(bad.cl2.df$z[i] -     min_smaller))*(good.cl2.df$z.cl2[good.cl2.df$z == min_greater]))/(abs(min_greater - min_smaller))
        bad.cl2.df$beta.cl2[i] <- bad.cl2.df$z.cl2[i]*bad.cl2.df$se[i]
      }

      summary_data_sig$beta.cl2[summary_data_sig$rsid == bad.cl2.df$rsid[i]] <- bad.cl2.df$beta.cl2[i]
      summary_data_sig$beta.cl3[summary_data_sig$rsid == bad.cl2.df$rsid[i]] <- (summary_data_sig$beta.cl1[summary_data_sig$rsid == bad.cl2.df$rsid[i]] +bad.cl2.df$beta.cl2[i])/2
    }
  }


  summary_data_sig <- dplyr::arrange(summary_data_sig,dplyr::desc(abs(summary_data_sig$z)))
  summary_data_sig <- cbind(summary_data_sig[,1:3],summary_data_sig[,6:8])
  return(summary_data_sig)

}
