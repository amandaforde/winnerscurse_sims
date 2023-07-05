## WINNER'S CURSE SIMULATION STUDY PLOTS:

## This script produces plots which demonstrate the results obtained from
## running pipeline.R. The first section includes plots which relate to the
## initial exploration of the various simulation scenarios, in which a simple
## correlation structure has been imposed on the SNPs. This is followed by the
## results of the application of the Winner's Curse correction methods to these
## simulated data sets. The second section illustrates the results obtained upon
## implementation of the methods using independent sets of SNPs.


library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(expm)
library("RColorBrewer")
col <- brewer.pal(8,"Dark2")
col1 <- brewer.pal(11,"RdYlBu")
col2 <- brewer.pal(11,"PRGn")


################################################################################
################################################################################

## PART 1A) LD SIMULATIONS INITIAL EXPLORATION

## S1 Fig: No. sig. SNPs vs prop. significantly overestimated

nsig_5e_8 <- read.csv("results/nsig_prop_bias_5e_8_LD_all.csv")
nsig_5e_8$Scenario <- c(rep("1",100), rep("2",100), rep("3",100), rep("4",100), rep("5",100), rep("6",100), rep("7",100), rep("8",100))

plot <- ggplot(nsig_5e_8,aes(x=n_sig,y=prop_x,colour=Scenario)) + geom_point(aes(color=Scenario), size=2.2) +
  scale_color_manual(values=c(col[1],col[2],col[3],col[4],col[5],col[6],col[7],col[8])) +
  xlab(expression(paste("No. sig SNPs (", 5%*%10^-8, ")")))+ ylab(expression(paste("Prop. sig SNPs", italic(" significantly "), "overestimated"))) +
  theme_bw() + theme(text = element_text(size=12), legend.position = "bottom", legend.box.background = element_rect(colour = "black"),legend.spacing.y = unit(0, "mm"), legend.background=element_blank()) +
  guides(color=guide_legend(nrow=1, byrow=TRUE))

scenarios <- data.frame(Scenario = c("Scenario",1:8), n = c("sample size",rep(c("30,000","300,000"),4)), h2 = c("heritability",rep(c(0.3, 0.3, 0.8, 0.8),2)), pi = c("polygenicity",rep(0.01,4),rep(0.001,4)))
scenarios <- t(scenarios)
rownames(scenarios) <- NULL
tab <- ggtexttable(scenarios, theme = ttheme("classic", base_size=10))
tab <- table_cell_font(tab, row = 1, column= 1:9, face = "bold", size=10)
tab <- table_cell_font(tab, column = 1, row = 1:4, face = "bold", size=10)
tab <- tab %>%
  tab_add_hline(at.row = c( 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = c(2), column.side = "left",linewidth=3, linetype = 1)

ggarrange(plot, tab, ncol=1, nrow=2, heights=c(1, 0.25))


################################################################################

## S2 Fig: z vs bias plots

## Function to simulate sets of summary statistics with LD

set.seed(1998)

n_snps <- 10^6
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
R_sqrt <- sqrtm(R)

## Function to simulate LD:
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

  ## R defined outside of function as it is fixed!
  for (i in 1:10000){
    D_1 <- diag(x = 1/(sqrt(nid*2*maf_re[((i-1)*100+1):(i*100)]*(1-maf_re[((i-1)*100+1):(i*100)]))), nrow=100, ncol=100)
    D_2 <- diag(x = (sqrt(nid*2*maf_re[((i-1)*100+1):(i*100)]*(1-maf_re[((i-1)*100+1):(i*100)]))), nrow=100, ncol=100)
    A <- D_1 %*% cormat %*% D_2
    sigma <- D_1 %*% cormat %*% D_1
    B <- cormatsq %*% D_1
    mu_all[((i-1)*100+1):(i*100)] <- A %*% true_beta_all[((i-1)*100+1):(i*100)]
    beta_hat_all[((i-1)*100+1):(i*100)] <- mu_all[((i-1)*100+1):(i*100)] + ( B %*% rnorm(100))
    se_all[((i-1)*100+1):(i*100)] <- diag(D_1)
  }

  summary_stats <- data.frame(rsid=seq(1,n_snps),beta=beta_hat_all,se=se_all,true_beta=mu_all)
  return(summary_stats)
}


## Scenario 1
out <- simulate_ss_ld(H=0.3,Pi=0.01,nid=30000,sc=0)
out$z <- out$beta/out$se
out$bias <- out$beta - out$true_beta
subout <- out[(abs(out$beta) > (abs(out$true_beta) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm((5e-4)/2, lower.tail=FALSE)),]
Scenario1 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                              aes(x=z,y=bias),
                                                              color='grey40',
                                                              size=1) + xlab("z") +
  ylab("bias") + ggtitle("Scenario 1", subtitle=expression(paste(italic(n), " = 30,000, " , italic(h)^2, " = 0.3, ", pi, " = 0.01")) ) +  theme_classic() + geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) +
  geom_hline(yintercept=0) +
  theme(text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, face="bold", colour=col[1]),
    plot.subtitle = element_text(hjust = 0.5)
  ) + ylim(-0.1, 0.1)


## Scenario 2
out <- simulate_ss_ld(H=0.3,Pi=0.01,nid=300000,sc=0)
out$z <- out$beta/out$se
out$bias <- out$beta - out$true_beta
subout <- out[(abs(out$beta) > (abs(out$true_beta) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm((5e-4)/2, lower.tail=FALSE)),]
Scenario2 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                              aes(x=z,y=bias),
                                                              color='grey40',
                                                              size=1) + xlab("z") +
  ylab("bias") + ggtitle("Scenario 2", subtitle=expression(paste(italic(n), " = 300,000, " , italic(h)^2, " = 0.3, ", pi, " = 0.01"))) +  theme_classic() + geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) +
  geom_hline(yintercept=0) +
  theme(text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, face="bold", colour=col[2]),
    plot.subtitle = element_text(hjust = 0.5, face="italic")
  ) + ylim(-0.1, 0.1)



## Scenario 3
out <- simulate_ss_ld(H=0.8,Pi=0.01,nid=30000,sc=0)
out$z <- out$beta/out$se
out$bias <- out$beta - out$true_beta
subout <- out[(abs(out$beta) > (abs(out$true_beta) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm((5e-4)/2, lower.tail=FALSE)),]
Scenario3 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                              aes(x=z,y=bias),
                                                              color='grey40',
                                                              size=1) + xlab("z") +
  ylab("bias") + ggtitle("Scenario 3", subtitle=expression(paste(italic(n), " = 30,000, " , italic(h)^2, " = 0.8, ", pi, " = 0.01"))) +  theme_classic() + geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) +
  geom_hline(yintercept=0) +
  theme(text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, face="bold", colour=col[3]),
    plot.subtitle = element_text(hjust = 0.5, face="italic")
  ) + ylim(-0.1, 0.1)



## Scenario 4
out <- simulate_ss_ld(H=0.8,Pi=0.01,nid=300000,sc=0)
out$z <- out$beta/out$se
out$bias <- out$beta - out$true_beta
subout <- out[(abs(out$beta) > (abs(out$true_beta) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm((5e-4)/2, lower.tail=FALSE)),]
Scenario4 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                              aes(x=z,y=bias),
                                                              color='grey40',
                                                              size=1) + xlab("z") +
  ylab("bias") + ggtitle("Scenario 4", subtitle=expression(paste(italic(n), " = 300,000, " , italic(h)^2, " = 0.8, ", pi, " = 0.01"))) +  theme_classic() + geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) +
  geom_hline(yintercept=0) +
  theme(text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, face="bold", colour=col[4]),
    plot.subtitle = element_text(hjust = 0.5, face="italic")
  ) + ylim(-0.1, 0.1)


## Scenario 5
out <- simulate_ss_ld(H=0.3,Pi=0.001,nid=30000,sc=0)
out$z <- out$beta/out$se
out$bias <- out$beta - out$true_beta
subout <- out[(abs(out$beta) > (abs(out$true_beta) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm((5e-4)/2, lower.tail=FALSE)),]
Scenario5 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                              aes(x=z,y=bias),
                                                              color='grey40',
                                                              size=1) + xlab("z") +
  ylab("bias") + ggtitle("Scenario 5", subtitle=expression(paste(italic(n), " = 30,000, " , italic(h)^2, " = 0.3, ", pi, " = 0.001"))) +  theme_classic() + geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) +
  geom_hline(yintercept=0) +
  theme(text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, face="bold", colour=col[5]),
    plot.subtitle = element_text(hjust = 0.5, face="italic")
  ) + ylim(-0.1, 0.1)


## Scenario 6
out <- simulate_ss_ld(H=0.3,Pi=0.001,nid=300000,sc=0)
out$z <- out$beta/out$se
out$bias <- out$beta - out$true_beta
subout <- out[(abs(out$beta) > (abs(out$true_beta) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm((5e-4)/2, lower.tail=FALSE)),]
Scenario6 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                              aes(x=z,y=bias),
                                                              color='grey40',
                                                              size=1) + xlab("z") +
  ylab("bias") + ggtitle("Scenario 6", subtitle=expression(paste(italic(n), " = 300,000, " , italic(h)^2, " = 0.3, ", pi, " = 0.001"))) +  theme_classic() + geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) +
  geom_hline(yintercept=0) +
  theme(text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, face="bold", colour=col[6]),
    plot.subtitle = element_text(hjust = 0.5, face="italic")
  ) + ylim(-0.1, 0.1)


## Scenario 7
out <- simulate_ss_ld(H=0.8,Pi=0.001,nid=30000,sc=0)
out$z <- out$beta/out$se
out$bias <- out$beta - out$true_beta
subout <- out[(abs(out$beta) > (abs(out$true_beta) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm((5e-4)/2, lower.tail=FALSE)),]
Scenario7 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                              aes(x=z,y=bias),
                                                              color='grey40',
                                                              size=1) + xlab("z") +
  ylab("bias") + ggtitle("Scenario 7", subtitle=expression(paste(italic(n), " = 30,000, " , italic(h)^2, " = 0.8, ", pi, " = 0.001"))) +  theme_classic() + geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) +
  geom_hline(yintercept=0) +
  theme(text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, face="bold", colour=col[7]),
    plot.subtitle = element_text(hjust = 0.5, face="italic")
  ) + ylim(-0.1, 0.1)

## Scenario 8
out <- simulate_ss_ld(H=0.8,Pi=0.001,nid=300000,sc=0)
out$z <- out$beta/out$se
out$bias <- out$beta - out$true_beta
subout <- out[(abs(out$beta) > (abs(out$true_beta) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm((5e-4)/2, lower.tail=FALSE)),]
Scenario8 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                              aes(x=z,y=bias),
                                                              color='grey40',
                                                              size=1) + xlab("z") +
  ylab("bias") + ggtitle("Scenario 8", subtitle=expression(paste(italic(n), " = 300,000, " , italic(h)^2, " = 0.8, ", pi, " = 0.001"))) +  theme_classic() + geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm((5e-4)/2, lower.tail=FALSE), colour="darkred", linetype="dashed",size=1) +
  geom_hline(yintercept=0) +
  theme(text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, face="bold", colour=col[8]),
    plot.subtitle = element_text(hjust = 0.5, face="italic")
  ) + ylim(-0.1, 0.1)


## Combine into one figure
ggarrange(Scenario1, Scenario2, Scenario3, Scenario4, Scenario5, Scenario6, Scenario7, Scenario8, labels = c("A", "B", "C", "D", "E", "F", "G", "H"), ncol = 2, nrow = 4)


################################################################################
################################################################################

## PART 1B) LD SIMULATIONS RESULTS

## Fig 1: change in RMSE of sig. SNPs (5e-8)

MSE_5e_8 <- read.csv("results/norm_5e-8_100sim_LD_all.csv")
MSE_5e_8 <- MSE_5e_8  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3', EB_df = "EB_df", EB_gam_nb = "EB_gam_nb", EB_gam_po = "EB_gam_po", EB_scam = "EB_scam"))
MSE_5e_8$method <- factor(MSE_5e_8$method, levels=c("CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))

MSE_5e_8a <- MSE_5e_8[which(MSE_5e_8$n_samples==30000),]
MSE_5e_8a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_8a)))
MSE_5e_8b <- MSE_5e_8[which(MSE_5e_8$n_samples==300000),]
MSE_5e_8b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_8b)))

plotA <- ggplot(MSE_5e_8a,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_8b,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

## S3 Fig: change in RMSE of sig. SNPs (5e-4)

MSE_5e_4 <- read.csv("results/norm_5e-4_100sim_LD_all.csv")
MSE_5e_4 <- MSE_5e_4  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3'))
MSE_5e_4$method <- factor(MSE_5e_4$method, levels=c("CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))

MSE_5e_4a <- MSE_5e_4[which(MSE_5e_4$n_samples==30000),]
MSE_5e_4a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_4a)))
MSE_5e_4b <- MSE_5e_4[which(MSE_5e_4$n_samples==300000),]
MSE_5e_4b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_4b)))

plotA <- ggplot(MSE_5e_4a,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_4b,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))


################################################################################
################################################################################

## PART 1C) INDEPENDENT SIMULATIONS RESULTS

## i) normal distribution

## S4 Fig: change in RMSE of sig. SNPs (5e-8)

MSE_5e_8 <- read.csv("results/norm_5e-8_100sim_all.csv")
MSE_5e_8 <- MSE_5e_8  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3'))
MSE_5e_8$method <- factor(MSE_5e_8$method, levels=c("CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT", "rep"))

MSE_5e_8 <- MSE_5e_8[which(MSE_5e_8$S==0),]
MSE_5e_8 <- MSE_5e_8[MSE_5e_8$flb != 100,]
MSE_5e_8a <- MSE_5e_8[which(MSE_5e_8$n_samples==30000),]
MSE_5e_8a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_8a)))
MSE_5e_8b <- MSE_5e_8[which(MSE_5e_8$n_samples==300000),]
MSE_5e_8b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_8b)))

plotA <- ggplot(MSE_5e_8a,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_8b,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

## S5 Fig: change in RMSE of sig. SNPs (5e-4)

MSE_5e_4 <- read.csv("results/norm_5e-4_100sim_all.csv")
MSE_5e_4 <- MSE_5e_4  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3'))
MSE_5e_4$method <- factor(MSE_5e_4$method, levels=c("CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT", "rep"))

MSE_5e_4 <- MSE_5e_4[which(MSE_5e_4$S==0),]
MSE_5e_4 <- MSE_5e_4[MSE_5e_4$flb != 100,]
MSE_5e_4a <- MSE_5e_4[which(MSE_5e_4$n_samples==30000),]
MSE_5e_4a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_4a)))
MSE_5e_4b <- MSE_5e_4[which(MSE_5e_4$n_samples==300000),]
MSE_5e_4b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_4b)))

plotA <- ggplot(MSE_5e_4a,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_4b,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

################################################################################

## ii) bimodal distribution

## S6 Fig: change in RMSE of sig. SNPs (5e-8)

MSE_5e_8 <- read.csv("results/bim_5e-8_100sim_all.csv")
MSE_5e_8 <- MSE_5e_8[MSE_5e_8$method != "naive",]
MSE_5e_8 <- MSE_5e_8  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3'))
MSE_5e_8$method <- factor(MSE_5e_8$method, levels=c("CL1", "CL2", "CL3", "EB", "boot", "FIQT", "rep"))

MSE_5e_8 <- MSE_5e_8[which(MSE_5e_8$S==0),]
MSE_5e_8 <- MSE_5e_8[MSE_5e_8$flb != 100,]
MSE_5e_8a <- MSE_5e_8[which(MSE_5e_8$n_samples==30000),]
MSE_5e_8a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_8a)))
MSE_5e_8b <- MSE_5e_8[which(MSE_5e_8$n_samples==300000),]
MSE_5e_8b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_8b)))

plotA <- ggplot(MSE_5e_8a,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_8b,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

## S7 Fig: change in RMSE of sig. SNPs (5e-4)

MSE_5e_4 <- read.csv("results/bim_5e-4_100sim_all.csv")
MSE_5e_4 <- MSE_5e_4[MSE_5e_4$method != "naive",]
MSE_5e_4 <- MSE_5e_4  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3'))
MSE_5e_4$method <- factor(MSE_5e_4$method, levels=c("CL1", "CL2", "CL3", "EB", "boot", "FIQT", "rep"))

MSE_5e_4 <- MSE_5e_4[which(MSE_5e_4$S==0),]
MSE_5e_4 <- MSE_5e_4[MSE_5e_4$flb != 100,]
MSE_5e_4a <- MSE_5e_4[which(MSE_5e_4$n_samples==30000),]
MSE_5e_4a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_4a)))
MSE_5e_4b <- MSE_5e_4[which(MSE_5e_4$n_samples==300000),]
MSE_5e_4b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_4b)))

plotA <- ggplot(MSE_5e_4a,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_4b,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

################################################################################

## iii) skewed distribution

## S8 Fig: change in RMSE of sig. SNPs (5e-8)

MSE_5e_8 <- read.csv("results/skew_5e-8_100sim_all.csv")
MSE_5e_8 <- MSE_5e_8[MSE_5e_8$method != "naive",]
MSE_5e_8 <- MSE_5e_8  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3'))
MSE_5e_8$method <- factor(MSE_5e_8$method, levels=c("CL1", "CL2", "CL3", "EB", "boot", "FIQT", "rep"))

MSE_5e_8 <- MSE_5e_8[which(MSE_5e_8$S==0),]
MSE_5e_8 <- MSE_5e_8[MSE_5e_8$flb != 100,]
MSE_5e_8a <- MSE_5e_8[which(MSE_5e_8$n_samples==30000),]
MSE_5e_8a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_8a)))
MSE_5e_8b <- MSE_5e_8[which(MSE_5e_8$n_samples==300000),]
MSE_5e_8b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_8b)))

plotA <- ggplot(MSE_5e_8a,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_8b,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

## S9 Fig: change in RMSE of sig. SNPs (5e-4)

MSE_5e_4 <- read.csv("results/skew_5e-4_100sim_all.csv")
MSE_5e_4 <- MSE_5e_4[MSE_5e_4$method != "naive",]
MSE_5e_4 <- MSE_5e_4  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3'))
MSE_5e_4$method <- factor(MSE_5e_4$method, levels=c("CL1", "CL2", "CL3", "EB", "boot", "FIQT", "rep"))

MSE_5e_4 <- MSE_5e_4[which(MSE_5e_4$S==0),]
MSE_5e_4 <- MSE_5e_4[MSE_5e_4$flb != 100,]
MSE_5e_4a <- MSE_5e_4[which(MSE_5e_4$n_samples==30000),]
MSE_5e_4a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_4a)))
MSE_5e_4b <- MSE_5e_4[which(MSE_5e_4$n_samples==300000),]
MSE_5e_4b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_4b)))

plotA <- ggplot(MSE_5e_4a,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_4b,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

################################################################################

## iv) binary trait

## S10 Fig: change in RMSE of sig. SNPs (5e-8)

MSE_5e_8 <- read.csv("results/bin_5e-8_100sim_all.csv")
MSE_5e_8 <- MSE_5e_8[MSE_5e_8$method != "naive",]
MSE_5e_8 <- MSE_5e_8  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3'))
MSE_5e_8$method <- factor(MSE_5e_8$method, levels=c("CL1", "CL2", "CL3", "EB", "boot", "FIQT", "rep"))

MSE_5e_8 <- MSE_5e_8[which(MSE_5e_8$S==0),]
MSE_5e_8 <- MSE_5e_8[MSE_5e_8$flb != 100,]
MSE_5e_8a <- MSE_5e_8[which(MSE_5e_8$n_samples==30000),]
MSE_5e_8a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_8a)))
MSE_5e_8b <- MSE_5e_8[which(MSE_5e_8$n_samples==300000),]
MSE_5e_8b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_8b)))

plotA <- ggplot(MSE_5e_8a,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_8b,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

## S11 Fig: change in RMSE of sig. SNPs (5e-4)

MSE_5e_4 <- read.csv("results/bin_5e-4_100sim_all.csv")
MSE_5e_4 <- MSE_5e_4[MSE_5e_4$method != "naive",]
MSE_5e_4 <- MSE_5e_4  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3'))
MSE_5e_4$method <- factor(MSE_5e_4$method, levels=c("CL1", "CL2", "CL3", "EB", "boot", "FIQT", "rep"))

MSE_5e_4 <- MSE_5e_4[which(MSE_5e_4$S==0),]
MSE_5e_4 <- MSE_5e_4[MSE_5e_4$flb != 100,]
MSE_5e_4a <- MSE_5e_4[which(MSE_5e_4$n_samples==30000),]
MSE_5e_4a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_4a)))
MSE_5e_4b <- MSE_5e_4[which(MSE_5e_4$n_samples==300000),]
MSE_5e_4b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_4b)))

plotA <- ggplot(MSE_5e_4a,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_4b,aes(x=method,y=rmse,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col1[11], col1[9], col2[11])) + xlab("Method") +
  ylab(expression(paste(italic("Change "), "in RMSE of sig. SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))


################################################################################

## Fig 1: Bias of positive sig. SNPs (5e-8)

MSE_5e_8 <- read.csv("results/norm_5e-8_100sim_LD_all_2.csv")
MSE_5e_8 <- MSE_5e_8  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3', EB_df = "EB_df", EB_gam_nb = "EB_gam_nb", EB_gam_po = "EB_gam_po", EB_scam = "EB_scam", naive = "naive"))
MSE_5e_8$method <- factor(MSE_5e_8$method, levels=c("naive","CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))

MSE_5e_8a <- MSE_5e_8[which(MSE_5e_8$n_samples==30000),]
MSE_5e_8a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_8a)))
MSE_5e_8b <- MSE_5e_8[which(MSE_5e_8$n_samples==300000),]
MSE_5e_8b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_8b)))

plotA <- ggplot(MSE_5e_8a,aes(x=method,y=bias_up,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + xlab("Method") +
  ylab(expression(paste(italic("Average "), "bias of positive SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_8b,aes(x=method,y=bias_up,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + xlab("Method") +
  ylab(expression(paste(italic("Average "), "bias of positive SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))


## Fig 1: Bias of negative sig. SNPs (5e-8)

MSE_5e_8 <- read.csv("results/norm_5e-8_100sim_LD_all_2.csv")
MSE_5e_8 <- MSE_5e_8  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3', EB_df = "EB_df", EB_gam_nb = "EB_gam_nb", EB_gam_po = "EB_gam_po", EB_scam = "EB_scam", naive = "naive"))
MSE_5e_8$method <- factor(MSE_5e_8$method, levels=c("naive","CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))

MSE_5e_8a <- MSE_5e_8[which(MSE_5e_8$n_samples==30000),]
MSE_5e_8a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_8a)))
MSE_5e_8b <- MSE_5e_8[which(MSE_5e_8$n_samples==300000),]
MSE_5e_8b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_8b)))

plotA <- ggplot(MSE_5e_8a,aes(x=method,y=bias_down,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + xlab("Method") +
  ylab(expression(paste(italic("Average "), "bias of negative SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_8b,aes(x=method,y=bias_down,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + xlab("Method") +
  ylab(expression(paste(italic("Average "), "bias of negative SNPs at  ", 5%*%10^-8))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))


## Fig 2: Bias of positive sig. SNPs (5e-4)

MSE_5e_8 <- read.csv("results/norm_5e-4_100sim_LD_all_2.csv")
MSE_5e_8 <- MSE_5e_8  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3', EB_df = "EB_df", EB_gam_nb = "EB_gam_nb", EB_gam_po = "EB_gam_po", EB_scam = "EB_scam", naive = "naive"))
MSE_5e_8$method <- factor(MSE_5e_8$method, levels=c("naive","CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))

MSE_5e_8a <- MSE_5e_8[which(MSE_5e_8$n_samples==30000),]
MSE_5e_8a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_8a)))
MSE_5e_8b <- MSE_5e_8[which(MSE_5e_8$n_samples==300000),]
MSE_5e_8b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_8b)))

plotA <- ggplot(MSE_5e_8a,aes(x=method,y=bias_up,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + xlab("Method") +
  ylab(expression(paste(italic("Average "), "bias of positive SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_8b,aes(x=method,y=bias_up,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + xlab("Method") +
  ylab(expression(paste(italic("Average "), "bias of positive SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))


## Fig 2: Bias of negative sig. SNPs (5e-4)

MSE_5e_8 <- read.csv("results/norm_5e-4_100sim_LD_all_2.csv")
MSE_5e_8 <- MSE_5e_8  %>% mutate(method = recode(method, BR = 'boot', cl1 = 'CL1', cl2 = 'CL2', cl3 = 'CL3', EB_df = "EB_df", EB_gam_nb = "EB_gam_nb", EB_gam_po = "EB_gam_po", EB_scam = "EB_scam", naive = "naive"))
MSE_5e_8$method <- factor(MSE_5e_8$method, levels=c("naive","CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))

MSE_5e_8a <- MSE_5e_8[which(MSE_5e_8$n_samples==30000),]
MSE_5e_8a$n_samples <- c(rep(c("30,000"),nrow(MSE_5e_8a)))
MSE_5e_8b <- MSE_5e_8[which(MSE_5e_8$n_samples==300000),]
MSE_5e_8b$n_samples <- c(rep(c("300,000"),nrow(MSE_5e_8b)))

plotA <- ggplot(MSE_5e_8a,aes(x=method,y=bias_down,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + xlab("Method") +
  ylab(expression(paste(italic("Average "), "bias of negative SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 30,000")
plotB <- ggplot(MSE_5e_8b,aes(x=method,y=bias_down,fill=method, color=method)) + geom_boxplot(size=0.7,aes(fill=method, color=method, alpha=0.2)) + facet_grid(h2~prop_effect,labeller=label_both) + scale_fill_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + scale_color_manual(values=c(col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11], col1[9], col2[4])) + xlab("Method") +
  ylab(expression(paste(italic("Avergae "), "bias of negative SNPs at  ", 5%*%10^-4))) + theme_bw() + geom_hline(yintercept=0, colour="black") +  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),text = element_text(size=12),legend.position = "none", strip.text = element_text(face="italic")) + ggtitle("Sample size: 300,000")

figure <- plotA + plotB
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))





################################################################################
################################################################################



mu <- rep(seq(from=4, to = 8, by = 0.05),10000)
z <- c()

for(i in 1:length(mu)){
  z <- c(z,rnorm(1,mu[i],1))
}

data <- data.frame(rsid = 1:length(mu), beta = z, se = rep(1,length(mu)), mu = mu)

data_sig <- data[abs(data$beta) > qnorm((5e-8)/2, lower.tail=FALSE),]

out_CL <- conditional_likelihood(data[,1:3])

out_CL$mu <- data$mu[out_CL$rsid]

mu1 <- unique(mu)
bias <- c()
bias_1 <-c()
se <- c()
se_1 <- c()
RMSE <- c()
RMSE_1 <- c()

for(i in 1:length(mu1)){
  bias <- c(bias,sum(out_CL$beta.cl1[out_CL$mu == mu1[i]] - mu1[i])/length(out_CL$beta.cl1[out_CL$mu == mu1[i]]))
  bias_1 <- c(bias_1,sum(out_CL$beta[out_CL$mu == mu1[i]] - mu1[i])/length(out_CL$beta[out_CL$mu == mu1[i]]))
  RMSE <- c(RMSE, sqrt(sum((out_CL$beta.cl1[out_CL$mu == mu1[i]] - mu1[i])^2)/length(out_CL$beta.cl1[out_CL$mu == mu1[i]])))
  RMSE_1 <- c(RMSE_1, sqrt(sum((out_CL$beta[out_CL$mu == mu1[i]] - mu1[i])^2)/length(out_CL$beta[out_CL$mu == mu1[i]])))
  se <- c(se, sqrt(sum((out_CL$beta.cl1[out_CL$mu == mu1[i]] - mean(out_CL$beta.cl1[out_CL$mu == mu1[i]]))^2)/length(out_CL$beta.cl1[out_CL$mu == mu1[i]])))
  se_1 <- c(se_1, sqrt(sum((out_CL$beta[out_CL$mu == mu1[i]] - mean(out_CL$beta[out_CL$mu == mu1[i]]))^2)/length(out_CL$beta[out_CL$mu == mu1[i]])))
}


results <- data.frame(mu = mu1, bias = bias, RMSE = RMSE, se = se, bias_1 = bias_1, RMSE_1 = RMSE_1, se_1 = se_1)

bias_plot <- ggplot(results, aes(x=mu, y=bias)) + geom_point(size=1.5) + geom_point(data=results,
                                                                                  aes(x=mu,y=bias_1),
                                                                                  color='blue3',
                                                                                  size=1.5,
                                                                                  shape=15) +
  geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_hline(yintercept=0) + ylab("bias") +
  xlab(expression(paste("True standardized effect size, ", mu )))

se_plot <- ggplot(results, aes(x=mu, y=se)) + geom_point(size=1.5) + geom_point(data=results,
                                                                              aes(x=mu,y=se_1),
                                                                              color='blue3',
                                                                              size=1.5,
                                                                              shape=15) +
  geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_hline(yintercept=1) + ylab("SE") +
  xlab(expression(paste("True standardized effect size, ", mu )))

RMSE_plot <- ggplot(results, aes(x=mu, y=RMSE)) + geom_point(size=1.5) + geom_point(data=results,
                                                                                  aes(x=mu,y=RMSE_1),
                                                                                  color='blue3',
                                                                                  size=1.5,
                                                                                  shape=15) +
  geom_vline(xintercept=qnorm((5e-8)/2, lower.tail=FALSE), colour="red", linetype="dashed",size=1) +
  geom_hline(yintercept=1) + ylab("RMSE") +
  xlab(expression(paste("True standardized effect size, ", mu )))

library(patchwork)
library(ggpubr)

figure <- bias_plot + se_plot + RMSE_plot
figure + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

