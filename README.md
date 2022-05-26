# Evaluating and comparing Winner's Curse correction methods 

This repository contains R code for simulations to evaluate and compare various Winner's Curse correction methods. These are methods which have been designed to appropriately adjust for the bias induced by Winner's Curse in GWAS SNP-trait association estimates. The adjustment methods are applied using the R package, `winnerscurse`. For more details regarding this package and its functionalities, please see [https://amandaforde.github.io/winnerscurse/](https://amandaforde.github.io/winnerscurse/). 

## Independent SNPs

The simulation study follows a factorial design in which GWAS summary statistics were simulated for a quantitative trait under 24 different genetic architectures, described by combinations of four parameters, namely sample size, heritability, polygenicity (proportion of effect SNPs) and selection coefficient. Assuming independence and a normal distribution of effect sizes, for a fixed array of 1,000,000 SNPs, our strategy entails first obtaining the expected standard error for the true effect size of each SNP and subsequently, sampling an estimated effect size based on this expected standard error. This provides an estimated effect size and corresponding standard error for each SNP. For each of these 24 different genetic architectures, 100 sets of summary statistics are simulated.  

The Winner’s Curse correction methods are applied to each data set, producing adjusted estimated effect sizes for each SNP. The performance of these methods are investigated at two different significance thresholds, namely `5e-8` and `5e-4`, with a stronger focus given to the more commonly used genome-wide significance threshold of `5e-8`.  In order to assess each method’s ability at producing less biased SNP-trait association estimates, the change in mean squared error (MSE) of significant SNPs due to method implementation is computed for each data set and method. Details of the procedure followed may be found `winnerscurse_sims_ind.R` with additional required functions contained in `useful_funs.R`.
 
The described simulation process is repeated in a similar fashion in order to produce GWAS summary statistics for a binary trait with a normal distribution of effect sizes and to evaluate methods in this setting. Furthermore, a quantitative phenotype with a bimodal effect size distribution as well as one with a skewed distribution were also considered. See `bim_skew_bin_10sim.R` for description of these alternative settings. 

## Non-independent SNPs

In addition to this set of simulations in which an independence assumption exists, an attempt was made to simulate more realistic summary statistics for a quantitative trait with a normal effect size distribution. In this case, a simple correlation structure was imposed on the SNPs in order to imitate the presence of linkage disequilibrium (LD) in real data. It was assumed that the same correlation structure exists in independent blocks of 100 SNPs. Details of this simulation process are contained in `winnerscurse_sims_LD.R`. 

**Note:** A rigorous description of all current progress to date, together with illustrations of results, can be viewed at [https://amandaforde.github.io/winnerscurse_sims/](https://amandaforde.github.io/winnerscurse_sims/). So far, the investigation of only those methods that use discovery GWAS summary statistics has been prioritised, rather than methods which combine information gained from both discovery and replication GWASs. 

