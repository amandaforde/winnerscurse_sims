# Evaluating and comparing Winner's Curse adjustment methods 

This repository contains R code for simulations to evaluate and compare various Winner's Curse correction methods. These are methods which have been designed to appropriately adjust for the bias induced by Winner's Curse in GWAS SNP-trait association estimates. The adjustment methods are applied using the R package, `winnerscurse`. For more details regarding this package and its functionalities, please see [https://amandaforde.github.io/winnerscurse/](https://amandaforde.github.io/winnerscurse/). 

GWAS summary statistics are simulated under 24 individual scenarios based on different combinations of values for heritability, polygenicity, sample size and selection coefficient. Winner's Curse correction methods are applied to each of these scenarios and their performance are evaluated under three bias evaluation metrics. 

Details of all current progress to date, together with illustrations of results, can be viewed at [https://amandaforde.github.io/winnerscurse_sims/](https://amandaforde.github.io/winnerscurse_sims/).
