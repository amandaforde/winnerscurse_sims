## WINNER'S CURSE SIMULATION STUDY:

## This script allows us to run the entire simulation study which evaluates
## and compares various winner's curse correction methods, from beginning to
## end.


## Load all required packages for simulations.
library(devtools)
devtools::install_github("amandaforde/winnerscurse")
library(winnerscurse)
library(tidyverse)
library(parallel)
library(scam)
library(mgcv)
library(expm)

## Run 'useful_funs.R' in order to define additional functions required for
## the simulations below.
source("useful_funs.R")


## Run each set of simulations, taking note of length of time taken for each
## script.
## NOTE: Simulations currently being run on Windows, hence mc.cores=1 in
## mclapply() in below scripts.


start_time <- Sys.time()
source("nsig_prop_bias_LD_100sim.R")
end_time <- Sys.time()
print("nsig_prop_bias_LD_100sim.R complete!")
end_time - start_time

start_time <- Sys.time()
source("winnerscurse_sims_LD.R")
end_time <- Sys.time()
print("winnerscurse_sims_LD.R complete!")
end_time - start_time

start_time <- Sys.time()
source("nsig_prop_bias_100sim.R")
end_time <- Sys.time()
print("nsig_prop_bias_100sim.R complete!")
end_time - start_time

start_time <- Sys.time()
source("winnerscurse_sims_ind.R")
end_time <- Sys.time()
print("winnerscurse_sims_ind.R complete!")
end_time - start_time

start_time <- Sys.time()
source("bim_skew_bin_100sim.R")
end_time <- Sys.time()
print("bim_skew_bin_100sim.R complete!")
end_time - start_time
