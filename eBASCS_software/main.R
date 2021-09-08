# eBASCS: Disentangling Overlapping Sources II, using Spatial, Spectral and Temporal Information.
# Code for Meyer et al. (2021) paper
# Antoine D. Meyer, Imperial College London (+ additional pieces of code from Luis Campos and David Jones)

# Load Required R packages

library("MASS")
library("MCMCpack")
library("cubature")
library("Rcpp")

# Initialize time
start <- proc.time()
initial_run = NULL
home = getwd()

# Choose and load dataset
file.name <- 'uv_data6.Rda'
file.name = gsub('./data/', '', file.name)
data_dir = '/data/'
load(paste(home, data_dir, file.name, sep = ''))

# Format output file
results <- paste(home,"/results/",sep="")
out.file.name <- 'myresults.Rda'

# Choose model to fit. Options are: "spatial", "BASCS", "spacetime", and "eBASCS".
model <- 'spacetime'

# Clean the dataset to prepare for MCMC input
source('_clean_data.R')

# Choose model options
source('_model_options.R')

# Choose prior distributions
source('_prior_options.R')

#Choose MCMC settings
source('_mcmc_options.R')

# Initialization of sampler
source('_init_time.R')

# Source required code
setwd(paste(home, "/additional_functions/", sep=""))
source(paste("extended_full_log_posterior",function_load_name,sep=""))
source(paste("full_log_posterior",function_load_name,sep=""))
source("logit.R")
source("inv.logit.R")
source(paste("mc.logext",function_load_name,sep=""))
source(paste("mc.spacetime",function_load_name,sep=""))
source(paste("mc.spatial",function_load_name,sep=""))
source(paste("mc.bext",function_load_name,sep=""))

if (model == 'eBASCS') {
  fixedk_eBASCS_mcmc <- mc.logext
}

if (model == 'spacetime') {
  fixedk_spacetime_mcmc <- mc.spacetime
}

if (model == 'spatial') {
  fixedk_spatial_mcmc <- mc.spatial
}

if (model == 'BASCS') {
  fixedk_bascs_mcmc <- mc.bext
}

if (spectral_model == 'extended') {
  source('spectral_post_extended_full.R')
}


setwd(home)
#############################
# Initial burn-in MCMC run
#############################

online_ordering <- "reference"
no_guess <- k_curr
mu_guess <- mu_curr
adapt_end <- 0  # Iteration number on which to stop adaptive MCMC location updates

paras_burnin <- vector('list', main_mcmc_iters)
store_alloc_burnin <- allocate_curr*0

for (tt in 1:burnin_mcmc_iters) {
  
  
  if (tt/print_interval == round(tt/print_interval)){
    print(paste0("Burn-In Iteration number: ",tt,'/',burnin_mcmc_iters))
    print(paste0('w: ', w))
  }
  
  if (model == 'eBASCS') {
    burnin_run_time <- fixedk_eBASCS_mcmc(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,eparas_all,ewt_all,k_curr,mix_num,bk,num_time_breaks, lambda, time_bin)
    
    store_alloc_burnin = store_alloc_burnin + burnin_run_time[[2]]
    
    # Extract parameters and allocation matrix after burnin
    new_paras <- burnin_run_time[[1]]
    k_curr <- new_paras[1]
    mix_num <- k_curr+1
    mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
    w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
    
    eparas_all <- burnin_run_time[[3]] 
    ewt_all <- burnin_run_time[[4]] 
    eparas_all <- burnin_run_time[[3]] 
    ewt_all <- burnin_run_time[[4]] 
    lambda <- burnin_run_time[[5]] 
    time_bin <- burnin_run_time[[6]]
    bk <- burnin_run_time[[7]]
    num_time_breaks <- burnin_run_time[[8]]
    
    allocate_curr <- burnin_run_time[[2]]
    
    
    predicted = factor((allocate_curr %*% c(1, 2, 3))[,1], levels = 1:3)
  }
  
  if (model == 'spacetime') {
    burnin_run_time <- fixedk_spacetime_mcmc(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,k_curr,mix_num,bk,num_time_breaks, lambda, time_bin)
    
    # extract allocation and add it to the overall allocation to get allocation rates
    store_alloc_burnin = store_alloc_burnin + burnin_run_time[[2]]
    
    # Extract parameters and allocation matrix after burnin
    new_paras <- burnin_run_time[[1]]
    k_curr <- new_paras[1]
    mix_num <- k_curr+1
    mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
    w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
    
    lambda <- burnin_run_time[[3]] 
    time_bin <- burnin_run_time[[4]]
    bk <- burnin_run_time[[5]]
    num_time_breaks <- burnin_run_time[[6]]
    
    
    allocate_curr <- burnin_run_time[[2]]
    accept_mu <- burnin_run_time[[7]]
    
    
    predicted = factor((allocate_curr %*% c(1:3))[,1], levels = 1:3)
    
  }
  
  if (model == 'BASCS') {
    burnin_run_time <- fixedk_bascs_mcmc(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,eparas_all,ewt_all,k_curr,mix_num)
    # extract allocation and add it to the overall allocation to get allocation rates
    store_alloc_burnin = store_alloc_burnin + burnin_run_time[[2]]
    
    # Extract parameters and allocation matrix after burnin
    new_paras <- burnin_run_time[[1]]
    k_curr <- new_paras[1]
    mix_num <- k_curr+1
    mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
    w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
    eparas <- burnin_run_time[[3]]
    ewt <- burnin_run_time[[4]]
    allocate_curr <- burnin_run_time[[2]]
    
    predicted = factor((allocate_curr %*% c(1, 2, 3))[,1], levels = 1:3)
  }
  
  if (model == 'spatial') {
     burnin_run_time <-fixedk_spatial_mcmc(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,k_curr,mix_num)
     
     # extract allocation and add it to the overall allocation to get allocation rates
     store_alloc_burnin = store_alloc_burnin + burnin_run_time[[2]]
     
     # Extract parameters and allocation matrix after burnin
     new_paras <- burnin_run_time[[1]]
     k_curr <- new_paras[1]
     mix_num <- k_curr+1
     mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
     w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
     
     allocate_curr <- burnin_run_time[[2]]
     accept_mu <- burnin_run_time[[3]]
     count_vector <- burnin_run_time[[4]]
     
     predicted = factor((allocate_curr %*% c(1:(k_curr+1)))[,1], levels = 1:k_curr+1)
   }
  
}

#############################
# Main MCMC run
#############################

online_ordering <- "reference"
no_guess <- k_curr
mu_guess <- mu_curr[,order(w[-mix_num],decreasing=TRUE)]
#mu_guess <- mu_curr
adapt_end <- 0  # Iteration number on which to stop adaptive MCMC location updates

paras <- vector('list', main_mcmc_iters)
store_alloc <- allocate_curr*0

accept_mu <- c(0,0)

for (tt in 1:main_mcmc_iters){
  
  # Print iteration at intervals
  if (tt/print_interval == round(tt/print_interval)){
    print(paste0("Main Run Iteration number: ",tt,'/',main_mcmc_iters))
    print(paste0('acceptance: ', accept_mu))
    print(paste0('w: ',w))
  }
  
  if (model == 'eBASCS') {
    main_run_time <- fixedk_eBASCS_mcmc(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,eparas_all,ewt_all,k_curr,mix_num,bk,num_time_breaks, lambda, time_bin)
    
    
    # extract allocation and add it to the overall allocation to get allocation rates
    store_alloc = store_alloc + main_run_time[[2]]
    
    # Extract parameters and allocation matrix after burnin
    new_paras <- main_run_time[[1]]
    k_curr <- new_paras[1]
    mix_num <- k_curr+1
    mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
    w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
    
    eparas_all <- main_run_time[[3]] 
    ewt_all <- main_run_time[[4]] 
    eparas_all <- main_run_time[[3]] 
    ewt_all <- main_run_time[[4]] 
    lambda <- main_run_time[[5]] 
    time_bin <- main_run_time[[6]]
    bk <- main_run_time[[7]]
    num_time_breaks <- main_run_time[[8]]
    
    allocate_curr <- main_run_time[[2]]
    accept_mu <- main_run_time[[9]]
    
    
    predicted = factor((allocate_curr %*% c(1:(k_curr+1)))[,1], levels = (1:(k_curr+1)))
    #confusion = table(predicted, truth = dat[,5])
    
    #paras[[tt]] <- list(par = new_paras, eparas_all = eparas_all, ewt_all = ewt_all,lambda = lambda, confusion = confusion)
    paras[[tt]] <- list(par = new_paras, eparas_all = eparas_all, ewt_all = ewt_all,lambda = lambda, predicted = predicted)
  }
  
  if (model == 'spacetime') {
    main_run_time <- fixedk_spacetime_mcmc(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,k_curr,mix_num,bk,num_time_breaks, lambda, time_bin)
    
    
    # extract allocation and add it to the overall allocation to get allocation rates
    store_alloc = store_alloc + main_run_time[[2]]
    
    # Extract parameters and allocation matrix after burnin
    new_paras <- main_run_time[[1]]
    k_curr <- new_paras[1]
    mix_num <- k_curr+1
    mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
    w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
    
    lambda <- main_run_time[[3]] 
    time_bin <- main_run_time[[4]]
    bk <- main_run_time[[5]]
    num_time_breaks <- main_run_time[[6]]
    
    allocate_curr <- main_run_time[[2]]
    accept_mu <- main_run_time[[7]]
    
    predicted = factor((allocate_curr %*% c(1:(k_curr+1)))[,1], levels = (1:(k_curr+1)))
    paras[[tt]] <- list(par = new_paras, lambda = lambda, predicted = predicted)
  }
  
  if (model == 'BASCS') {
    main_run_time <- fixedk_bascs_mcmc(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,eparas_all,ewt_all,k_curr,mix_num)
    
    
    # extract allocation and add it to the overall allocation to get allocation rates
    store_alloc = store_alloc + main_run_time[[2]]
    
    # Extract parameters and allocation matrix after burnin
    new_paras <- main_run_time[[1]]
    k_curr <- new_paras[1]
    mix_num <- k_curr+1
    mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
    w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
    eparas_all <- main_run_time[[3]]
    ewt_all <- main_run_time[[4]]
    allocate_curr <- main_run_time[[2]]
    
    predicted = factor((allocate_curr %*% c(1:(k_curr+1)))[,1], levels = 1:(k_curr+1))
    
    paras[[tt]] <- list(par = new_paras, eparas_all = eparas_all, ewt_all = ewt_all, predicted = predicted)
  }
  
  if (model == 'spatial') {
    main_run_time <-fixedk_spatial_mcmc(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,k_curr,mix_num)
    
    # extract allocation and add it to the overall allocation to get allocation rates
    store_alloc = store_alloc + main_run_time[[2]]
    
    # Extract parameters and allocation matrix after burnin
    new_paras <- main_run_time[[1]]
    k_curr <- new_paras[1]
    mix_num <- k_curr+1
    mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
    w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
    
    allocate_curr <- main_run_time[[2]]
    accept_mu <- main_run_time[[3]]
    count_vector <- main_run_time[[4]]
    
    predicted = factor((allocate_curr %*% c(1:(k_curr+1)))[,1], levels = 1:(k_curr+1))
    paras[[tt]] <- list(new_paras, allocate_curr)
  }
}

main_run_time <- list(paras, allocate_curr, alloc_rate = store_alloc/main_mcmc_iters)

# Save output
save(list=ls(),file=paste(results,out.file.name,sep=""))

# Show time taken 
end <- proc.time()
print(end-start)
