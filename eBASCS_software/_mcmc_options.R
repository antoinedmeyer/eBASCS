# Choose options for MCMC

# Number of Iterations and Thinning
burnin_mcmc_iters <- 500  # Number of RJMCMC iterations for the burnin period
main_mcmc_iters <- 4000 # Number of RJMCMC iterations for the main run

# Thinning:
mcmc_runs <- 10  # Number of steps of MCMC to run between each recording or parameters
print_interval <- 50  # Number of iterations between print of current iteration number

# Adaptive location updates
adapt_end <- 1500  # Iteration number on which to stop adaptive MCMC location updates
mu_adapt_prop_sd <- 1e-2  # Adaptive MCMC location updates (first adapt_end iterations):
# the standard deviation of the proposal distribution for 
# updating the location of source i is mu_adapt_prop_var/sqrt(ni)
# where ni is the number of photons currently assigned to 
# source i.

# Standard Deviation Location Proposal
mu_fixed_prop_sd <- 1e-3  # Location fixed updates (after the first adapt_end iterations):
# the standard deviation of the proposal distribution for 
# updating the location of a source is mu_fixed_prop_var

# Standard Deviation Energy Distributions Proposal
if (model == 'BASCS' | model == 'eBASCS') {
  specm_sd <- 10  # Proposal: e ~ N(e_current,10*i), where e is the mean parameter of the spectral model, and
  # i is the source number (in the inferred sources are ordered so that source 1 is the brightest).
  # The proposal is rejected if the proposed value of e is outside (emean.min,emean.max), see (3). 
  speca_sd <- 1  # Proposal: a ~ N(a_current,speca_sd), where a is the shape parameter of the spectral model
  specwt_sd <- 0.1  # Proposal: ewt ~ inv.logit(rnorm(1,logit(ewt_current),specwt_sd)), where ewt is the weight parameter
  # in the extended full model (whcih uses a two gamma spectral model)
}

# Initialization for Location Parameters. Options: "modes" (KDE), "top.corners", "rand.corners", ...
location_init <- 'modes'