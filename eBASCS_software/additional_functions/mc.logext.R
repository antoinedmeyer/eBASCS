#--------------------------------------------------------------------------#
#VERSION OF MCMC THAT ALLOWS FOR CONSTANT SPECTRUM ACROSS TIME AND SOURCES
#HENCE ONLY ONE SPECTRAL MODEL PER SOURCE IS FITTED
#--------------------------------------------------------------------------#
no_guess <-2
# mu_guess <- mu_curr

mc.logext <- function(fix_runs,online_ordering,rjmcmc_run,w,allocate_curr,mu_curr,eparas_all,ewt_all,k_curr,mix_num,bk,num_time_breaks, lambda, time_bin){
  
  # Number of time breaks is given for all sources
  if(length(num_time_breaks) == 1){
    num_time_breaks = rep(num_time_breaks, k_curr)
  }
  # # checks
  # length(num_time_breaks) == k_curr
  # length(bk) == num_time_breaks + 1
  
  # Standard MCMC updates
  for (t2 in 1:fix_runs){
    # Update photon allocations 
    probs <- matrix(NA,mix_num,obs_num)
    
    for(i in 1:k_curr){
      dlambda = lambda[[i]][time_bin[[i]]]
      # calculate likelihoods
      #dpsf = psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial,mu_curr[,i])
      dpsf = psf(spatial,exp(mu_curr[,i]))
      dE = ewt_all[i]*dgamma(energy,exp(eparas_all[i,3]),exp(eparas_all[i,3])/exp(eparas_all[i,1]))+(1-ewt_all[i])*dgamma(energy,exp(eparas_all[i,4]),exp(eparas_all[i,4])/exp(eparas_all[i,2]))
      #dpsf*dE*dlambda = likelihood given assignment s, and w = p(s). Hence probs[,i] is posterior probability
      #that photons belong to source i.
      probs[i,] = w[i]*dpsf*dE*dlambda
    }
    probs[is.na(probs)] <- 0
    
    # update background probabilities
    # subset background lambda accourding to photon time arrival then if enery is bounded above by the max allowed energy.
    dlambda_back = lambda[[mix_num]][time_bin[[mix_num]]][energy <= max_back_energy]
    
    # prob(background photon) = (mixture weight) * (psf equivalent unif) * (uniform energy dist'n) * (relative time intensity)
    probs[mix_num,energy <= max_back_energy] <- w[mix_num]*(1/img_area)*(1/max_back_energy)*dlambda_back
    probs[mix_num,energy > max_back_energy] <- 0
    
    #allocate the photons given the probabilities.  
    allocate_curr <- t(matrix(unlist(lapply(1:obs_num, function(i) rmultinom(1, 1, probs[,i]))),ncol=obs_num))  # Don't need to normalize as rmultinom function does it automatically
    
    # Counts vector
    mix_num <- k_curr+1
    count_vector <- matrix(NA,mix_num,1)
    count_vector[1:mix_num] <- apply(allocate_curr[,1:mix_num],2,sum)
    
    # Update positions
    mu_prop <- mu_curr
    for (i in 1:k_curr){      
      index <- allocate_curr[,i]==1
      if (count_vector[i]>0){
        # Adaptive version (eventually ended to ensure convegence)
        if (rjmcmc_run < adapt_end){
          mu_prop[,i] <- rnorm(2,mu_curr[,i],mu_adapt_prop_sd/sqrt(count_vector[i]))
        }
        # Non-adaptive version
        if (rjmcmc_run >= adapt_end){
          mu_prop[,i] <- rnorm(2,mu_curr[,i],mu_fixed_prop_sd)
        }
        logr <- sum(log(psf(spatial[index,],exp(mu_prop[,i]))))+sum(mu_prop[,i])-(sum(log(psf(spatial[index,],exp(mu_curr[,i]))))+sum(mu_curr[,i]))
        u <- runif(1,0,1)
        if(is.na(logr)==0){
          if (log(u) < logr){
            mu_curr[,i] <- mu_prop[,i]
            accept_mu[i] <- accept_mu[i] +1 
          }
        }
      }
      # Have to make sure that sources without phoons assigned move (doesn't effect likelihood) 
      if (count_vector[i]==0){
        if (rjmcmc_run < adapt_end){
          mu_curr[,i] <- c(runif(1,log(min(dat[dat$x > 0,]$x)),log(max(dat$x))),
                           runif(1,log(min(dat[dat$y > 0,]$y)),log(max(dat$y))))
        } else {
          mu_curr[,i] <- rnorm(2,mu_curr[,i],mu_fixed_prop_sd)
        }
      }
    }
    
    # Order parameters by suspected sources intensities (associated with particular positions)
    if (k_curr > 1 & online_ordering =="reference"){
      to_order <- min(no_guess,k_curr)
      next_index <- which.min(apply((mu_curr-mu_guess[,1])^2,2,sum))
      next_index_store <- next_index
      for (i in 2:to_order){
        next_order <- order(apply((mu_curr-mu_guess[,i])^2,2,sum))
        next_index_all <- setdiff(next_order,next_index_store)
        next_index_store <- c(next_index_store,next_index_all[1])
      }
      indexmu <- c(next_index_store,setdiff(1:k_curr,next_index_store))
      mu_curr <- mu_curr[,indexmu]
      count_vector[1:k_curr] <- count_vector[indexmu]
      allocate_curr[,1:k_curr] <- allocate_curr[,indexmu]
      eparas_all <- eparas_all[indexmu,]
      ewt_all <- ewt_all[indexmu]
      num_time_breaks[1:k_curr] <- num_time_breaks[indexmu]
      bk[1:k_curr] <- bk[indexmu]
      time_bin[1:k_curr] <- time_bin[indexmu]
    }
    
    # Update weights
    alpha <- rep(wprior,mix_num)
    w <- rdirichlet(1,alpha+count_vector)
    
    # Update lambda weights
    for(i in 1:mix_num){
      lambda0 <- rep(lambdaprior,num_time_breaks[i])
      count_vector_lambda = table(cut(arrival_time[which(allocate_curr[,i] == 1)], bk[[i]]))
      lambda[[i]] = rdirichlet(1,lambda0 + count_vector_lambda)
    }
    
    # Update spectral parameters (full model) 
    # Update spectral parameters (extended full model) 
    for (i in 1:k_curr){
      index <- which(allocate_curr[,i] == 1)
      cspatial <- energy[index]
      gm1curr <- eparas_all[i,1]
      gm2curr <- eparas_all[i,2]
      ga1curr <- eparas_all[i,3]
      ga2curr <- eparas_all[i,4]
      ewtcurr <- ewt_all[i]
      ewtprop <- inv.logit(rnorm(1,logit(ewtcurr),specwt_sd))
      gm1prop <- rnorm(1,gm1curr,i*specm_sd)
      ga1prop <- rnorm(1,ga1curr,speca_sd)
      if ((exp(gm1prop) > emean.min) & (exp(gm1prop) < emean.max) & (exp(ga1prop) > 0)){
        logr <- spectral_post(exp(gm1prop),exp(ga1prop),exp(gm2curr),exp(ga2curr),ewtprop,cspatial) - spectral_post(exp(gm1curr),exp(ga1curr),exp(gm2curr),exp(ga2curr),ewtcurr,cspatial)
        u <- runif(1,0,1)
        if (log(u) < logr){
          eparas_all[i,c(1,3)] <- c(gm1prop,ga1prop)
          ewt_all[i] <- ewtprop
          gm1curr <- eparas_all[i,1]
          ga1curr <- eparas_all[i,3]
          ewtcurr <- ewt_all[i]
        }
      }
      gm2prop <- rnorm(1,gm2curr,i*specm_sd)
      ga2prop <- rnorm(1,ga2curr,speca_sd)
      if ((exp(gm2prop) > emean.min) & (exp(gm2prop) < emean.max) & (exp(ga2prop) > 0)){
        logr <- spectral_post(exp(gm1curr),exp(ga1curr),exp(gm2prop),exp(ga2prop),ewtcurr,cspatial)- spectral_post(exp(gm1curr),exp(ga1curr),exp(gm2curr),exp(ga2curr),ewtcurr,cspatial)
        u <- runif(1,0,1)
        if (log(u) < logr){
          eparas_all[i,c(2,4)] <- c(gm2prop,ga2prop)
        }
      }
    }
    # Order Gammas for identifiability and combine conditions
    if (k_curr > 1){
      epara_order <- apply(eparas_all[,c(1,2)],1,order)
      epara_order_all <- rbind(epara_order,epara_order+2)
      for (i in 1:k_curr){
        eparas_all[i,] <- eparas_all[i,epara_order_all[,i]] 
      }
      for (i in 1:k_curr){
        if (sum(epara_order[,i] != c(1,2))>0){
          ewt_all[i] <- 1- ewt_all[i]
        }
      }
    }
  }
  
  # Output parameters and log-posterior
  value <- list(c(k_curr,c(mu_curr),c(w)),allocate_curr, eparas_all, ewt_all,lambda, time_bin, bk, num_time_breaks, accept_mu)
  return(value)
}

