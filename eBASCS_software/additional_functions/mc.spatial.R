no_guess <-2
#mu_guess <- mu_curr

mc.spatial <- function(num_samples,online_ordering,rjmcmc_run,w,allocate_curr,mu_curr,k_curr,mix_num){
  
  for (t1 in 1:num_samples){
    
    # Print iteration at intervals
    if (t1/print_interval == round(t1/print_interval)){
      print(paste("Iteration number: ",t1,sep=""))
    }
    
    # Update photon allocations 
    probs <- matrix(NA,mix_num,obs_num)
    #probs[1:k_curr,] <- t(matrix(unlist(lapply(1:k_curr,function(i) w[i]*psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial,mu_curr[,i]))),ncol=k_curr))
    probs[1:k_curr,] <- t(matrix(unlist(lapply(1:k_curr,function(i) w[i]*psf(spatial,exp(mu_curr[,i])))),ncol=k_curr))
    #probs[mix_num,energy <= max_back_energy] <- w[mix_num]*(1/img_area)
    probs[mix_num,] <- w[mix_num]*(1/img_area)
    #probs[mix_num,energy > max_back_energy] <- 0
    probs[is.na(probs)] <- 0
    allocate_curr <- t(matrix(unlist(lapply(1:obs_num,function(i) rmultinom(1, 1, probs[,i]))),ncol=obs_num))  # Don't need to normalize as rmultinom function does it automatically
    
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
      next_index <- which.min(apply((exp(mu_curr)-exp(mu_guess[,1]))^2,2,sum))
      next_index_store <- next_index
      if (to_order > 1){
        for (i in 2:to_order){
          next_order <- order(apply((exp(mu_curr)-exp(mu_guess[,i]))^2,2,sum))
          next_index_all <- setdiff(next_order,next_index_store)
          next_index_store <- c(next_index_store,next_index_all[1])
        }
      }
      indexmu <- c(next_index_store,setdiff(1:k_curr,next_index_store))
      mu_curr <- mu_curr[,indexmu]
      count_vector[1:k_curr] <- count_vector[indexmu]
      allocate_curr[,1:k_curr] <- allocate_curr[,indexmu]
    }
    
    # Update weights
    alpha <- rep(wprior,mix_num)
    w <- rdirichlet(1,alpha+count_vector)
    
  }
  #value <- list(c(k_curr,c(mu_curr),c(w)),allocate_curr,log_posterior_spatial(mu_curr,allocate_curr,w,k_curr))
  value <- list(c(k_curr,c(mu_curr),c(w)),allocate_curr, accept_mu, count_vector)
  return(value)
}
