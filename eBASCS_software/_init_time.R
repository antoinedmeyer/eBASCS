source(paste(home, "/_initialization_time.R",sep=""))
eng.setting <- "const"
k_curr <- theta  # Number of sources
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Initialize breakpoints and time_bins
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

bk = list()
time_bin = list()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Time spacing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

if(time_spacing == 'equal'){
  
  for(i in 1:(k_curr + 1)){
    bk[[i]] = seq(min(arrival_time), max(arrival_time), length.out = num_time_breaks + 1)
    time_bin[[i]] = cut(arrival_time, bk[[i]], include.lowest = TRUE)
    levels(time_bin[[i]]) = 1:num_time_breaks
    time_bin[[i]] = as.numeric(time_bin[[i]])
  }
  num_time_breaks = rep(num_time_breaks, k_curr + 1)
  
  if (spectral_model == 'full') {
    
    # rearrange initial (alpha, gamma)
    eparas_all_tmp = list()
    lambda_tmp = list()
    for(i in 1:(k_curr)){
      if (eng.setting=="equalE") {
        eparas_all_tmp[[i]] = list()
        for(k in 1:num_time_breaks[i]){
          eparas_all_tmp[[i]][[k]] = matrix(eparas_all[[k]][i,], nrow = 2, ncol = 1)
        }
      }
      if (eng.setting=="const") {
        eparas_all_tmp[[i]] = list()
        eparas_all_tmp[[i]] = matrix(eparas_all[i,], nrow = 2, ncol = 1)
      }
      if (eng.setting=='notconst') {
        eparas_all_tmp[[i]] = list()
        for(k in 1:num_time_breaks[i]){
          eparas_all_tmp[[i]][[k]] = matrix(eparas_all[[k]][i,], nrow = 2, ncol = 1)
        }
      }
    }
    eparas_all = eparas_all_tmp
    # initialize lambda
    # for(i in 1:(k_curr+1)){
    #   lambda_tmp[[i]] = prop.table(table(time_bin[[i]][allocate_curr[,i]==1]))
    # }
    
    eparas_all = eparas_all_tmp
    # lambda = lambda_tmp
  }
  
}

if(time_spacing == 'custom.HBC515'){
  
  num_time_breaks = 3
  
  for(i in 1:(k_curr + 1)){
    #bk[[i]] = seq(min(arrival_time), max(arrival_time), length.out = num_time_breaks + 1)
    bk[[i]] = c(0,5000,20000,max(arrival_time))
    time_bin[[i]] = cut(arrival_time, bk[[i]], include.lowest = TRUE)
    levels(time_bin[[i]]) = 1:num_time_breaks
    time_bin[[i]] = as.numeric(time_bin[[i]])
  }
  num_time_breaks = rep(num_time_breaks, k_curr + 1)
  
  # rearrange initial (alpha, gamma)
  eparas_all_tmp = list()
  lambda_tmp = list()
  for(i in 1:(k_curr)){
    if (eng.setting=="equalE") {
      eparas_all_tmp[[i]] = list()
      for(k in 1:num_time_breaks[i]){
        eparas_all_tmp[[i]][[k]] = matrix(eparas_all[[k]][i,], nrow = 2, ncol = 1)
      }
    }
    if (eng.setting=="const") {
      eparas_all_tmp[[i]] = list()
      eparas_all_tmp[[i]] = matrix(eparas_all[i,], nrow = 2, ncol = 1)
    }
  }
  eparas_all = eparas_all_tmp
  # initialize lambda
  # for(i in 1:(k_curr+1)){
  #   lambda_tmp[[i]] = prop.table(table(time_bin[[i]][allocate_curr[,i]==1]))
  # }
  
  eparas_all = eparas_all_tmp
  # lambda = lambda_tmp
  
}

if (time_spacing == 'algo.UVCet') {
  segmentation <- function(dat,breaks,max) {
    bins <- cut(dat$t, breaks = breaks,labels = FALSE)
    data <- cbind(dat, bins)
    bin_values <- c()
    for (i in 1:breaks) {
      bin_values[i] <- nrow(data[data$bins == i,])
    }
    #Compute first order differences
    fod <- abs(bin_values[2:breaks] - bin_values[1:breaks-1])
    #Get topmost values
    maxes <- tail(sort(fod),max)
    #Get corresponding time bins
    indices <- sort(match(maxes,fod))
    #Delete middle ones
    for (j in 1:(length(indices)-1)) {
      if(j < length(indices)) {
        if(indices[j+1] == (indices[j]+1)) {
          k <- j+1
          indices <- indices[-k]
        }
      }
      else {
        break
      }
    }
    breakpoints <- c()
    for (l in 1:length(indices)) {
      # if ((indices[l+1] - indices[l]) <= breaks/10) {
      #   breakpoints <- c(breakpoints, min(data[data$bins==indices[l],]$t),max(data[data$bins==indices[l+1],]$t))
      #   l <- l+1
      # }
      # else {
      #   breakpoints <- c(breakpoints, min(data[data$bins==indices[l],]$t))
      #   l <- l+1
      # }
      if ((l %% 2) == 1) {
        breakpoints <- c(breakpoints, min(data[data$bins==indices[l],]$t))
      }
      else {
        breakpoints <- c(breakpoints, max(data[data$bins==indices[l],]$t))
      }
    }
    hist(dat$t, breaks = breaks, xlab = 'Time', main = '')
    for (i in 1:length(breakpoints)){
      abline(v = breakpoints[i] , col = 'red', lwd = 2)
    }
    return(breakpoints)
  }
  num_time_breaks = length(segmentation(dat,50,8))+1
  for(i in 1:(k_curr + 1)){
    bk[[i]] = c(0,segmentation(dat,50,8),max(arrival_time))
    time_bin[[i]] = cut(arrival_time, bk[[i]], include.lowest = TRUE)
    levels(time_bin[[i]]) = 1:num_time_breaks
    time_bin[[i]] = as.numeric(time_bin[[i]])
  }
  num_time_breaks = rep(length(segmentation(dat,50,8))+1,k_curr+1)
}


initialize_list <- initialization.time(k_curr, num_time_breaks)  

mix_num <- k_curr+1  # Number of mixture components (sources + background)
mu_curr <- log(t(initialize_list[[1]]))  # Locations of sources (2 x k_curr matrix) - locations initially on edge of image
w <- initialize_list[[2]]  # Relative intensities (vector of length mix_num)
allocate_curr <- initialize_list[[3]]  # Allocation of photons to sources (and background) - initially obs_num x mix_num matrix of zeros


if (model == 'BASCS' || model == 'eBASCS') {
  eparas_all <- log(initialize_list[[4]])  # Shape and mean spectral parameters - disregard if spectral_model=="none"
  #eparas_all <- initialize_list[[4]]
  ewt_all <- initialize_list[[5]]  # Spectral model weights (extended full model) - disregard unless spectral_model=="extended_full"
}

if (model == 'spacetime' || model == 'eBASCS') {
  lambda <- initialize_list[[6]]  # Relative Intensities of time arrival distribution
}
