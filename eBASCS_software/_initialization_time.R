initialization.time <- function(k_curr, num_time_breaks){
  
  mix_num <- k_curr+1
  img_area <- 4*xl*yl
  
  allocate_curr <- matrix(0,obs_num,mix_num)
  
  library(rlist)
  library('ks')
  
  if(location_init == 'modes') {
    modefinder <- function(data) {
      kde <- kde(dat[,1:2], xmin = c(min(dat$x),min(dat$y)), xmax = c(max(dat$x),max(dat$y)))
      split <- floor((35/max(dat$x))*length(kde$eval.points[[1]]))
      m1 <- max(kde$estimate[1:split,])
      m2 <- max(kde$estimate[split:length(kde$eval.points[[1]]),])
      coord1 <- which(kde$estimate[1:split,] == m1, arr.ind = TRUE)
      coord2 <- which(kde$estimate[split:length(kde$eval.points[[1]]),] == m2, arr.ind = TRUE)
      mode1 <- c(kde$eval.points[[1]][1:split][coord1[1]],kde$eval.points[[2]][coord1[2]])
      mode2 <- c(kde$eval.points[[1]][split:length(kde$eval.points[[1]])][coord2[1]],kde$eval.points[[2]][coord2[2]])
      return(list(mode1,mode2))
    }
    m <- modefinder(dat)
    mu_curr <- rbind(m[[1]],m[[2]])
  }
  
  if (location_init == 'top.corners') {
    if (k_curr == 2) {
      #mu_curr <- rbind(c(xlow,yup),c(xup,yup))
      mu_curr <- rbind(c(min(dat[dat$x > 0,]$x), max(dat$y)), c(max(dat$x),max(dat$y)))
    }
    else {
      mu_curr <- rbind(c(xlow,yup),c(xup,yup),c(xup,ylow))
    }
  }
  if (location_init == 'rand.corners') {
    corners <- list(c(xlow,ylow),c(xlow,yup),c(xup,ylow),c(xup,yup))
    samples <- list.sample(corners,2,replace=FALSE)
    mu_curr <- rbind(samples[[1]],samples[[2]])
  }
  if (location_init == 'split.image') {
    mu_curr <- rbind(c(runif(1,xlow,xl),runif(1,ylow,yup)),c(runif(1,xl,xup),runif(1,ylow,yup)))
  }
  if (location_init == 'smart.split') {
    density <- density(dat[,1])
    xsplit <- density$x[match(max(density$y),density$y)]
    mu_curr <- rbind(c(runif(1,xlow,xsplit),runif(1,ylow,yup)),c(runif(1,xsplit,xup),runif(1,ylow,yup)))
  }
  
  w <- rep(1/mix_num,mix_num)
  
  if (model == 'spacetime' || model == 'eBASCS') {
    lambda <- replicate(mix_num, {rep(1/num_time_breaks[1],num_time_breaks[1])}, simplify = FALSE)
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # Energy Parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  if (model == 'BASCS' || model == 'eBASCS') {
    if (eng.setting == 'const') {
      ewt_all <- rep(0.5,k_curr)
      no.eparas <- ifelse(spectral_model == "full",2,4)
      if (spectral_model != "none"){
        eparas_all <- matrix(NA,k_curr,no.eparas)
        for (i in 1:k_curr){
          gmeans <- runif(2,emean.min,emean.max)
          if (spectral_model == "full"){
            eparas_all[i,1] <- gmeans[1]
          } else {
            eparas_all[i,1] <- min(gmeans)
            eparas_all[i,2] <- max(gmeans)
          }
        }
        if (spectral_model == "full"){
          eparas_all[,2] <- 5
        } else {
          eparas_all[,3:4] <- 5
        }
      }
    }
    else {
      energy_parameters_all = replicate(num_time_breaks, {
        ewt <- rep(0.5,k_curr)
        no.eparas <- ifelse(spectral_model == "full",2,4)
        if (spectral_model != "none"){
          eparas <- matrix(NA,k_curr,no.eparas)
          for (i in 1:k_curr){
            gmeans <- runif(2,emean.min,emean.max)
            if (spectral_model == "full"){
              eparas[i,1] <- gmeans[1]
            } else {
              eparas[i,1] <- min(gmeans)
              eparas[i,2] <- max(gmeans)
            }
          }
          if (spectral_model == "full"){
            eparas[,2] <- 5
          } else {
            eparas[,3:4] <- 5
          }
        }
        list(eparas = eparas, ewt = ewt)
      }, simplify = FALSE)
      allocate_curr <- matrix(0,obs_num,mix_num)
      
      eparas_all = lapply(energy_parameters_all, function(x) x[[1]])
      ewt_all = lapply(energy_parameters_all, function(x) x[[2]])
    }
  }
  
  if (model == 'spatial') {
    value <-list(mu_curr,w,allocate_curr,mix_num,img_area)
  }
  if (model == 'spacetime') {
    value <- list(mu_curr,w,allocate_curr,mix_num,img_area,lambda)
  }
  if (model == 'BASCS') {
    value <- list(mu_curr,w,allocate_curr,eparas_all,ewt_all,mix_num,img_area)
  }
  if (model == 'eBASCS') {
    value <- list(mu_curr,w,allocate_curr,eparas_all,ewt_all,lambda,mix_num,img_area)
  }
  return(value)
}

