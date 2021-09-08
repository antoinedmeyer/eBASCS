# Format data for MCMC input and extract basic information

# Spatial Data
spatial <- dat[,1:2]
spatial <- as.matrix(spatial)
colnames(spatial) <- NULL
rownames(spatial) <- NULL
obs_num <- length(spatial[,1])
xlow <- min(spatial[,1])
xup <- max(spatial[,1])
ylow <- min(spatial[,2])
yup <- max(spatial[,2])
yl <- (yup - ylow)/2
xl <- (xup - xlow)/2
img_area <- 4*xl*yl  # Image area - rectangle assumed

# Spectral Data
if (model == 'BASCS' | model == 'eBASCS') {
  energy <- dat$energy
  max_back_energy <- max(energy)
}

#Temporal Data
if (model == 'eBASCS' | model == 'spacetime') {
  arrival_time <- dat$t
} 


