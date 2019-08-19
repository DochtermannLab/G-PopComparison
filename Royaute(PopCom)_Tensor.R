library(gdata);library(matrixcalc);library(MCMCglmm);library(magrittr)

########################################################################
# Eigentensor analysis following code provided by Aguire et al. (2014) #
########################################################################

#### 0. Storage of MCMC samples ####
#number of MCMC samples
MCMCsamp <- 1000 

#number of traits 
n <- 7 

#number of matrices to compare
m <- 4



#number of random effects specified in the model. In our analyses these were animal and residual effects.
r <- 2

#trait names
traitnames <- c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo")

#matrix labels 
Gnames <- c("G_Dunn","G_AG","G_SOC","G_LC")

# Import Covariance matrices posterior values
G.MCMC.Dunn=as.matrix(Gmat.MCMC.Dunn$VCV)
G.MCMC.AG=as.matrix(Gmat.MCMC.AG$VCV)
G.MCMC.SOC=as.matrix(Gmat.MCMC.SOC$VCV)
G.MCMC.LC=as.matrix(Gmat.MCMC.LC$VCV)


MCMCarray <- array(,c(MCMCsamp,(n^2)*r,m)) 
#empty array 

#G_Dunn stored as the 1st element of dim[4] 
MCMCarray[,,1] <- G.MCMC.Dunn

#G_AG stored as the 2nd element of dim[4]
MCMCarray[,,2] <- G.MCMC.AG

#G_SOC stored as the 3rd element of dim[4]
MCMCarray[,,3] <- G.MCMC.SOC

#G_LC stored as the 4th element of dim[4]
MCMCarray[,,4] <- G.MCMC.LC

# Reshaping the array and standardizing G
Garray <- array(,c(n,n,m,MCMCsamp))
dimnames(Garray) <- list(traitnames,traitnames,Gnames)

Parray <- array(,c(n,n,m,MCMCsamp))
dimnames(Parray) <- list(traitnames,traitnames,Gnames)

for (i in 1:m){
  for (j in 1:MCMCsamp){
    G <- matrix(MCMCarray[j,1:(n^2),i],ncol= n)
    R <- matrix(MCMCarray[j,((n^2)+1):((n^2)*2),i],ncol= n)
    Garray[,,i,j] <- G
    Parray[,,i,j] <- G + R
  }
}

Garray[1,1,3,5] # Va of first trait (Latency, element 1,1) and first population (ABQ)
Garray[,,2,13] # Gmatrix for second population (AG) and 13th MCMC sample
Garray[,,,42] # Gmatrices of all populations and 42nd MCMC sample

# Apply multivariate standardization (Gp = P^(-1/2) * G * P^(-1/2)) to the Garray
# Step 1: function to calculate P^(-1/2)
inv.rootP <- function (P){
  rootP <- matrix(0,n, n)  
  for (i in 1:n){
    val <- eigen(P)$values
    vec <- eigen(P)$vectors
    rootP <- rootP + (vec[,i] %*% t(vec[,i]))*sqrt(val[i])
  }
  solve(rootP)
}

# Step 2: Apply function to each MCMC sample of each population
HHGarray <- array(,c(n,n,m,MCMCsamp)) # Store standardized G arrays
for (k in 1:MCMCsamp){
  for (j in 1:m){
    P <- inv.rootP(Parray[,,j,k])
    HHGarray[,,j,k] <- P %*% Garray[,,j,k] %*% P #Standardized G 
  }
}

# Randomized G matrices for hypothesis tests
rand.Garray <- array(,c(n,n,m,MCMCsamp))
dimnames(rand.Garray) <- list(traitnames,traitnames,Gnames)
for (i in 1:MCMCsamp){
  G1.bv<-rbv(PooledPedigree.analysis.Dunn,Garray[,,1,i])
  G2.bv<-rbv(PooledPedigree.analysis.AG,Garray[,,2,i])
  G3.bv<-rbv(PooledPedigree.analysis.SOC,Garray[,,3,i])
  G4.bv<-rbv(PooledPedigree.analysis.LC,Garray[,,4,i])
  a.pop <- cumsum(c(dim(PooledPedigree.analysis.Dunn)[1],
                    dim(PooledPedigree.analysis.AG)[1],
                    dim(PooledPedigree.analysis.SOC)[1],
                    dim(PooledPedigree.analysis.LC)[1]))
  pop.bv <- rbind(G1.bv,G2.bv,G3.bv,G4.bv)
  rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1],replace=F),]
  rand.Garray[,,1,i] <- cov(rand.pop.bv[1:a.pop[1],])
  rand.Garray[,,2,i] <- cov(rand.pop.bv[(a.pop[1] + 1):a.pop[2],])
  rand.Garray[,,3,i] <- cov(rand.pop.bv[(a.pop[2] + 1):a.pop[3],])
  rand.Garray[,,4,i] <- cov(rand.pop.bv[(a.pop[3] + 1):a.pop[4],])
  }

#Function to do projection
proj<- function(G, b) t(b) %*% G %*% (b)

#### 0.1 Covariance tensor function ####
covtensor <- function(Gs){
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]] 
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  neigten <- n*(n+1)/2 
  #Number of eigentensors
  MCMC.S <- array(,c(neigten, neigten, MCMCsamp))
  dimnames(MCMC.S) <- list(paste("e", 1:neigten, sep=""), paste("e", 1:neigten, sep=""))
  for (k in 1:MCMCsamp){
    MCMCG <- Gs[,,,k] 
    MCMCvarmat <- t(apply(MCMCG, 3, diag)) 
    #find the variances of the kth G and store them 
    MCMCcovmat <- t(apply(MCMCG, 3, lowerTriangle)) 
    #find the covariances of the kth G and store them
    MCMC.S[1:n,1:n, k] <- cov(MCMCvarmat, MCMCvarmat) 
    #fill the upper left quadrant of the kth S
    MCMC.S[(n+1):neigten,(n+1):neigten, k] <- 2*cov(MCMCcovmat, MCMCcovmat)
    #fill the lower right quadrant of the kth S
    MCMC.S[1:n,(n+1):neigten, k] <- sqrt(2)*cov(MCMCvarmat, MCMCcovmat)
    #fill the upper right quadrant of the kth S
    MCMC.S[(n+1):neigten,1:n, k] <- sqrt(2)*cov(MCMCcovmat, MCMCvarmat)
    #fill the lower left quadrant of the kthS
  }  
  av.S <- apply(MCMC.S, 1:2, mean)
  #posterior mean S
  av.S.val <- eigen(av.S)$values
  #eigenvalues of posterior mean S 
  av.S.vec <- eigen(av.S)$vectors
  #eigenvalues of posterior mean S
  eTmat <- array(, c(n, n, neigten))
  dimnames(eTmat) <- list(traitnames, traitnames, paste("E", 1:neigten, sep=""))  
  for (i in 1:neigten){
    emat <- matrix(0, n, n) 
    lowerTriangle(emat) <- 1/sqrt(2)*av.S.vec[(n+1):neigten,i]
    emat <- emat + t(emat)
    diag(emat) <- av.S.vec[1:n,i]
    eTmat[,,i] <- emat 
  }
  #construct the second-order eigentensors of posterior mean S
  eT.eigen <- array(, c(n+1, n, neigten))
  for (i in 1:neigten){
    eT.eigen[1,,i] <- t(eigen(eTmat[,,i])$values) 
    #Eigenvalues of the ith eigentensor
    eT.eigen[2:(n+1),,i] <- eigen(eTmat[,,i])$vectors 
    #Eigenvectors of the ith eigentensor
    eT.eigen[,,i] <- eT.eigen[,order(abs(eT.eigen[1,,i]), decreasing = T), i]
  }
  MCMC.S.val <- matrix(, MCMCsamp, neigten)
  colnames(MCMC.S.val) <- paste("E", 1:neigten, sep="")
  for (i in 1:MCMCsamp){
    for(j in 1:neigten){
      MCMC.S.val[i,j] <- t(av.S.vec[,j]) %*% MCMC.S[,,i] %*% av.S.vec[,j]
    }
  }
  #posterior distribution of the genetic variance for the eigenvectors of posterior mean S
  av.G.coord <- array(, c(m, neigten, 1))
  dimnames(av.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    av.G.coord[,i,] <- apply((apply(Gs, 1:3, mean)) , 3, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the jth avG for the eigentensors of posterior mean S
  MCMC.G.coord <- array(, c(m, neigten, MCMCsamp))
  dimnames(MCMC.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    MCMC.G.coord[,i,] <- apply(Gs, 3:4, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the kth MCMC sample of the jth G for the eigentensors of posterior mean S
  tensor.summary <- data.frame(rep(av.S.val,each=n), t(data.frame(eT.eigen)))
  colnames(tensor.summary) <- c("S.eigval", "eT.val", traitnames)
  rownames(tensor.summary)<- paste(paste("e", rep(1:neigten, each=n), sep=""), rep(1:n,neigten), sep=".")
  list(tensor.summary = tensor.summary, av.S = av.S, eTmat = eTmat, av.G.coord = av.G.coord, MCMC.S = MCMC.S, MCMC.S.val = MCMC.S.val, MCMC.G.coord = MCMC.G.coord)
}


##### 1. Covariance Tensor results ####
MCMC.covtensor <- covtensor(Garray)
nnonzero <- min(n*(n+1)/2,m-1)
MCMC.covtensor.rand <- covtensor(rand.Garray)
HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero]), prob=0.95), 
                    HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero]), prob=0.95))
mean.eT.val <- cbind(colMeans(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero])), 
                     colMeans(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero])))
round(mode.eT.val, 3)
round(HPD.eT.val, 3)

#### 2. Significance test ####
# E1
quantile(as.mcmc(MCMC.covtensor$MCMC.S.val[,1]), c(.0275,.975))
quantile(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1]), c(.0275,.975)) #overlap
quantile(as.mcmc(MCMC.covtensor$MCMC.S.val[,1]), c(.05,.95))
quantile(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1]), c(.05,.95)) #Pmcmc=0.9
quantile(as.mcmc(MCMC.covtensor$MCMC.S.val[,1]), c(.075,.925))
quantile(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1]), c(.075,.925)) #Pmcmc=0.85

#E2
quantile(as.mcmc(MCMC.covtensor$MCMC.S.val[,2]), c(.0275,.975))
quantile(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,2]), c(.0275,.975)) #overlap
quantile(as.mcmc(MCMC.covtensor$MCMC.S.val[,2]), c(.05,.95))
quantile(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,2]), c(.05,.95)) #overlap 0.9
quantile(as.mcmc(MCMC.covtensor$MCMC.S.val[,2]), c(.075,.925))
quantile(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,2]), c(.075,.925)) #Pmcmc=0.85

#E3
quantile(as.mcmc(MCMC.covtensor$MCMC.S.val[,3]), c(.0275,.975))
quantile(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,3]), c(.0275,.975)) #overlap
quantile(as.mcmc(MCMC.covtensor$MCMC.S.val[,3]), c(.05,.95))
quantile(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,3]), c(.05,.95)) #Pmcmc=0.9
quantile(as.mcmc(MCMC.covtensor$MCMC.S.val[,3]), c(.075,.925))
quantile(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,3]), c(.075,.925)) #Pmcmc=0.85


#### 3. Alignement with gmax (Table Sx) ####
# Critical values of random distribution
quantile(rand.vec.cor,c(.9,.95,0.99,.999))




