library(gdata);library(matrixcalc);library(MCMCglmm);library(magrittr)

################################################################################
# Krzanowski Subspace analysis following code provided by Aguire et al. (2014) #
################################################################################

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


#### 0.1 Krzaanowski subspace function ####
#Function to do projection
proj<- function(G, b) t(b) %*% G %*% (b)

kr.subspace <- function(Gs, vec){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]] 
  if(length(vec) != m){stop("vec must have length = m")}
  h <- function (g, v){
    AA <- array(, c(n, n, m))  
    for (k in 1:m){
      g.vec <- eigen(g[,,k])$vectors[,1:(v[k])] 
      AA[,,k] <- g.vec %*% t(g.vec)
    }
    H <- apply(AA, 1:2, sum)
    list(H = H, AA = AA)
  }
  #internal function to calculate AA and H
  MCMC.H <- array(, c(n, n, MCMCsamp))
  dimnames(MCMC.H) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[4]])      
  MCMC.AA <- array(, c(n, n, m, MCMCsamp))
  dimnames(MCMC.AA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]], dimnames(Gs)[[4]])
  for (i in 1:MCMCsamp){
    kr <- h(Gs[,,,i], v = vec)
    MCMC.H[,,i] <- kr$H
    MCMC.AA[,,,i] <- kr$AA
  }	
  #calculate AA and H for the ith MCMC sample of the G array		
  avH <- apply(MCMC.H, 1:2, mean)
  rownames(avH) <- dimnames(Gs)[[1]]
  colnames(avH) <- dimnames(Gs)[[1]]
  #calculate the posterior mean H
  avAA <- apply(MCMC.AA, 1:3, mean)
  dimnames(avAA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]])
  #calculate the posterior mean AA
  avH.vec <- eigen(avH)$vectors
  #eigenanalysis of posterior mean H	
  proj<- function(a, b) t(b) %*% a %*% b
  #internal function to do projection
  avH.theta <- matrix(, n, m)
  for (i in 1:n){
    for (i in 1:n){
      avH.theta[i,] <- acos(sqrt(apply(avAA, 3, proj, b = avH.vec[,i]))) * (180/pi)
    }
  }
  #angles between the eigenvectors posterior mean H and the posterior mean subspaces of each population
  MCMC.H.val <- matrix(, MCMCsamp, n)
  colnames(MCMC.H.val) <- paste("h", 1:n, sep="")
  for (i in 1:n){
    MCMC.H.val[,i] <- apply(MCMC.H, 3, proj, b = avH.vec[,i])
  }
  #posterior distribution of the genetic variance for the eigenvectors of posterior mean H 
  MCMC.H.theta <- array(, c(n, m, MCMCsamp))
  rownames(MCMC.H.theta) <- paste("h", 1:n, sep="")
  colnames(MCMC.H.theta) <- dimnames(Gs)[[3]]
  for(i in 1:n){
    for(j in 1:MCMCsamp){
      MCMC.H.theta[i,,j] <- acos(sqrt(apply(MCMC.AA[,,,j], 3, proj, b = avH.vec[,i]))) * (180/pi)
    }
  }
  #posterior distribution of the angles between the eigenvectors of posterior mean H and the MCMC samples of the subspaces of each population
  list(avAA = avAA, avH = avH, MCMC.AA = MCMC.AA, MCMC.H = MCMC.H, MCMC.H.val = MCMC.H.val, MCMC.H.theta = MCMC.H.theta)
}


#### 1. Krzanowski's common subspaces results ####
MCMCG.kr <- kr.subspace(Garray, vec = rep(n/2,m))
MCMCG.kr.rand <- kr.subspace(rand.Garray, vec = rep(n/2,m))
round(eigen(MCMCG.kr$avH)$vectors, 3)
round(apply(MCMCG.kr$MCMC.H.theta, 1:2, mean), 1)


val <- matrix(, n, m)
for (i in 1:m){
  avG <- apply(Garray, 1:3, mean)
  val[,i] <- round(cumsum(t(eigen(avG[,,i])$values))/sum(eigen(avG[,,i])$values)*100)
}

# Number of subspaces containing 90 % of genetic variation
n.vec <- apply(ifelse(round(val,1) < 90, 1, 0), 2, sum)+1
MCMCG.kr1 <- kr.subspace(Garray, vec = n.vec)
MCMCG.kr.rand1 <- kr.subspace(rand.Garray, vec = n.vec)

round(eigen(MCMCG.kr1$avH)$vectors, 3)
round(apply(MCMCG.kr1$MCMC.H.theta, 1:2, mean), 1)

round(eigen(MCMCG.kr1$avH)$values/sum(eigen(MCMCG.kr1$avH)$values), 5)*100


#### 2. Significance test: calculation of Pmcmc ####
# h1
quantile(as.mcmc(MCMCG.kr1$MCMC.H.val[,1]), c(.0275,.975))
quantile(as.mcmc(MCMCG.kr.rand1$MCMC.H.val[,1]), c(.0275,.975)) #0.95 overlap
quantile(as.mcmc(MCMCG.kr1$MCMC.H.val[,1]), c(.05,.95))
quantile(as.mcmc(MCMCG.kr.rand1$MCMC.H.val[,1]), c(.05,.95)) #0.9 overlap
quantile(as.mcmc(MCMCG.kr1$MCMC.H.val[,1]), c(.075,.925))
quantile(as.mcmc(MCMCG.kr.rand1$MCMC.H.val[,1]), c(.075,.925)) #0.85 overlap
quantile(as.mcmc(MCMCG.kr1$MCMC.H.val[,1]), c(.1,.9))
quantile(as.mcmc(MCMCG.kr.rand1$MCMC.H.val[,1]), c(.1,.9)) #0.8 overlap
quantile(as.mcmc(MCMCG.kr1$MCMC.H.val[,1]), c(.125,.875))
quantile(as.mcmc(MCMCG.kr.rand1$MCMC.H.val[,1]), c(.125,.875)) #0.75 overlap
quantile(as.mcmc(MCMCG.kr1$MCMC.H.val[,1]), c(.15,.85))
quantile(as.mcmc(MCMCG.kr.rand1$MCMC.H.val[,1]), c(.15,.85)) #0.7 overlap
quantile(as.mcmc(MCMCG.kr1$MCMC.H.val[,1]), c(.175,.825))
quantile(as.mcmc(MCMCG.kr.rand1$MCMC.H.val[,1]), c(.175,.825)) # Pmcmc = 0.65


#### 3. Alignement with gmax (Third row of diagonal elements for Figure 2) ####
# Critical values of random distribution
quantile(rand.vec.cor,c(.9,.95,0.99,.999))

# gmax x h1
abs(vec.corr(eigen(G.Dunn.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,1]))
abs(vec.corr(eigen(G.AG.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,1]))
abs(vec.corr(eigen(G.SOC.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,1]))
abs(vec.corr(eigen(G.LC.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,1]))

# gmax x h2
abs(vec.corr(eigen(G.Dunn.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,2]))
abs(vec.corr(eigen(G.AG.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,2]))
abs(vec.corr(eigen(G.SOC.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,2]))
abs(vec.corr(eigen(G.LC.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,2]))

# gmax x h3
abs(vec.corr(eigen(G.Dunn.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,3]))
abs(vec.corr(eigen(G.AG.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,3]))
abs(vec.corr(eigen(G.SOC.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,3]))
abs(vec.corr(eigen(G.LC.Cov)$vectors[,1],eigen(MCMCG.kr1$avH)$vectors[,3]))

