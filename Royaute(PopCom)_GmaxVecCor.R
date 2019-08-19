#### 0. Vector correlation function and random distribution for vector correlations ####
vec.corr<-function(z1=z1,z2=z2){
  (sum(z1*z2) / ( sqrt(sum(z1 * z1)) * sqrt(sum(z2 * z2)) ) )
}


n.vec=1e6
rand.vec.cor=vector("numeric",n.vec)
for (i in 1:n.vec) {
  vec1=runif(7,-1,1)
  vec2=runif(7,-1,1)
  rand.vec.cor[i]=abs(vec.angle(vec1,vec2))
}

# original object can be also be imported:
rand.vec.cor=readRDS(file = "rand.vec.cor.rds")
hist(rand.vec.cor)
# Critical values
quantile(rand.vec.cor,c(.9,.95,0.99,.999))



#### 1. Eigen decomposition of G ####
#### 1.1 Get posterior means of G matrices ####
G.Dunn.Cov=matrix(colMeans(Gmat.MCMC.Dunn$VCV[,1:49]),nrow=7,
                  dimnames = list(c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo"),
                                  c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo")))
G.AG.Cov=matrix(colMeans(Gmat.MCMC.AG$VCV[,1:49]),nrow=7,
                dimnames = list(c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo"),
                                c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo")))
G.SOC.Cov=matrix(colMeans(Gmat.MCMC.SOC$VCV[,1:49]),nrow=7,
                 dimnames = list(c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo"),
                                 c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo")))
G.LC.Cov=matrix(colMeans(Gmat.MCMC.LC$VCV[,1:49]),nrow=7,
                dimnames = list(c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo"),
                                c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo")))

#### 1.2 Eigenvectors and eigenvalues of G by populations (Table S5) ####
eigen(G.Dunn.Cov)
eigen(G.AG.Cov)
eigen(G.SOC.Cov)
eigen(G.LC.Cov)

#### 1.3 Shape of G matrices (variance explained by gmax, first diagonal elements of Figure 2) ####
eigen(G.Dunn.Cov)$values/sum(eigen(G.Dunn.Cov)$values)
eigen(G.AG.Cov)$values/sum(eigen(G.AG.Cov)$values)
eigen(G.SOC.Cov)$values/sum(eigen(G.SOC.Cov)$values)
eigen(G.LC.Cov)$values/sum(eigen(G.LC.Cov)$values)

#### 1.4 Calcluate Size of G matrices (Total genetic variance: sum of eigenvalues, Second diagonal elements of Figure 2) ####
sum(eigen(G.Dunn.Cov)$values)
sum(eigen(G.AG.Cov)$values)
sum(eigen(G.SOC.Cov)$values)
sum(eigen(G.LC.Cov)$values)

#### 2. Calculate vector correlations of gmax among populations (Figure 2, lower diagonal elements) ####
(AGxDunn.angle<-abs(vec.corr(eigen(G.AG.Cov)$vectors[,1],eigen(G.Dunn.Cov)$vectors[,1])))
(SOCxDunn.angle<-abs(vec.corr(eigen(G.SOC.Cov)$vectors[,1],eigen(G.Dunn.Cov)$vectors[,1])))
(DunnxLC.angle<-abs(vec.corr(eigen(G.Dunn.Cov)$vectors[,1],eigen(G.LC.Cov)$vectors[,1])))
(SOCxAG.angle<-abs(vec.corr(eigen(G.SOC.Cov)$vectors[,1],eigen(G.AG.Cov)$vectors[,1])))
(AGxLC.angle<-abs(vec.corr(eigen(G.AG.Cov)$vectors[,1],eigen(G.LC.Cov)$vectors[,1])))
(SOCxLC.angle<-abs(vec.corr(eigen(G.SOC.Cov)$vectors[,1],eigen(G.LC.Cov)$vectors[,1])))

# Critical values of random distribution
quantile(rand.vec.cor,c(.9,.95,0.99,.999))

