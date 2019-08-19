library(MCMCglmm);library(ggthemes);library(patchwork);library(cowplot);library(ggbeeswarm);library(tidyverse);library(plyr)

#### Figure 2: Size, Shape, Orientation of G and alignment of gmax  ####
# Critical values of random distribution
quantile(rand.vec.cor,c(.9,.95,0.99,.999))
#### 1. Shape of G matrices (variance explained by gmax, first diagonal elements of Figure 2) ####
eigen(G.Dunn.Cov)$values/sum(eigen(G.Dunn.Cov)$values)
eigen(G.AG.Cov)$values/sum(eigen(G.AG.Cov)$values)
eigen(G.SOC.Cov)$values/sum(eigen(G.SOC.Cov)$values)
eigen(G.LC.Cov)$values/sum(eigen(G.LC.Cov)$values)

#### 2. Size of G matrices (Total genetic variance: sum of eigenvalues, Second diagonal elements of Figure 2) ####
sum(eigen(G.Dunn.Cov)$values)
sum(eigen(G.AG.Cov)$values)
sum(eigen(G.SOC.Cov)$values)
sum(eigen(G.LC.Cov)$values)

#### 4. Orientation of G matrices  (Alignement of h with gmax, third row of diagonal elements for Figure 2) ####
#### 5. Vector correlations of gmax among populations (Figure 2, lower diagonal elements) ####
(AGxDunn.angle<-abs(vec.corr(eigen(G.AG.Cov)$vectors[,1],eigen(G.Dunn.Cov)$vectors[,1])))
(SOCxDunn.angle<-abs(vec.corr(eigen(G.SOC.Cov)$vectors[,1],eigen(G.Dunn.Cov)$vectors[,1])))
(DunnxLC.angle<-abs(vec.corr(eigen(G.Dunn.Cov)$vectors[,1],eigen(G.LC.Cov)$vectors[,1])))
(SOCxAG.angle<-abs(vec.corr(eigen(G.SOC.Cov)$vectors[,1],eigen(G.AG.Cov)$vectors[,1])))
(AGxLC.angle<-abs(vec.corr(eigen(G.AG.Cov)$vectors[,1],eigen(G.LC.Cov)$vectors[,1])))
(SOCxLC.angle<-abs(vec.corr(eigen(G.SOC.Cov)$vectors[,1],eigen(G.LC.Cov)$vectors[,1])))


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
# Critical values of random distribution
quantile(rand.vec.cor,c(.9,.95,0.99,.999))


#### Figure 3: Change in average genetic and phenotypic correlations over generations ####
#### Figure 3a  Observed  vs Expected rA ####
# Observed mean rA
X=posterior.cor(Gmat.MCMC.F1$VCV[,1:49])
F1.ra.bar=rowMeans(abs(X[,c(2,3,4,5,6,7,10,11,12,13,14,
                            18,19,20,21,26,27,28,34,35,42)]))
posterior.mode(as.mcmc(F1.ra.bar))
HPDinterval(as.mcmc(F1.ra.bar))

X=posterior.cor(Gmat.MCMC.F2$VCV[,1:49])
F2.ra.bar=rowMeans(abs(X[,c(2,3,4,5,6,7,10,11,12,13,14,
                            18,19,20,21,26,27,28,34,35,42)]))
posterior.mode(as.mcmc(F2.ra.bar))
HPDinterval(as.mcmc(F2.ra.bar))
F2.ra.bar.exp=.5*F1.ra.bar
posterior.mode(as.mcmc(F2.ra.bar.exp))
HPDinterval(as.mcmc(F2.ra.bar.exp))

F2.dist=data.frame(c(as.numeric(F2.ra.bar),as.numeric(F2.ra.bar.exp)))
names(F2.dist)="Val"
F2.dist$Cat=factor(c(rep("Observed",1000),
                     rep("Expected",1000)),
                   levels = c("Observed", "Expected"))
# preliminary plots
F2.dist.plot=ggplot(F2.dist, aes(y=Val, x=Cat, fill=Cat, color=Cat)) + geom_beeswarm() + scale_color_wsj() + theme_test() +
  theme(legend.position = "none")

data_summary <- function(x) {
  m <- mean(x)
  ymin <- as.numeric(HPDinterval(as.mcmc(x))[1])
  ymax <- as.numeric(HPDinterval(as.mcmc(x))[2])
  return(c(y=m,ymin=ymin,ymax=ymax))
}

F2.dist.plot + 
  stat_summary(fun.data=data_summary, size=2, color="black")



mean.ra=rbind(mean(F2.ra.bar),
              mean(F2.ra.bar.exp))

HPD.ra=rbind(HPDinterval(as.mcmc(F2.ra.bar)),
             HPDinterval(as.mcmc(F2.ra.bar.exp)))

ra.bar=data.frame(cbind(mean.ra,HPD.ra))
names(ra.bar)=c("mean","low.ci","up.ci")
ra.bar$Cat=factor(c("Observed","Expected"), 
                  levels = c("Observed", "Expected"))

pd <- position_dodge(0.2)
ra.barPlot=ggplot(ra.bar, aes(y=mean, x=Cat, color=Cat)) +
  geom_point(size=2,position=pd) + geom_line(position=pd) +
  geom_errorbar(aes(ymin=low.ci,ymax=up.ci,width=.2),size=.75, position=pd) +
  scale_color_wsj() + ylim(0.1,.54) +
  ylab(expression(Mean~italic("r")["A"]~" "~F[2])) + xlab("") + theme_test() +
  theme(legend.position="none")
ra.barPlot

#### Figure 3b  Observed  vs Expected rP ####
# Observed mean rP
X=posterior.cor(Pmat.MCMC.F0$VCV[,1:49])
F0.r.bar=rowMeans(abs(X[,c(2,3,4,5,6,7,10,11,12,13,14,
                           18,19,20,21,26,27,28,34,35,42)]))
posterior.mode(as.mcmc(F0.r.bar))
HPDinterval(as.mcmc(F0.r.bar))

X=posterior.cor(Pmat.MCMC.F1$VCV[,1:49])
F1.r.bar=rowMeans(abs(X[,c(2,3,4,5,6,7,10,11,12,13,14,
                           18,19,20,21,26,27,28,34,35,42)]))
posterior.mode(as.mcmc(F1.r.bar))
HPDinterval(as.mcmc(F1.r.bar))

X=posterior.cor(Pmat.MCMC.F2$VCV[,1:49])
F2.r.bar=rowMeans(abs(X[,c(2,3,4,5,6,7,10,11,12,13,14,
                           18,19,20,21,26,27,28,34,35,42)]))
posterior.mode(as.mcmc(F2.r.bar))
HPDinterval(as.mcmc(F2.r.bar))

# Expected rP, done on Gmat.MCMC.ABQ$VCV, run on full correlation matrix
# rP.F0 = 2 G.F1 + E
# rP.F1 = G.F1 + E
# rP.F2 = 1/2 G.F1 + E

rp.F0.list<- vector("list",1000)
for(i in 1:1000){
  VG=diag(matrix(Gmat.MCMC.F1$VCV[,1:49][i,], nrow=7))
  G <- 2*matrix(Gmat.MCMC.F1$VCV[,1:49][i,], nrow=7) 
  diag(G)=VG
  E<- matrix(Gmat.MCMC.AllPop$VCV[,50:98][i,], nrow=7)
  rp.F0.list[[i]] <- cov2cor(G + E)
}

rp.F1.list<- vector("list",1000)
for(i in 1:1000){
  VG=diag(matrix(Gmat.MCMC.F1$VCV[,1:49][i,], nrow=7))
  G <- matrix(Gmat.MCMC.F1$VCV[,1:49][i,], nrow=7) 
  diag(G)=VG
  E<- matrix(Gmat.MCMC.AllPop$VCV[,50:98][i,], nrow=7)
  rp.F1.list[[i]] <- cov2cor(G + E)
}

rp.F2.list<- vector("list",1000)
for(i in 1:1000){
  VG=diag(matrix(Gmat.MCMC.F1$VCV[,1:49][i,], nrow=7))
  G <- 1/2*matrix(Gmat.MCMC.F1$VCV[,1:49][i,], nrow=7) 
  diag(G)=VG
  E<- matrix(Gmat.MCMC.AllPop$VCV[,50:98][i,], nrow=7)
  rp.F2.list[[i]] <- cov2cor(G + E)
}


F0.r.bar.exp <- vector("numeric",1000)
for(i in 1:1000){F0.r.bar.exp[[i]] <- 
  mean(abs(rp.F0.list[[i]][lower.tri(rp.F0.list[[i]])]))
}
posterior.mode(as.mcmc(F0.r.bar.exp))
HPDinterval(as.mcmc(F0.r.bar.exp))

F1.r.bar.exp <- vector("numeric",1000)
for(i in 1:1000){F1.r.bar.exp[[i]] <- 
  mean(abs(rp.F1.list[[i]][lower.tri(rp.F1.list[[i]])]))
}
posterior.mode(as.mcmc(F1.r.bar.exp))
HPDinterval(as.mcmc(F1.r.bar.exp))

F2.r.bar.exp <- vector("numeric",1000)
for(i in 1:1000){F2.r.bar.exp[[i]] <- 
  mean(abs(rp.F2.list[[i]][lower.tri(rp.F2.list[[i]])]))
}
posterior.mode(as.mcmc(F2.r.bar.exp))
HPDinterval(as.mcmc(F2.r.bar.exp))


mean.rP=rbind(mean(F0.r.bar),
              mean(F1.r.bar),
              mean(F2.r.bar),
              mean(F0.r.bar.exp),
              mean(F1.r.bar.exp),
              mean(F2.r.bar.exp))

HPD.rP=rbind(HPDinterval(as.mcmc(F0.r.bar)),
             HPDinterval(as.mcmc(F1.r.bar)),
             HPDinterval(as.mcmc(F2.r.bar)),
             HPDinterval(as.mcmc(F0.r.bar.exp)),
             HPDinterval(as.mcmc(F1.r.bar.exp)),
             HPDinterval(as.mcmc(F2.r.bar.exp)))

rp.bar=data.frame(cbind(mean.rP,HPD.rP))
names(rp.bar)=c("mean","low.ci","up.ci")
rp.bar$Generation=factor(rep(c("F0","F1","F2"),2))
rp.bar$Cat=factor(c(rep("Observed",3),
                    rep("Expected",3)), 
                  levels = c("Observed", "Expected"))
rp.bar2=subset(rp.bar, Generation!="F1")

pd <- position_dodge(0.2)
rp.barPlot=ggplot(rp.bar2, aes(y=mean, x=Generation, group=Cat, color=Cat)) +
  geom_point(size=2,position=pd) + geom_line(position=pd) +
  geom_errorbar(aes(ymin=low.ci,ymax=up.ci,width=.2),size=.75, position=pd) +
  scale_color_wsj() + #ylim(.2,.54) +
  ylab(expression(Mean~italic("r")["P"])) + theme_test() +
  theme(legend.title=element_blank(),
        legend.position=c(.15, .9))
rp.barPlot

Fig3=ra.barPlot/rp.barPlot
Fig3



#### Figure 4: Projecting BLUPS on D ####
# Original model runs had the Socorro population (SOC) coded as ABQ (Albuquerque) in the original datafram
# The Socorro population is therefore set as the intercept for these models


# Extract model intercepts, which indicates the average phenotype as deviation from the global intercep (Socorro population)
summary(Gmat.MCMC.AllPop.blups)
Mu.MCMC=data.frame(colMeans(Gmat.MCMC.AllPop.blups$Sol))
names(Mu.MCMC)="Mean"
Mu.MCMC$coef=factor(rownames(Mu.MCMC))
Mu.MCMC=Mu.MCMC[1:28,]
Mu.MCMC$Pop=factor(c(rep("SOC",7),
                     rep("AG",7),
                     rep("DUN",7),
                     rep("LC",7)),
                   levels = c("SOC","AG","DUN","LC"))
Mu.MCMC$Trait=factor(rep(c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo"),4),
                     levels = c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist",
                                "AP.Lat.Mov","AP.Var.Velo"))
Mu.MCMC=Mu.MCMC[,c(1,3,4)]
Mu.MCMC=spread(Mu.MCMC, Trait, Mean)
Mu.MCMC=Mu.MCMC[,2:8]
rownames(Mu.MCMC)=c("SOC","AG","DUN","LC")

# Add the intercept value (SOC column) to next three populations to recover the population average phenotypes
Mu.MCMC[2,]=Mu.MCMC[2,]+Mu.MCMC[1,]
Mu.MCMC[3,]=Mu.MCMC[3,]+Mu.MCMC[1,]
Mu.MCMC[4,]=Mu.MCMC[4,]+Mu.MCMC[1,]

# Transform Mu.MCMC into the phenotypic divergence matrix D (covariance matrix in mean pehnotypes among populationss)
D.MCMC=cov(Mu.MCMC)
D.MCMC

# Eigen decomposition of D
eigen(D.MCMC)
eigen(cov2cor(D.MCMC))$values/sum(eigen(cov2cor(D.MCMC))$values)*100
#eigen(D.MCMC)$values/sum(eigen(D.MCMC)$values)*100

#### 1. Calculate population corrdinates on d1 and d2 ####
D_Pop = as.matrix(Mu.MCMC) %*% eigen(D.MCMC)$vectors[,1:2]
D_Pop

# Store BLUPS into dataframe
Data_BLUPS=data.frame(posterior.mode(Gmat.MCMC.AllPop.blups$Sol))
str(Data_BLUPS)
Data_BLUPS$animal=factor(rownames(Data_BLUPS))
names(Data_BLUPS)=c("BLUPS","animal")
Data_BLUPS=Data_BLUPS[-c(1:109),]  # remove rows 1-109 cntaining all fixed effects coefficients
#Data_BLUPS$animal=factor(as.integer(c(1:965),7))
# Code each for its corresponding behavioral trait type
Data_BLUPS$Trait=factor(c(rep("Latency",1036),
                          rep("OF.Dist",2073-1036),
                          rep("UZ",3110-2073),
                          rep("OF.Var.Velo",4147-3110),
                          rep("AP.Dist",5183-4147),
                          rep("AP.Lat.Mov",6221-5183),
                          rep("AP.Var.Velo",7258-6221)),
                        levels = c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist",
                                   "AP.Lat.Mov","AP.Var.Velo"))

# Recover only the individul unique ID code from the MCMCglmm BLUPS object
animal_split=str_split_fixed(Data_BLUPS$animal,"animal.", n=2) %>% data.frame()
names(animal_split)=c("Trait","animal")

animal_split$Trait=str_sub(animal_split$Trait,6, str_length(animal_split$Trait)-1) %>% as.factor()
animal_split$Trait=revalue(animal_split$Trait, c("Latency.cen"="Latency"))
animal_split$Trait=factor(animal_split$Trait, 
                          levels=c("Latency","OF.Dist","UZ",
                                   "OF.Var.Velo","AP.Dist","AP.Lat.Mov",
                                   "AP.Var.Velo"))
str(animal_split)
#merge(animal_split,Data_BLUPS)


Data_BLUPS$animal=animal_split$animal
Data_BLUPS$Trait=animal_split$Trait

# Dataframe containing the BLUPS of all individuals in the dataset for each behavioral trait
Data_BLUPS=spread(Data_BLUPS, Trait, BLUPS)
str(Data_BLUPS)

# Assign populations code to each individual
# Rename Dunn as DUN for [lotting legend
Data2=Data
Data2$Pop=revalue(Data2$Pop, c("Dunn"="DUN"))
Data2$Pop=factor(Data2$Pop,levels = c("DUN", "AG", "SOC", "LC"))
Data_BLUPS=merge(Data2[,c("Pop","animal")],Data_BLUPS)

# Project BLUPS on the leading eigenvectors of phenotypic divergence: d1 and d2
D_BLUPS=as.data.frame(as.matrix(Data_BLUPS[,c(3:9)]) %*% eigen(D.MCMC)$vectors[,1:2])
names(D_BLUPS)=c("d1","d2")

# Assign populations code to each individual
D_BLUPS=cbind(Data_BLUPS$Pop,D_BLUPS)
names(D_BLUPS)=c("Pop","d1","d2")

# Values are centered around 0, need to add population coordinates to center pbservations around population means along d1 and d2
D_Pop_df=as.data.frame(D_Pop)
names(D_Pop_df)=c("d1_popmean","d2_popmean")
D_Pop_df$Pop=factor(rownames(D_Pop_df))

D_BLUPS=merge(D_BLUPS,D_Pop_df)
D_BLUPS$d1=D_BLUPS$d1+D_BLUPS$d1_popmean
D_BLUPS$d2=D_BLUPS$d2+D_BLUPS$d2_popmean

#### 2. Plot Figure ####
#### Fig 4a ####
#D_Pop_df=as.data.frame(D_Pop)
#names(D_Pop_df)=c("d_1","d_2")
#D_Pop_df$Pop=factor(rownames(D_Pop_df))

#D_BLUPS=merge(D_BLUPS,D_Pop_df)
#D_BLUPS$d1=D_BLUPS$d1+D_BLUPS$d_1
#D_BLUPS$d2=D_BLUPS$d2+D_BLUPS$d_2

Mu.MCMC_df=as.data.frame(Mu.MCMC)
names(Mu.MCMC_df)=c("Latency_mu","OF.Dist_mu","UZ_mu","OF.Var.Velo_mu","AP.Dist_mu","AP.Lat.Mov_mu","AP.Var.Velo_mu")
Mu.MCMC_df$Pop=factor(rownames(Mu.MCMC_df))
Data_BLUPS_mu=merge(Data_BLUPS,Mu.MCMC_df)
Data_BLUPS_mu[,c(3:9)]=Data_BLUPS_mu[,c(3:9)]+Data_BLUPS_mu[,c(10:16)]

Data_BLUPS_mu$Pop=factor(Data_BLUPS_mu$Pop,levels=c("DUN","AG","SOC","LC"))


BLUPS.APxLat.Plot=ggplot(Data_BLUPS_mu, aes(x= AP.Dist^2, y=Latency^2, group=Pop, fill=Pop, color=Pop)) +
  stat_ellipse(geom = "polygon", alpha = .7, size=1, color ="black") + geom_point(size=3, alpha = .5) + 
  ylab("Latency (s)") + xlab("Antipredator response (cm)") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() + theme(legend.position="none",
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16))
BLUPS.APxLat.Plot = BLUPS.APxLat.Plot + theme(aspect.ratio=1)

BLUPS.OFxLat.Plot=ggplot(Data_BLUPS_mu, aes(x= OF.Dist^2, y=Latency^2, group=Pop, fill=Pop, color=Pop)) +
  stat_ellipse(geom = "polygon", alpha = .7, size=1, color ="black") + geom_point(size=3, alpha = .5) + 
  ylab("Latency (s)") + xlab("Open-field distance (cm)") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() + theme(legend.position="none",
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16))
BLUPS.OFxLat.Plot = BLUPS.OFxLat.Plot + theme(aspect.ratio=1)


BLUPS.OFxAP.Plot=ggplot(Data_BLUPS_mu, aes(x= OF.Dist^2, y=AP.Dist^2, group=Pop, fill=Pop, color=Pop)) +
  stat_ellipse(geom = "polygon", alpha = .7, size=1, color ="black") + geom_point(size=3, alpha = .5) +
  ylab("Antipredator response (cm)") + xlab("Open-field distance (cm)") + scale_color_wsj() + scale_fill_wsj() +
  theme_test() + theme(legend.position="none",
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16))
#BLUPS.OFxAP.Plot + xlim(0,615) + ylim(0,615) + theme(aspect.ratio=1)
BLUPS.OFxAP.Plot = BLUPS.OFxAP.Plot + xlim(0,615) + ylim(0,615) + theme(aspect.ratio=1)


#### Fig 4b ####
#D_BLUPS$Pop=factor(D_BLUPS$Pop,levels=c("DUN","AG","SOC","LC"))

Fig4b=ggplot(D_BLUPS, aes(x= d1, y=d2, group=Pop, fill=Pop, color=Pop)) +
  stat_ellipse(geom = "polygon", alpha = .7, size=1, color ="black") + geom_point(size=5, alpha = .5) + 
  ylab(expression(paste(Divergence~axis~2, " (", bold(d[2]), ", ", "31.4%)"))) + 
  xlab(expression(paste(Divergence~axis~1, " (", bold(d[1]), ", ", "58.2%)"))) + 
  scale_color_wsj() + scale_fill_wsj() +
  theme_test() + theme(legend.position=c(.1, .9), legend.title = element_blank(),
                       axis.title = element_text(size=20),
                       axis.text = element_text(size=16))
Fig4b


#### Fig 4 compiled ###
Fig4 = BLUPS.APxLat.Plot + BLUPS.OFxLat.Plot + BLUPS.OFxAP.Plot + plot_layout(nrow = 3) | Fig4b
Fig4 =Fig4 + plot_layout(nrow=1,ncol=2, widths = c(1,3), heights = c(1,3))


ggsave("Fig4.jpeg", Fig1, device = "jpeg", width = 12*1.32, height = 12, units = "in")


#### Figure S1: Population divergence in covariance structure ####
#### Figure S1a. Population coordinates on E1 and E2 ####
H.coord.plot=data.frame(cbind(rbind(matrix(MCMC.covtensor$av.G.coord[,1,],ncol=1),
                                    matrix(MCMC.covtensor$av.G.coord[,2,],ncol=1)),
                              rbind(HPD.tensor.coord[,,1],HPD.tensor.coord[,,2])))
names(H.coord.plot)=c("mean","low.ci","up.ci")
H.coord.plot$tensor=factor(c(rep("E1",4),rep("E2",4)))
H.coord.plot$pop=factor(rep(c("SOC","AG","DUN","LC"),2),levels=c("DUN","AG","SOC","LC"))
#levels=c("DUN","AG","SOC","LC")


pd <- position_dodge(0.5)
FigS2a=ggplot(H.coord.plot,
             aes(y=mean, x=pop, color=pop)) +
  geom_point(size=4,position=pd) + 
  geom_errorbar(aes(ymin=low.ci,ymax=up.ci,width=.2),size=1, position=pd) +
  scale_color_wsj() + facet_wrap(~tensor, scales = "free") + theme_test() +
  theme(legend.position = "none",
        axis.title = element_text(size=20),
        axis.text = element_text(size=16),
        strip.text.x = element_text(size=21, face="bold"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  ylab("Coordinates in the eigentensor") + xlab("")
FigS2a


#### Figure 1b. Populations genetic variance in the direction of e11, e21 and e22 ####
HPD.e11 <- HPDinterval(t(as.mcmc(e11.proj)),prob = 0.95)
HPD.e21 <- HPDinterval(t(as.mcmc(e21.proj)),prob = 0.95)
HPD.e22 <- HPDinterval(t(as.mcmc(e22.proj)),prob = 0.95)

mu.e11=rowMeans(e11.proj)
mu.e21=rowMeans(e21.proj)
mu.e22=rowMeans(e22.proj)

H.Va.plot=data.frame(cbind(rbind(matrix(rowMeans(e11.proj),ncol=1),
                                 matrix(rowMeans(e21.proj),ncol=1),
                                 matrix(rowMeans(e22.proj),ncol=1)),
                           rbind(HPDinterval(t(as.mcmc(e11.proj)),prob = 0.95),
                                 HPDinterval(t(as.mcmc(e21.proj)),prob = 0.95),
                                 HPDinterval(t(as.mcmc(e22.proj)),prob = 0.95))))
names(H.Va.plot)=c("mean","low.ci","up.ci")
H.Va.plot$vector=factor(c(rep("e11",4),
                          rep("e21",4),
                          rep("e22",4)))
H.Va.plot$pop=factor(rep(c("SOC","AG","DUN","LC"),3),levels=c("DUN","AG","SOC","LC"))

FigS2b=ggplot(H.Va.plot,
             aes(y=mean, x=pop, color=pop)) +
  geom_point(size=4,position=pd) + 
  geom_errorbar(aes(ymin=low.ci,ymax=up.ci,width=.2),size=1, position=pd) +
  scale_color_wsj() + facet_wrap(~vector, scales = "free") + theme_test() +
  theme(legend.position = "none",
        axis.title = element_text(size=20),
        axis.text = element_text(size=16),
        strip.text.x = element_text(size=21, face="bold"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  ylab("Additive genetic variance") + xlab("")
FigS2b


FigS2=FigS2a/FigS2b



#### Figure S2: Krzanowski shared subspaces compared with random expectation ####
HPD.H.val <- cbind(HPDinterval(as.mcmc(MCMCG.kr1$MCMC.H.val), prob=0.95), 
                   HPDinterval(as.mcmc(MCMCG.kr.rand1$MCMC.H.val), prob=0.95))
means.H.val <- cbind(colMeans(as.mcmc(MCMCG.kr1$MCMC.H.val)), 
                     colMeans(as.mcmc(MCMCG.kr.rand1$MCMC.H.val)))
round(means.H.val, 3)
round(HPD.H.val, 3)

H.val.plot=data.frame(rbind(cbind(means.H.val[,1],HPD.H.val[,1:2]),
                            cbind(means.H.val[,2],HPD.H.val[,3:4])))
names(H.val.plot)=c("mean","low.ci","up.ci")
H.val.plot$eigen=factor(c("h1","h2","h3","h4","h5","h6","h7",
                          "h1","h2","h3","h4","h5","h6","h7"))
H.val.plot$type=factor(c(rep("Observed",7),rep("Random Expectation",7)))


pd <- position_dodge(0.5)
FigS2=ggplot(subset(H.val.plot,
                    eigen=="h1"|eigen=="h2"|eigen=="h3"),
             aes(y=mean, x=eigen, color=type)) +
  geom_point(size=4,position=pd) + 
  geom_errorbar(aes(ymin=low.ci,ymax=up.ci,width=.2),size=1, position=pd) +
  scale_color_wsj() + 
  ylab(expression(Eigenvalues~of~bold(H))) + xlab(expression(Eigenvectors~of~bold(H))) + theme_test()

FigS2=FigS2 + theme(legend.position = c(.222,.1), legend.title = element_blank(),
                    axis.title = element_text(size=14),
                    axis.text = element_text(size=12))
FigS2

#### Figure S3: Krzanowski shared subspaces compared with random expectation ####
HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero]), prob=0.95), 
                    HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero]), prob=0.95))
means.eT.val <- cbind(colMeans(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero])), 
                      colMeans(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero])))
round(means.eT.val, 3)
round(HPD.eT.val, 3)

eT.val.plot=data.frame(rbind(cbind(means.eT.val[,1],HPD.eT.val[,1:2]),
                             cbind(means.eT.val[,2],HPD.eT.val[,3:4])))
names(eT.val.plot)=c("mean","low.ci","up.ci")
eT.val.plot$eigen=factor(c("E1","E2","E3","E1","E2","E3"))
eT.val.plot$type=factor(c(rep("Observed",3),rep("Random Expectation",3)))


pd <- position_dodge(0.5)
FigS3=ggplot(eT.val.plot, aes(y=mean, x=eigen, color=type)) +
  geom_point(size=4,position=pd) + 
  geom_errorbar(aes(ymin=low.ci,ymax=up.ci,width=.2),size=1, position=pd) +
  scale_color_wsj() + 
  ylab(expression(Eigenvalues~of~bold(S))) + xlab(expression(Eigentensors~of~bold(S))) + theme_test()
FigS3=FigS3 + theme(legend.position = c(.72,.9), legend.title = element_blank(),
                    axis.title = element_text(size=14),
                    axis.text = element_text(size=12))
FigS3





#### Table S1: Pedigree information ####
# Nb of Dams
Data %>% group_by(Pop,Generation) %>% distinct(Dam) %>% tally()

# Nb of Sires
Data %>% group_by(Pop,Generation) %>% distinct(Sire) %>% tally()

# Nb of individuals in pedigree
Data %>% group_by(Pop,Generation) %>% tally() %>% data.frame()

# Nb of phenotyped individuals
# exclude non-phenotyped individuals
Data_NoNA=Data[-c(947:965),]
Data_NoNA %>% group_by(Pop,Generation) %>% tally() %>% data.frame()


#### Table S2: Primary eigenvectors of D, H and E ####
#### 1. Eigenvectors of D ####
eigen(D.MCMC)
eigen(cov2cor(D.MCMC))$values/sum(eigen(cov2cor(D.MCMC))$values)*100

#### 2. Eigenvectors of H ####
eigen(MCMCG.kr1$avH)
eigen(MCMCG.kr1$avH)$values/sum(eigen(MCMCG.kr1$avH)$values)*100


#### 3. Eigenvectors of E ####
eigen(MCMCG.kr1$avH)
eigen(MCMCG.kr1$avH)$values/sum(eigen(MCMCG.kr1$avH)$values)*100




#### Table S3: Alignment of D, H and E ####
#### 1. Alignment of D and H ####
#### 2. Alignment of D and E ####
#### 1. Alignment of E and H ####

#### Table S4: Genetic covariance matrices ####
G.Dunn.Cov;G.AG.Cov;G.SOC.Cov;G.LC.Cov
# get credible intervals 
HPDinterval(Gmat.MCMC.Dunn$VCV[,1:49])
HPDinterval(Gmat.MCMC.AG$VCV[,1:49])
HPDinterval(Gmat.MCMC.SOC$VCV[,1:49])
HPDinterval(Gmat.MCMC.LC$VCV[,1:49])


#### Table S5 Eigen decomposition of G by populations ####
eigen(G.Dunn.Cov)
eigen(G.AG.Cov)
eigen(G.SOC.Cov)
eigen(G.LC.Cov)

eigen(G.Dunn.Cov)$values/sum(eigen(G.Dunn.Cov)$values)*100
eigen(G.AG.Cov)$values/sum(eigen(G.AG.Cov)$values)*100
eigen(G.SOC.Cov)$values/sum(eigen(G.SOC.Cov)$values)*100
eigen(G.LC.Cov)$values/sum(eigen(G.LC.Cov)$values)*100
