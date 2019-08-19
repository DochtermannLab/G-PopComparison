#rm(list=ls())

library(tidyverse); library(MCMCglmm);library(pedantics)


#### 0. Data import ####
#### 0.1 MCMCglmm models import ####
# G matrices
Gmat.MCMC.Dunn=readRDS(file = "Gmat.MCMC.Dunn.rds")
Gmat.MCMC.AG=readRDS(file = "Gmat.MCMC.AG.rds")
Gmat.MCMC.SOC=readRDS(file = "Gmat.MCMC.SOC.rds")
Gmat.MCMC.LC=readRDS(file = "Gmat.MCMC.LC.rds")
Gmat.MCMC.AllPop=readRDS(file = "Gmat.MCMC.AllPop.rds")
Gmat.MCMC.AllPop.blups=readRDS(file = "Gmat.MCMC.AllPop.blups.rds")

# P matrices per generation
Pmat.MCMC.F0=readRDS(file = "Pmat.MCMC.F0.rds")
Pmat.MCMC.F1=readRDS(file = "Pmat.MCMC.F1.rds")
Pmat.MCMC.F2=readRDS(file = "Pmat.MCMC.F2.rds")

#G matrices per generation
Gmat.MCMC.F1=readRDS(file = "Gmat.F1.MCMC.rds")
Gmat.MCMC.F2=readRDS(file = "Gmat.F2.MCMC.rds")


#### 0.2 Raw Data Import ####
#Data<- read.csv("G.integer.Aim1.Full.csv",header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
Data<- read.csv("Royaute(PopComp)_Data.csv",header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
str(Data)


Data_Dunn=subset(Data, Pop=="Dunn")
Data_AG=subset(Data, Pop=="AG")
Data_SOC=subset(Data, Pop=="SOC")
Data_LC=subset(Data, Pop=="LC")


#### 0.3 Pedigree generation ####
PooledPedigree2<-Data[,c(1,3,2)]
PooledPedigree.analysis<-fixPedigree(PooledPedigree2)
drawPedigree(PooledPedigree.analysis)
View(PooledPedigree.analysis)

# Pedigree by pop
PooledPedigree.analysis.SOC=fixPedigree(Data_SOC[,c(1,3,2)])
PooledPedigree.analysis.AG=fixPedigree(Data_AG[,c(1,3,2)])
PooledPedigree.analysis.Dunn=fixPedigree(Data_Dunn[,c(1,3,2)])
PooledPedigree.analysis.LC=fixPedigree(Data_LC[,c(1,3,2)])

#### 0.4 Summary stats for table S1 ####
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



#############################################################################################
#### The following sections (1-6) Show the MCMCglmm code necessary to run the animal models #
#### Because MCMCglmm takes a long time to run these data, these sections can be skipped    #
#### All additional R files use the original MCMCglmm object imported in sectin 0.1         #
#############################################################################################

#### 1. Run Models ####

#Define MCMC chain iterations 
Mult=1;NITT=120000;THIN=100;BURN=20000
#Mult=.5
Mult=40

#### 2. Gmat by Populations ####
#### 2.1 PRIOR BASED ON LL estimates ####
library(lme4)
summary(lmer(sqrt(Latency.Cens)~Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev+(1|Sire),Data,na.action=na.omit)) # Va=4*2.858, Vr=31.434
summary(lmer(sqrt(OF.Dist)~Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev+(1|Sire),Data,na.action=na.omit)) # Va=4*3.912, Vr=94.845
summary(glmer(UZ~Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev+(1|Sire),Data,na.action=na.omit,family = "poisson")) # Va=4*0.1912, Vr=1
summary(lmer(log(OF.Var.Velo+1)~Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev+(1|Sire),Data,na.action=na.omit)) # Va=4*0.04876, Vr=1.47720
summary(lmer(sqrt(AP.Dist)~Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev+(1|Sire),Data,na.action=na.omit)) # Va=4*5.483, Vr=83.950
summary(glmer(round(AP.Lat.Mov,0)~Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev+(1|Sire),Data,na.action=na.omit, family = "poisson")) # Va=4*0.4556, Vr=1
summary(lmer(log(AP.Var.Velo+1)~Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev+(1|Sire),Data,na.action=na.omit)) # Va=4*0.03946, Vr=0.69877

Va=c(4*2.858,4*3.912,4*0.1912,4*0.04876,4*5.483,4*0.4556,4*0.03946)
Vr=c(31.434,94.845,1,1.47720,83.950,1,0.69877)

prior.7v<-
  list(R=list(V=diag(7)*Vr,nu=7),
       G=list(G1=list(V=diag(7)*Va,nu=7)))

#### 2.2 run Models ####
# Note: This portion of the script takes about 1 week to run with the full number of iteration.
# We recommend lowering the MULT argument for a quicker inspection of the models (MULT<10)
# The original models are are also provided in the Royaute(PopComp)_ModelImport.R script
Gmat.MCMC.Dunn<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           random = ~us(trait):animal,rcov=~us(trait):units,
           pedigree = PooledPedigree.analysis.Dunn,data=Data_Dunn, 
           prior=prior.7v, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult)

Gmat.MCMC.AG<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           random = ~us(trait):animal,rcov=~us(trait):units,
           pedigree = PooledPedigree.analysis.AG,data=Data_AG, 
           prior=prior.7v, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult)

Gmat.MCMC.SOC<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           random = ~us(trait):animal,rcov=~us(trait):units,
           pedigree = PooledPedigree.analysis.SOC,data=Data_SOC, 
           prior=prior.7v, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult)

Gmat.MCMC.LC<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           random = ~us(trait):animal,rcov=~us(trait):units,
           pedigree = PooledPedigree.analysis.LC,data=Data_LC, 
           prior=prior.7v, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult)#,



#### 3. Gmat All Populations ####
Gmat.MCMC<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           random = ~us(trait):animal,rcov=~us(trait):units,
           pedigree = PooledPedigree.analysis,data=Data, 
           prior=prior.7v, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult)


# The following version was used to extract the BLUPs for plotting Figure 4:
Gmat.MCMC.blups<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           random = ~us(trait):animal,rcov=~us(trait):units,
           pedigree = PooledPedigree.analysis,data=Data, 
           prior=prior.7v, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult, pr=T)



#### 4. Phenotypic correlations by generation ####
prior.7v.P<-list(R=list(V=diag(7),nu=7))

Pmat.MCMC.F0<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           rcov=~us(trait):units,
           data=subset(Data, Generation=="F0"), 
           prior=prior.7v.P, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult, pr = T)
summary(Pmat.MCMC.F0)
Pmat.F0<-cov2cor(matrix(posterior.mode(posterior.cor(Pmat.MCMC.F0$VCV[,1:49])),nrow=7,
                        dimnames = list(c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo"),
                                        c("Latency","OF.Dist","UZ","OF.Var.Velo","AP.Dist","AP.Lat.Mov","AP.Var.Velo"))))
Pmat.F0
HPDinterval(posterior.cor(Pmat.MCMC.F0$VCV[,1:49]))


Pmat.MCMC.F1<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           rcov=~us(trait):units,
           data=subset(Data, Generation=="F1"), 
           prior=prior.7v.P, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult, pr = T)


Pmat.MCMC.F2<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           rcov=~us(trait):units,
           data=subset(Data, Generation=="F2"), 
           prior=prior.7v.P, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult, pr = T)


#### 5. Genetic correlations by generation ####
# Gmat F0-F1
Gmat.F1.MCMC<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           random = ~us(trait):animal,rcov=~us(trait):units,
           pedigree = PooledPedigree.analysis,data=subset(Data, Generation!="F2"), 
           prior=prior.7v, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult)

Gmat.F2.MCMC<-
  MCMCglmm(cbind(Latency.cen,Latency.up,
                 sqrt(OF.Dist),UZ,log(OF.Var.Velo+1),
                 sqrt(AP.Dist),as.integer(AP.Lat.Mov),log(AP.Var.Velo+1))~(trait-1)+
             trait:(Pop+Mass+Generation+Sex+Stage+Injured+Temp+Time+JDate_dev)+
             at.level(trait,5):(Understimate.AP)+
             at.level(trait,6):(Understimate.AP)+
             at.level(trait,7):(Understimate.AP),
           random = ~us(trait):animal,rcov=~us(trait):units,
           pedigree = PooledPedigree.analysis,data=subset(Data, Generation!="F0"), 
           prior=prior.7v, verbose=F,
           family=c("cengaussian","gaussian","poisson","gaussian",
                    "gaussian","poisson","gaussian"),
           nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult)

#### 6. Extract Analyses files ####
save.image(file = "Royaute(PopComp)_Image.RData")

# Extract each model
#Gmat
saveRDS(Gmat.MCMC.SOC, file = "Gmat.MCMC.SOC.rds")
saveRDS(Gmat.MCMC.AG, file = "Gmat.MCMC.AG.rds")
saveRDS(Gmat.MCMC.Dunn, file = "Gmat.MCMC.Dunn.rds")
saveRDS(Gmat.MCMC.LC, file = "Gmat.MCMC.LC.rds")
saveRDS(Gmat.MCMC, file = "Gmat.MCMC.AllPop.rds")
saveRDS(Gmat.MCMC.blups, file = "Gmat.MCMC.AllPop.blups.rds")

#Pmat per generation
saveRDS(Pmat.MCMC.F0, file = "Pmat.MCMC.F0.rds")
saveRDS(Pmat.MCMC.F1, file = "Pmat.MCMC.F1.rds")
saveRDS(Pmat.MCMC.F2, file = "Pmat.MCMC.F2.rds")
#Gmat per generation
saveRDS(Gmat.F1.MCMC, file = "Gmat.F1.MCMC.rds")
saveRDS(Gmat.F2.MCMC, file = "Gmat.F2.MCMC.rds")

