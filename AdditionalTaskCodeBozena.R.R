# load needed packages 
#install.packages("pacman")
library(pacman)
p_load("dplyr", "ggplot2","foreach", "ggplot2", "datasets", "MASS", 
       "e1071", "class", "C50", "caret", "rpart.plot", "klaR", "corrplot", "VIM")  # should install and load all the needed packages
p_load("vegan", "reshape", "foreach")
library(dplyr);library(ggplot2);library(foreach);library(datasets);library(MASS)    
library(e1071);library(class);library(C50);library(foreach);library(caret);library(xlsx);library(plyr);library(dplyr);   
library(mice);library(VIM); library(corrplot)
library(pca3d);library(caret);library("rpart.plot");library(klaR);library(randomForest);library(mlbench)

set.seed(23456775)

miceINITIAL <- read.csv("Data_Cortex_Nuclear.csv")
miceINITIAL <- miceINITIAL[,-1]
aggr(miceINITIAL,sortVars=F, combined=T,bars=F,numbers=T,prop=F,sortCombs=T)#Checking Na structure
#removing NAs mice
miceINITIAL <- miceINITIAL[- (which(rowSums(is.na(miceINITIAL))>10)),]
aggr(miceINITIAL,sortVars=F, combined=T,bars=F,numbers=T,prop=F,sortCombs=T)#Checking Na structure
summary(miceINITIAL) #checking proteins with NAs
###REMOVING PROTEINS WITH MORE THAN 15% NAs
which( colnames(miceINITIAL)=="EGR1_N" )
miceINITIAL <- miceINITIAL[,-75]
which( colnames(miceINITIAL)=="H3MeK4_N" )
miceINITIAL <- miceINITIAL[,-75]
which( colnames(miceINITIAL)=="H3AcK18_N" )
miceINITIAL <- miceINITIAL[,-74]
which( colnames(miceINITIAL)=="BCL2_N" )
miceINITIAL <- miceINITIAL[,-70]
which( colnames(miceINITIAL)=="BAD_N" )
miceINITIAL <- miceINITIAL[,-69]

aggr(miceINITIAL,sortVars=F, combined=T,bars=F,numbers=T,prop=F,sortCombs=T)#Checking Na structure

##FILLING NAs with mean
for(i in 1:ncol(miceINITIAL)){
  miceINITIAL[is.na(miceINITIAL[,i]), i] <- mean(miceINITIAL[,i], na.rm = TRUE)
}
aggr(miceINITIAL,sortVars=F, combined=T,bars=F,numbers=T,prop=F,sortCombs=T) #Checking Na structure if its filled
#all NAs FILLED with means

miceINITIAL <- miceINITIAL[,-73]
miceINITIAL <- miceINITIAL[,-73]
miceINITIAL <- miceINITIAL[,-73] #removing genotype, behavior and treatment

#normalized data

daten<-miceINITIAL #just copy
miceINITIAL <- as.data.frame(apply(miceINITIAL[, 1:72], 2, function(x) (x - min(x))/(max(x)-min(x))))
miceINITIAL$class <- daten$class

daten<-miceINITIAL #just copy



###1
cl<-as.numeric(daten$class)
daten$numclass<-cl
attach(daten)
mice_ANOVA<-aov(numclass ~ DYRK1A_N + ITSN1_N + BDNF_N+ NR1_N + NR2A_N + 
                   pAKT_N + pBRAF_N + pCAMKII_N + pCREB_N + pELK_N +
                   pERK_N + pJNK_N + PKCA_N + pMEK_N + pNR1_N + pNR2A_N +
                   pNR2B_N + pPKCAB_N + pRSK_N + AKT_N + BRAF_N +
                   CAMKII_N + CREB_N + ELK_N + ERK_N + GSK3B_N + JNK_N + 
                   MEK_N + TRKA_N + RSK_N + APP_N + Bcatenin_N + SOD1_N +
                   MTOR_N + P38_N + pMTOR_N + DSCR1_N + AMPKA_N + NR2B_N +
                   pNUMB_N + RAPTOR_N + TIAM1_N + pP70S6_N + NUMB_N + 
                   pP70S6_N + pGSK3B_N + pPKCG_N + CDK5_N + S6_N + ADARB1_N + 
                   AcetylH3K9_N + RRP1_N + BAX_N + ARC_N + ERBB4_N + nNOS_N + 
                   Tau_N + GFAP_N + GluR3_N + GluR4_N + IL1B_N + P3525_N + 
                   pCASP9_N + PSD95_N + SNCA_N + Ubiquitin_N + pGSK3B_Tyr216_N + 
                   SHH_N + pS6_N + pCFOS_N + SYP_N + CaNA_N)

summary(mice_ANOVA)


###2
attach(miceINITIAL)
b <- ggplot(miceINITIAL, aes(x = ITSN1_N, y = pERK_N))
b + geom_point()
c <- ggplot(miceINITIAL, aes(x = BRAF_N, y = pERK_N))
c + geom_point()
d <- ggplot(miceINITIAL, aes(x = APP_N, y = pERK_N))
d + geom_point()
sod <- ggplot(miceINITIAL, aes(x = SOD1_N, y = pERK_N))
sod + geom_point()
b + geom_point(aes(shape = class, color = class)) +  scale_shape_manual(values=seq(0,15)) +ylim(0, 0.5) + xlim(0, 0.5)
c + geom_point(aes(shape = class, color = class))  +  scale_shape_manual(values=seq(0,15)) +ylim(0, 0.4) + xlim(0, 0.4)
d + geom_point(aes(shape = class, color = class))  +  scale_shape_manual(values=seq(0,15)) +ylim(0, 0.5)
sod + geom_point(aes(shape = class, color = class))  +  scale_shape_manual(values=seq(0,15)) + ylim(0, 0.4)


cor(pERK_N, ITSN1_N)
cor(pERK_N, BRAF_N)
cor(pERK_N, APP_N)
cor(pERK_N, SOD1_N)


fit<-lm(pERK_N ~ ITSN1_N + BRAF_N + APP_N) 
summary(fit)
plot(fit)
