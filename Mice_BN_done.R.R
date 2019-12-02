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

daten2 <- daten[sample(nrow(daten)),] #shuffled rows

#splitting to test and train
sample <- sample.int(n = nrow(daten), size = floor(.5*nrow(daten)), replace = F)
train <- daten[sample, ]
test  <- daten[-sample, ]
attach(daten)
#daten <- select(daten,-c 73)
test3<-daten #just another copy for myself


#######################
#PCA
#########################

dataPCA2 <- prcomp(test3[, 1:72])
biplot(dataPCA2)
pcaSummary <- summary(dataPCA2)$importance %>% as.data.frame()
pcaSummary$names <- rownames(pcaSummary)
pcaSummary %>%
  reshape::melt() %>%
  filter(names != "Standard deviation") %>%
  ggplot(aes(variable, value, color = names, group = names)) +
  geom_point() +
  geom_line() +
  ylab("% of variance") +
  xlab("principal component")

data.frame(dataPCA2$x, class = test3$class) %>%
  ggplot(., aes(PC1, PC2, color = class)) +
  geom_point()
data.frame(dataPCA2$x, class = test3$class) %>%
  ggplot(., aes(PC2, PC3, color = class)) +
  geom_point()
pca3d(dataPCA2,group=class, legend = "topleft")

######################################################

###KMEANS##

#######################################################


##############################

#SEARCH FOR K
searchKres <- foreach( i = 1:15, .combine = rbind) %do% {
  res <- kmeans(test3[, 1:4], centers = i, iter.max = 2000)
  data.frame(k = i, tot.withinss = res$tot.withinss)
}
plot(searchKres, type = "b")

# try k=4 because of praphics
res <- kmeans(test3[, 1:72], centers = 4, iter.max = 100000)
table(res$cluster, test3$class)
kmeanstable <- prop.table(table(res$cluster, test3$class),1)
kmeanstable
mean(diag(kmeanstable))

# try k=8 because we have 8 groups
res2 <- kmeans(test3[, 1:72], centers = 8, iter.max = 100000)
table(res2$cluster, test3$class)
kmeanstable2 <- prop.table(table(res2$cluster, test3$class),1)
kmeanstable2
mean(diag(kmeanstable2))
# 8 groups is better but still very bad accuracy

#use PCA before?
micekmeanspca2 <- kmeans(dataPCA2$x[, 1:5], 8, 10000)
table(micekmeanspca2$cluster, daten$class)
kmeanstable <- prop.table(table(micekmeanspca2$cluster, test3$class),1)
kmeanstable
mean(diag(kmeanstable))
# very similar... low accuracy

#####################################################

#Random Forest

#####################################################

for(i in 1:5){
    sample <- sample.int(n = nrow(daten), size = floor(.5*nrow(daten)), replace = F)
    train <- daten[sample, ]
    test  <- daten[-sample, ]
    df.rf <- randomForest(class ~ DYRK1A_N + ITSN1_N + BDNF_N+ NR1_N + NR2A_N + 
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
                        SHH_N + pS6_N + pCFOS_N + SYP_N  + 
                        CaNA_N, data = train,
                      importance = TRUE)
    varImpPlot(df.rf) #importance of proteins 
    test$pred.class.rf <- predict(df.rf, test)
    table(test$class, test$pred.class.rf)
    Rftable <- prop.table(table(test$class, test$pred.class.rf),1)
    Rftable
    print(mean(diag(Rftable)))
}
#Based on the mean decrease in accuracy plot, SOD1_N, pPKCG_N, Tau_N, APP_N, pNUMB_N, and pCAMKII_N are the most important proteins for classification

######################################################
#LDA
########################################################

for(i in 1:5){
  sample <- sample.int(n = nrow(daten), size = floor(.5*nrow(daten)), replace = F)
  train <- daten[sample, ]
  test  <- daten[-sample, ]
df.lda <- lda(class ~ DYRK1A_N + ITSN1_N + BDNF_N+ NR1_N + NR2A_N + 
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
                SHH_N + pS6_N + pCFOS_N + SYP_N +  
                CaNA_N, data = train)
test$pred.class.lda <- predict(df.lda, test)$class
table(test$class, test$pred.class.lda)
LDAtable <- prop.table(table(test$class, test$pred.class.lda),1)
LDAtable
print(mean(diag(LDAtable)))
}


##############IMPORTANCE OF PROTEINS
# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(class~., data=daten, method="lvq", preProcess="scale", trControl=control)
importance <- varImp(model, scale=FALSE)
print(importance)
plot(importance)
#########GOT SAME SOD1_N  pERK_N  CaNA_N Ubiquitin_N  DYRK1A_N     

#######################
#NAIVE BAYES
######################


for(i in 1:5){
  sample <- sample.int(n = nrow(daten), size = floor(.5*nrow(daten)), replace = F)
  train <- daten[sample, ]
  test  <- daten[-sample, ]
  modelNBClass <- naiveBayes(formula(class ~ .), data = train)
  predictionMBClass <- predict(modelNBClass, test)
  print(paste0("foldID ", i))
  print(caret::confusionMatrix(predictionMBClass, test$class))
  data.frame(
      foldID = i, 
      validationAccuracy = caret::confusionMatrix(
        predictionMBClass, test$class)$overall[1])
  }
  

##############################
#KNN
#############################

cvKNN <- function(test3, kcv = 5, kNN = 1) {
  set.seed(1) 
  idx <- sample(rep_len(1:kcv, nrow(test3)))
  res <- foreach (i = 1:kcv, .combine = rbind) %do% {
    dTrain <- test3[idx != i, which(names(test3) != "class")]
    dTest <- test3[idx == i, which(names(test3) != "class")]
    vTrain <- test3[idx != i, which(names(test3) == "class"), drop = TRUE] 
    vTest <- test3[idx == i, which(names(test3) == "class"), drop = TRUE]
    predictionKNN <- knn(train = dTrain, test = dTest, cl = vTrain, k = kNN)
    # it's also good to look at confusion table
    # print(caret::confusionMatrix(predictionKNN, vTest)$overall[1])
    data.frame(
      foldID = i, 
      kNN = kNN, 
      validationAccuracy = caret::confusionMatrix(predictionKNN, vTest)$overall[1])
  }
  return(res)
}
resultCVKNN <- foreach(kNN = 1:10, .combine = rbind) %do% {
  print(kNN)
  cvKNN(test3, kcv = 5, kNN = kNN)
}
resultCVKNN %>%
  group_by(kNN) %>% 
  summarise(meanAcc = mean(validationAccuracy)) 

######CHOOSING KNN 1, accuracy 0.993

agg = aggregate(resultCVKNN,
                by = list(resultCVKNN$kNN),
                FUN = mean)
ggplot(aes(x = kNN, y = validationAccuracy), data = agg) +  
  geom_point() + geom_line()


#########################
#nMDS JUST FOR TESTING IT OUT
###########################

d <- dist(test3[, -73]) 
fit <- vegan::metaMDS(d, k = 2) 
stressplot(fit)
plot(fit)
plotdata <- data.frame(fit$points, class = test3$class)
ggplot(plotdata, aes(MDS1, MDS2, color = class)) +
  geom_point()

#thank you :)
