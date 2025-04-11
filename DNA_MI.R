# rm (list = ls( )) 
options(expressions = 500000)
setwd("/Users/rifatzahan/Desktop/Working Files/CMPT 830/Project/GSE88890_RAW/")
update.packages(ask = FALSE, dependencies = c ( 'Suggests' ))


library(foreign)
library(varhandle)
library(InformationValue)
library(gplots)
library(ROCR)    
library(Rtsne)
library(e1071)
library(rpart)
library(survival)
library(epiR)
library(lattice)
library(ggplot2)
library(caret)
library(kernlab)
library(tsne)
library(magrittr)
library(ggpubr)
library(dplyr)
library(pca3d)
library(rgl)
library(Matrix)
library(foreach)
library(glmnet)
library(e1071)
library(rgl)
library(misc3d)


id_dat <- read.csv("Individuals_Description.csv", header = T)
genetic_dat <- read.csv("GSE88890_normalisedBetas.csv", header = T)
genetic_dat <- sample_n(genetic_dat, 15000)
genetic_dat <- as.data.frame(t(genetic_dat))
colnames(genetic_dat) <- as.character(unlist(genetic_dat[1,]))
genetic_dat <- genetic_dat[-1, ]
genetic_dat$ID <- factor(row.names(genetic_dat))

total <- merge(id_dat,genetic_dat,by="ID")

ba11_dat <- subset(total, Cortex == "BA11")


##########
## BA11 ##
##########
ba11_dat_genes <- ba11_dat
ba11_dat_genes[1:6] <- list(NULL)
ba11_dat_genes <- unfactor(ba11_dat_genes)

ba11_dat_genes_group <- cbind(ba11_dat$Group, ba11_dat_genes)
colnames(ba11_dat_genes_group)[1] <- "Group"
deathClasses <- factor(ba11_dat_genes_group$Group)


ba11_dat_genes_group$group_coded[ba11_dat_genes_group$Group == "MDD Suicide"] <- 1 
ba11_dat_genes_group$group_coded[ba11_dat_genes_group$Group == "NPS Death"] <- 0 
ba11_dat_genes_group$group_coded

####################################
### Principal Component Analysis ###
####################################
# perform PCA
# https://www.r-bloggers.com/principal-component-analysis-in-r/


BA11PCA <- prcomp(scale(ba11_dat_genes_group[,-1]))

ggplot(as.data.frame(BA11PCA$x[,1:2]))+ theme_bw() + 
  geom_point(aes(x=PC1,y=PC2, colour = ba11_dat_genes_group$Group))

scoresGenes <- BA11PCA$x

std_dev<- BA11PCA$sdev
BA11PCA.var<- std_dev^2
round(BA11PCA.var)

#proportion of variance explained
prop_varex <- BA11PCA.var/sum(BA11PCA.var)
cumsum(round(prop_varex,3)*100)


plot(cumsum(prop_varex), xlab = "Principal Component",ylab = "Cumulative Proportion of Variance Explained",type = "b")


#############
### T-SNE ###
#############

trainBA11<- ba11_dat_genes_group 
Death <- trainBA11$Group
trainBA11$label <- Death
colorsBA11 = rainbow(length(unique(trainBA11$label)))
names(colorsBA11) = unique(trainBA11$label)

tsne_BA11_1 <- Rtsne(trainBA11[,-1], dims = 1, perplexity=5, epoch = 5, costs = 20, pca = FALSE) 
tsne_BA11_2 <- Rtsne(trainBA11[,-1], dims = 2, perplexity=5, epoch = 5, costs = 20, pca = FALSE) 
tsne_BA11_3 <- Rtsne(trainBA11[,-1], dims = 3, perplexity=5, epoch = 5, costs = 20, pca = FALSE) 
tsne_BA11_4 <- Rtsne(trainBA11[,-1], dims = 4, perplexity=5, epoch = 5, costs = 20, pca = FALSE) 


tsneplot <- ggplot(as.data.frame(tsne_BA11_2$Y[,1:2]))+ theme_bw() + 
  geom_point(aes(x=tsne_BA11_2$Y[,1],y=tsne_BA11_2$Y[,2], colour = Death)) +
  xlab("t-SNE1") +
  ylab("t-SNE2") + 
  scale_fill_discrete(name="Death", breaks=c("MDD Suicide", "NPS Death")) + 
  ggtitle("t-SNE") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

pcaplot <- ggplot(as.data.frame(BA11PCA$x[,1:2]))+ theme_bw() + 
  geom_point(aes(x=PC1,y=PC2, colour = Death)) + 
  xlab("PC1") +
  ylab("PC2") + 
  ggtitle("PCA") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggarrange(pcaplot, tsneplot,
          ncol = 1, nrow = 2)



plot3d(x=tsne_BA11_3$Y[,1],y=tsne_BA11_3$Y[,2],z=tsne_BA11_3$Y[,3],
       xlab = "t-SNE1", ylab = "t-SNE2", zlab = "t-SNE3",
       col=c("red","green")[trainBA11[,1]],
       type="s",radius=2.5)
legend3d("topright", legend = paste(c('MDD Suicide', 'NPS Death')), 
         pch = 16, col = c("red","green"), cex=2, inset=c(0.02))
snapshot3d(filename = 'tsne3dba11.png', fmt = 'png')

plot3d(x=BA11PCA$x[,1],y=BA11PCA$x[,2],z=BA11PCA$x[,3],
       xlab = "PC1", ylab = "PC2", zlab = "PC3",
       col=c("red","green")[trainBA11[,1]],
       type="s",radius=6)
legend3d("topright", legend = paste(c('MDD Suicide', 'NPS Death')), 
         pch = 16, col = c("red","green"), cex=2, inset=c(0.02))
snapshot3d(filename = 'pca3d11.png', fmt = 'png')




###########################
### Logistic regression ###
###########################

### With PCA ###

modGenes1 <- glm(ba11_dat_genes_group$group_coded ~ scoresGenes[,1], family = binomial(link = "logit"))
modGenes2 <- glm(ba11_dat_genes_group$group_coded ~ scoresGenes[,1:2], family = binomial(link = "logit"))
modGenes3 <- glm(ba11_dat_genes_group$group_coded ~ scoresGenes[,1:3], family = binomial(link = "logit"))

# http://r-statistics.co/Logistic-Regression-With-R.html
predicted1 = predict(modGenes1,ba11_dat_genes_group$Group, type='response')
predicted2 = predict(modGenes2,ba11_dat_genes_group$Group, type='response')
predicted3 = predict(modGenes3,ba11_dat_genes_group$Group, type='response')

optCutOff1 <- optimalCutoff(ba11_dat_genes_group$group_coded, predicted1)[1] 
optCutOff2 <- optimalCutoff(ba11_dat_genes_group$group_coded, predicted2)[1] 
optCutOff3 <- optimalCutoff(ba11_dat_genes_group$group_coded, predicted3)[1]


misClassError(ba11_dat_genes_group$group_coded, predicted1, threshold = optCutOff1)
misClassError(ba11_dat_genes_group$group_coded, predicted2, threshold = optCutOff2)
misClassError(ba11_dat_genes_group$group_coded, predicted3, threshold = optCutOff3)




pred1 <- prediction(predicted1, ba11_dat_genes_group$group_coded)    
pred2 <- prediction(predicted2, ba11_dat_genes_group$group_coded)    
pred3 <- prediction(predicted3, ba11_dat_genes_group$group_coded) 


perf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")    
perf2 <- performance(pred2, measure = "tpr", x.measure = "fpr")     
perf3 <- performance(pred3, measure = "tpr", x.measure = "fpr")  

plot(perf1, col=1, main="ROC curves for Logistic Regression Model", xlab="False Positive Rate", 
     ylab="True Positive Rate")    
plot(perf2, add = TRUE, col = 2)
plot(perf3, add = TRUE, col = 3)
legend(0.7, 0.3, c('PC1', 'PC2', 'PC3'), 1:3)
abline(0, 1)



### Logistic With t-SNE ###


# dim = 2 
ba11_dat_genes_group$tsne1 <- tsne_BA11_2$Y[,1]
ba11_dat_genes_group$tsne2 <- tsne_BA11_2$Y[,2]

ba11regdata <- subset(ba11_dat_genes_group, select = c(Group, group_coded, tsne1, tsne2))
indexes = sample(1:nrow(ba11regdata), size=0.25*nrow(ba11regdata))
test2 =  ba11regdata[indexes,]
testtsne2 <- subset(test2, select = -c(Group, group_coded))
dim(testtsne2)  # 10 2
traindat2 = ba11regdata[-indexes,]
dim(traindat2) # 30 4


# https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html


modGenestsne2 <- cv.glmnet(x = as.matrix(traindat2[,-(1:2)]), y = traindat2[,1], 
                           family = 'binomial', type.measure = 'auc')
testtsne2$lasso.prob <- predict(modGenestsne2,type="response", 
                               newx = as.matrix(testtsne2), s = 'lambda.min')


# dim = 3 
ba11_dat_genes_group$tsne1 <- tsne_BA11_3$Y[,1]
ba11_dat_genes_group$tsne2 <- tsne_BA11_3$Y[,2]
ba11_dat_genes_group$tsne3 <- tsne_BA11_3$Y[,3]

ba11regdata <- subset(ba11_dat_genes_group, select = c(Group, group_coded, tsne1, tsne2, tsne3))
indexes = sample(1:nrow(ba11regdata), size=0.25*nrow(ba11regdata))
test3 =  ba11regdata[indexes,]
testtsne3 <- subset(test3, select = -c(Group, group_coded))
dim(testtsne3)  # 10 2
traindat3 = ba11regdata[-indexes,]
dim(traindat3) # 30 4


modGenestsne3 <- cv.glmnet(x = as.matrix(traindat3[,-(1:2)]), y = traindat2[,1], 
                           family = 'binomial', type.measure = 'auc')
testtsne3$lasso.prob <- predict(modGenestsne3,type="response", 
                                newx = as.matrix(testtsne3), s = 'lambda.min')


pred2 <- prediction(testtsne2$lasso.prob, test2$group_coded)
pred3 <- prediction(testtsne3$lasso.prob, test3$group_coded)


perf2 <- performance(pred2,"tpr","fpr")
perf3 <- performance(pred3,"tpr","fpr")



plot(perf2, col=1, main="ROC curves for Logistic Regression Model", xlab="False Positive Rate", 
     ylab="True Positive Rate")
plot(perf3, col = 2, add = TRUE)
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
legend(0.7, 0.3, c('Lambda = 0.1', 'Lamba = '), 1:5)


#####################################
## Support Vector Machine with PCA ##
#####################################

testSet <- cbind(ba11_dat_genes_group$group_coded, scoresGenes[,1:3])

svm.model <- svm(Group ~ scoresGenes[,1], data = ba11_dat_genes_group, cost = 1, gamma = 2)
svm.pred  <- predict(svm.model, scoresGenes[,1])
tableMat <- table(pred = svm.pred, true = ba11_dat_genes_group$Group)
epi.tests(tableMat)



svm.model <- svm(Group ~ scoresGenes[,1], data = ba11_dat_genes_group, cost = 1, gamma = 2,  probability=TRUE)
svm.pred  <- predict(svm.model, scoresGenes[,1],probability=TRUE)

svm.model1 <- svm(Group ~ scoresGenes[,1:2], data = ba11_dat_genes_group, cost = 1, gamma = 2,  probability=TRUE)
svm.pred1  <- predict(svm.model1, scoresGenes[,1:2],probability=TRUE)

svm.model2 <- svm(Group ~ scoresGenes[,1:3], data = ba11_dat_genes_group, cost = 1, gamma = 2,  probability=TRUE)
svm.pred2  <- predict(svm.model2, scoresGenes[,1:3],probability=TRUE)


svmmod.rocr<-prediction(attr(svm.pred,"probabilities")[,2], ba11_dat_genes_group$Group)
svmmod1.rocr<-prediction(attr(svm.pred1,"probabilities")[,2], ba11_dat_genes_group$Group)
svmmod2.rocr<-prediction(attr(svm.pred2,"probabilities")[,2], ba11_dat_genes_group$Group)

svm.perf<-performance(svmmod.rocr, measure = "tpr", x.measure = "fpr")
svmmod1.perf<-performance(svmmod1.rocr, measure = "tpr", x.measure = "fpr")
svmmod2.perf<-performance(svmmod2.rocr, measure = "tpr", x.measure = "fpr")

plot(svm.perf,col=1)
plot(svmmod1.perf,add=TRUE,col=2)
plot(svmmod2.perf,add=TRUE,col=3)
legend(0.7, 0.3, c('PC1', 'PC2', 'PC3'), 1:3)
abline(0, 1)

svm.auc<-as.numeric(performance(svmmod.rocr, measure = "auc", x.measure = "cutoff")@ y.values)
svmmod1.auc<-as.numeric(performance(svmmod1.rocr, measure = "auc", x.measure = "cutoff")@ y.values)
svmmod2.auc<-as.numeric(performance(svmmod2.rocr, measure = "auc", x.measure = "cutoff")@ y.values)


svm.auc
svmmod1.auc
svmmod2.auc

#######################################
## Support Vector Machine with t-SNE ##
#######################################

# dim = 2
ba11_dat_genes_group$tsne1 <- tsne_BA11_2$Y[,1]
ba11_dat_genes_group$tsne2 <- tsne_BA11_2$Y[,2]
ba11regdata2 <- subset(ba11_dat_genes_group, select = c(Group, group_coded, tsne1, tsne2))
indexes = sample(1:nrow(ba11regdata2), size=0.25*nrow(ba11regdata2))
test2 =  ba11regdata2[indexes,]
testtsne2 <- subset(test2, select = -c(Group, group_coded))
dim(testtsne2)  # 10 2
traindat2 = ba11regdata2[-indexes,]
dim(traindat2) # 30 4

# dim = 3
ba11_dat_genes_group$tsne1 <- tsne_BA11_3$Y[,1]
ba11_dat_genes_group$tsne2 <- tsne_BA11_3$Y[,2]
ba11_dat_genes_group$tsne3 <- tsne_BA11_3$Y[,3]
ba11regdata3 <- subset(ba11_dat_genes_group, select = c(Group, group_coded, tsne1, tsne2, tsne3))
indexes = sample(1:nrow(ba11regdata3), size=0.25*nrow(ba11regdata3))
test3 =  ba11regdata3[indexes,]
testtsne3 <- subset(test3, select = -c(Group, group_coded))
dim(testtsne3)  # 10 3
traindat3 = ba11regdata3[-indexes,]
dim(traindat3) # 30 5


svm_tune2 <- tune(svm, train.x=traindat2$tsne1 + traindat2$tsne2, train.y=traindat2$Group, 
                 kernel="polynomial", probability=TRUE, ranges=list(cost=10^(-1:3), gamma=c(.5,1,2)))
summary(svm_tune2)
tuned2 <- tune.svm(Group ~ tsne1 + tsne2, data = traindat2, gamma = 0.5, cost = 0.1, kernel = "polynomial",
                 tunecontrol=tune.control(cross=2))
summary(tuned2)
col <- c("tsne1", "tsne2", "Group")
svm.model2 <- svm(Group ~ tsne1 + tsne2, data = traindat2, cost = 0.1, gamma = 0.5, kernel="polynomial",
                  probability=TRUE)
plot(svm.model2, traindat2[,col])
svm.pred2  <- predict(svm.model2, testtsne2,probability=TRUE)
plot(svm.pred2)
tableMat <- table(pred = svm.pred2, true = test2$Group)
mean(svm.pred2==traindat2$Group)
epi.tests(tableMat)




svm_tune3 <- tune(svm, train.x=traindat3$tsne1 + traindat3$tsne2 + traindat3$tsne3, 
                  train.y=traindat3$Group, 
                  kernel="polynomial", probability=TRUE, ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))
summary(svm_tune3)
tuned3 <- tune.svm(Group ~ tsne1 + tsne2 + tsne3, data = traindat3, gamma = 0.5, cost = 0.1, 
                   kernel = "polynomial", probability=TRUE,
                   tunecontrol=tune.control(cross=2))
summary(tuned3)
col <- c("Group", "tsne1", "tsne2", "tsne3")
svm.model3 <- svm(Group ~ tsne1 + tsne2 + tsne3, data = traindat3, cost = 0.1, gamma = 0.5, 
                  kernel="polynomial", probability=TRUE)


plot3d(x=traindat3$tsne1,y=traindat3$tsne2,z=traindat3$tsne3,
       xlab = "t-SNE1", ylab = "t-SNE2", zlab = "t-SNE3",
       col=c("red","green")[traindat3[,1]],
       type="s",radius=2.5)
legend3d("topright", legend = paste(c('MDD Suicide', 'NPS Death')), 
         pch = 16, col = c("red","green"), cex=2, inset=c(0.02))
#snapshot3d(filename = 'tsne3d.png', fmt = 'png')

newdat.list = lapply(testtsne3, function(x) seq(min(x), max(x), len=50))
newdat      = expand.grid(newdat.list)
newdat.pred = predict(svm.model3, newdata=newdat, decision.values=T)
newdat.dv   = attr(newdat.pred, 'decision.values')
newdat.dv   = array(newdat.dv, dim=rep(50, 3))

contour3d(newdat.dv, level=0, x=newdat.list$tsne1, y=newdat.list$tsne2, z=newdat.list$tsne3, add=T)

svm.pred3  <- predict(svm.model3, testtsne3, probability=TRUE)
plot(svm.pred3)
tableMat <- table(pred = svm.pred3, true = test3$Group)
mean(svm.pred3==test3$Group)
epi.tests(tableMat)



svmmod2.rocr<-prediction(attr(svm.pred2,"probabilities")[,2], test2$Group)
svmmod3.rocr<-prediction(attr(svm.pred3,"probabilities")[,2], test3$Group)


svm.perf2<-performance(svmmod2.rocr, measure = "tpr", x.measure = "fpr")
svm.perf3<-performance(svmmod3.rocr, measure = "tpr", x.measure = "fpr")

plot(svm.perf2,col=1)
plot(svm.perf3,add=TRUE,col=2)
legend(0.4, 0.3, c('2 dimensional t-SNE: AUC = 95%', '3 dimensional t-SNE: AUC = 100%'), 1:2)


svm.auc2<-as.numeric(performance(svmmod2.rocr, measure = "auc", x.measure = "cutoff")@ y.values)
svm.auc3<-as.numeric(performance(svmmod3.rocr, measure = "auc", x.measure = "cutoff")@ y.values)


svm.auc2
svm.auc3


