# rm (list = ls( )) 
setwd("/Users/rifatzahan/Desktop/Working Files/CMPT 830/Project/GSE88890_RAW/")
update.packages(ask = FALSE, dependencies = c ( 'Suggests' ))


library(foreign)
library(varhandle)
library(manhattanly)
library(ggplot2)


id_dat <- read.csv("Individuals_Description.csv", header = T)
genetic_dat <- read.csv("genetic_normalisedBetas.csv", header = T, nrows = 100)
#genetic_dat <- read.csv("GSE88890_normalisedBetas.csv", header = T, nrows = 100)
genetic_dat <- as.data.frame(t(genetic_dat))
colnames(genetic_dat) <- as.character(unlist(genetic_dat[1,]))
genetic_dat <- genetic_dat[-1, ]
genetic_dat$ID <- factor(row.names(genetic_dat))

total <- merge(id_dat,genetic_dat,by="ID")
ba11_dat <- subset(total, Cortex == "BA11")
ba25_dat <- subset(total, Cortex == "BA25")


##########
## BA11 ##
##########
ba11_dat_genes <- ba11_dat
ba11_dat_genes[1:6] <- list(NULL)
ba11_dat_genes <- unfactor(ba11_dat_genes)

ba11_dat_genes_group <- cbind(ba11_dat$Group, ba11_dat_genes)
colnames(ba11_dat_genes_group)[1] <- "Group"
#deathClasses <- factor(ba11_dat_genes_group$Group)


suicideBA11 <- subset(ba11_dat_genes_group, Group == "MDD Suicide")
nonSuicideBA11 <- subset(ba11_dat_genes_group, Group == "NPS Death")

pvaluesBA11 = c()
suicidemean = c()
normalmean = c()

for (i in 2:dim(suicideBA11)[2]){
  suicidemean <- append(suicidemean, (t.test(suicideBA11[i], nonSuicideBA11[i]))$estimate[[1]])
  normalmean <- append(normalmean, (t.test(suicideBA11[i], nonSuicideBA11[i]))$estimate[[2]])
  pvaluesBA11 <-  append(pvaluesBA11, (t.test(suicideBA11[i], nonSuicideBA11[i]))$p.value)
}

Genes <- names(suicideBA11)[-1]
ba11_dat_volcano <- as.data.frame(cbind(Genes, suicidemean, normalmean, pvaluesBA11))
ba11_dat_volcano$Diff <- (unfactor(ba11_dat_volcano$suicidemean)-unfactor(ba11_dat_volcano$normalmean))
ba11_dat_volcano$pvaluesBA11 <- unfactor(ba11_dat_volcano$pvaluesBA11)

ba11_dat_volcano[which(ba11_dat_volcano['pvaluesBA11'] <= 0.05),"Group"] <- "Significant"
ba11_dat_volcano[which(ba11_dat_volcano['pvaluesBA11'] > 0.05),"Group"] <- "Insignificant"


ggplot(ba11_dat_volcano, aes(x=Diff, y=pvaluesBA11)) +
  geom_point(aes(colour = Group), size=2.5) +
  scale_colour_manual(values = c("Significant"= "green", "Insignificant"="grey"))


#pvaluesBA11 = append(ttestresultBA11, (t.test(suicideBA11[i], nonSuicideBA11[i]))$p.value)
#diffmean = append()




ttestresultBA11 <- round(ttestresultBA11, 2)
sig_genes <- as.factor(ttestresultBA11 <= 0.05)
sig_index <- which(sig_genes == "TRUE")
sig_index

dfBA11 <- as.data.frame(cbind(ba11_dat$ID, ba11_dat$Age, ba11_dat$Sex, ba11_dat_genes_group$Group, 
                              ba11_dat_genes_group[sig_index[1]], ba11_dat_genes_group[sig_index[2]], 
                              ba11_dat_genes_group[sig_index[3]],ba11_dat_genes_group[sig_index[4]]))

colnames(dfBA11)[1] <- "ID"
colnames(dfBA11)[2] <- "Age"
colnames(dfBA11)[3] <- "Sex"
colnames(dfBA11)[4] <- "Group"

logisticBA11 <- glm(Group ~ Age + as.factor(Sex) + cg00001446 + cg00001593 + cg00003014 + cg26894438, data = dfBA11, 
                    family = binomial("logit"))
odds_ratioBA11 <- exp(cbind(OR = coef(logisticBA11), confint(logisticBA11)))
odds_ratioBA11 <- round(odds_ratioBA11, 3)
odds_ratioBA11
#



##########
## BA25 ##
##########
ba25_dat_genes <- ba25_dat
ba25_dat_genes[1:6] <- list(NULL)
ba25_dat_genes <- unfactor(ba25_dat_genes)

ba25_dat_genes_group <- cbind(ba25_dat$Group, ba25_dat_genes)
colnames(ba25_dat_genes_group)[1] <- "Group"
deathClasses <- factor(ba25_dat_genes_group$Group)


suicideBA25 <- subset(ba25_dat_genes_group, Group == "MDD Suicide")
nonSuicideBA25 <- subset(ba25_dat_genes_group, Group == "NPS Death")

ttestresultBA25 = c()
pvaluesBA25 = c()

for (i in 2:dim(suicideBA25)[2])
  ttestresultBA25 = append(ttestresultBA25, (t.test(suicideBA25[i], nonSuicideBA25[i]))$p.value)

ttestresultBA25 <- round(ttestresultBA25, 2)
sig_genes <- as.factor(ttestresult <= 0.05)
sig_index <- which(sig_genes == "TRUE")
sig_index

dfBA25 <- as.data.frame(cbind(ba25_dat$ID, ba25_dat$Age, ba25_dat$Sex, ba25_dat_genes_group$Group, 
                              ba25_dat_genes_group[sig_index[1]], ba25_dat_genes_group[sig_index[2]], 
                              ba25_dat_genes_group[sig_index[3]],ba25_dat_genes_group[sig_index[4]]))

colnames(dfBA25)[1] <- "ID"
colnames(dfBA25)[2] <- "Age"
colnames(dfBA25)[3] <- "Sex"
colnames(dfBA25)[4] <- "Group"

logisticBA25 <- glm(Group ~ Age + as.factor(Sex) + cg00001245 + cg00001349 + cg00002769 + cg26894259, 
                    data = dfBA25, 
                    family = binomial("logit"))
odds_ratioBA25 <- exp(cbind(OR = coef(logisticBA25), confint(logisticBA25)))
odds_ratioBA25 <- round(odds_ratioBA25, 3)
odds_ratioBA25
#






# https://www.r-bloggers.com/predicting-employment-related-factors-in-malaysia-a-regression-analysis-approach/
