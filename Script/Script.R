library(phytools)
library(geiger)
library(caper)
library(ape)
library(nlme)
library(mvMORPH)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(nortest)
library(multcomp)
library(RColorBrewer)

#Load tree file
tree = read.nexus("TreeTE.tre")
tree = ladderize(tree)
plot(tree)

#Load data
data = read.csv(file="Data.csv",row.names=1,header=TRUE,stringsAsFactors=TRUE)

#Check correspondence between tree tips and data
name.check(tree,data)

################################
#Relationship between Bm and Q.#
################################

#OLS regression
OLS=lm(LogQ.~LogBm,data)
summary(OLS)
ggplot(data, aes(x = LogBm, y = LogQ.)) +
  geom_point(aes(x = LogBm, y = LogQ., color=Group, shape=Group)) +
  geom_smooth(method="lm" , color="skyblue4", lty=5, se=TRUE) +
  theme_classic() +
  scale_colour_brewer(palette="Set2")

#Fit evolutionary models to Bm and Q. data
fitBM = fitContinuous(tree,data[c(1,4)],model="BM")
fitOU = fitContinuous(tree,data[c(1,4)],model="OU")
fitWN = fitContinuous(tree,data[c(1,4)],model="white")

#Check best fit of different evolutionary models for Q. and Bm using Akaike weights
aic.valsBm<-setNames(c(fitBM$LogBm$opt$aicc,fitOU$LogBm$opt$aicc,fitWN$LogBm$opt$aicc),c("BM","OU","WN"))
aic.valsBm
aic.w(aic.valsBm)
aic.valsQ.<-setNames(c(fitBM$LogQ.$opt$aicc,fitOU$LogQ.$opt$aicc,fitWN$LogQ.$opt$aicc),c("BM","OU","WN"))
aic.valsQ.
aic.w(aic.valsQ.)

#PGLS using OU model
model1<-gls(LogQ.~LogBm,data=data,correlation = corMartins(value=2.5, phy=tree), method = "ML")
anova(model1)
summary(model1)
model2<-gls(LogQ.~LogBm*Group,data=data,correlation = corMartins(value=2.5, phy=tree), method = "ML")
anova(model2)
summary(model2)
model3<-gls(LogQ.~LogBm+Group,data=data,correlation = corMartins(value=2.5, phy=tree), method = "ML")
anova(model3)
summary(model3)

#Check significance of different models
anova(model1,model2,model3)

#Check for normality of transformed residuals according to Butler et al. (2000)
transfResiduals<-chol(solve(vcv(tree)))%*%residuals(model3)
hist <- data %>%
  ggplot( aes(x=transfResiduals)) +
  geom_histogram(fill="cyan4", color="lightgray", alpha=0.9) +
  ggtitle("Transformed Residuals") +
  theme_ipsum() +
  theme(plot.title = element_text(size=15))
plot(hist)
lillie.test(transfResiduals)
#Res-Outlier <- transfResiduals[-which(rownames(transfResiduals) == "Lobodon_carcinophaga"), ]
#lillie.test(Res-Outlier)

#Post-Hoc pairwise tests for Phylogenetic ANCOVA
post.hoc<-glht(model3,linfct=mcp(Group="Tukey"))
summary(post.hoc, test = adjusted("bonferroni")) #P-corrected according to the Bonferroni method

#Plot PGLS regressions
pred = predict(model3, interval="confidence")
ggplot(data, aes(x = LogBm, y = LogQ., color = Group, shape=Group) ) +
  geom_point() +
  geom_line(aes(y = pred), linewidth = 1) +
  theme_classic() +
  scale_colour_brewer(palette="Set2")

###############################################
#Estimation of MMR values using Q. and Bm data#
###############################################

#Data preparation
data <- as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)
data[,1] = as.numeric(data[,1])
data[,4] = as.numeric(data[,4])
data[,6] = as.numeric(data[,6])

#Use only taxa with MMR data
data.MMR = data[!is.na(data$LogMMR) | data$Group == "Xenarthra" | data$Group == "Xenarthra fossil (Glyptodontinae+Folivora)",]
spmatch <- match(tree$tip.label, row.names(data.MMR))
spmatch
tree.prun <- drop.tip(tree,tree$tip.label[is.na(spmatch)])

#Fit multivariate evolutionary models to data
mvBM = mvBM(tree=tree.prun, data=data.MMR[c(1,4,6)], error = NULL, model = "BM1")
mvOU = mvOU(tree=tree.prun, data=data.MMR[c(1,4,6)], error = NULL, model = "OU1", param = list(root=FALSE, sigma = NULL, alpha = NULL, vcv = "fixedRoot", decomp = "spherical", decompSigma = "spherical"), method = "rpf", scale.height = TRUE, optimization = "L-BFGS-B", control = list(maxit = 100000), precalc = NULL,
            diagnostic = TRUE, echo = TRUE)

aicw(c(mvBM$AICc,mvOU$AICc))

#Estimate MMR values 
data.estim = data.MMR
estim = estim(tree=tree.prun, data=data.estim[c(1,4,6)], object=mvOU, error=NULL, asr=FALSE)
estim$estimates
Lower = estim$estimates - estim$se*1.96
Upper = estim$estimates + estim$se*1.96

#Cross-validate
data.cross = data.MMR
data.cross[,6] = NA
estim.cross = estim(tree=tree.prun, data=data.cross[c(1,4,6)], object=mvOU, error=NULL, asr=FALSE)
Lower.cross = estim.cross$estimates - estim.cross$se*1.96
Upper.cross = estim.cross$estimates + estim.cross$se*1.96

cross = cbind.data.frame(data.MMR[,6L],estim.cross$estimates[,3L])
cross <- setNames(cross,c("Observed", "Predicted"))
cross <- na.omit(cross)
lm = lm(cross$Predicted~cross$Observed)
summary(lm)
cross.pred = predict(lm)
ggplot(cross, aes(x = Observed, y = Predicted) ) +
  geom_point() +
  geom_line(aes(y = cross.pred), linewidth = 1) +
  scale_colour_brewer(palette="Set2")

#Ancestral state reconstruction of MMR values
taxa = row.names(data.estim)
estimMMR = cbind.data.frame(taxa,estim$estimates[,3],Lower[,3],Upper[,3])

anc = estim(tree=tree.prun, data=data.estim[c(1,4,6)], object=mvOU, error=NULL, asr=TRUE)
fit <- contMap(tree.prun,(log10(10^estim$estimates[,3]/(10^estim$estimates[,2])^0.67)),plot=FALSE,method="user",anc.states=(log10(10^anc$estimates[,3]/(10^anc$estimates[,2])^0.67)))
cols <- c("darkorchid4","dodgerblue3","cyan3","lightskyblue","lightcyan","white","lightyellow","lightsalmon","coral","orangered2","firebrick")
fit <- setMap(fit,cols)
plot(fit,fsize=0.7)

#Plot estimated MMR values and confidence intervals
estimMMR[,1] = as.factor(estimMMR[,1])
estimMMR[,1] = factor(estimMMR[,1], levels = row.names(estimMMR))

ggplot(estimMMR) +
  geom_segment( aes(x=taxa, xend=taxa, y=(10^Lower.cross[,3]/(10^estim.cross$estimates[,2])^0.67), yend=(10^Upper.cross[,3]/(10^estim.cross$estimates[,2])^0.67)), color="cadetblue", alpha=0.3, linewidth=2) +
  geom_point( aes(x=taxa, y=(10^Lower.cross[,3]/(10^estim.cross$estimates[,2])^0.67)), color="lightskyblue4", alpha=0.8, size=1) +
  geom_point( aes(x=taxa, y=(10^Upper.cross[,3]/(10^estim.cross$estimates[,2])^0.67)), color="lightskyblue4", alpha=0.8, size=1) +
  geom_point( aes(x=taxa, y=(10^estim.cross$estimates[,3]/(10^estim.cross$estimates[,2])^0.67)), color="firebrick4", alpha=0.5, size=3) +
  geom_point( aes(x=taxa, y=(10^data.MMR[,6]/(10^data.MMR[,4])^0.67)), color="gray0", alpha=1, size=2 ) +
  coord_flip()+
  theme_ipsum() +
  theme(legend.position = "none",) +
  xlab("") +
  scale_y_continuous(trans='log10') +
  ylab("MMR")


#--------------------------------------------#

###################################################
#Duplicating analysis using the morphological tree#
###################################################

#Load morpho tree file
treeMorpho = read.nexus("TreeMorph.tre")
treeMorpho = ladderize(treeMorpho)
plot(treeMorpho)

#Check correspondence between tree tips and data
name.check(treeMorpho,data)

################################
#Relationship between Bm and Q.#
################################

#OLS regression
OLS=lm(LogQ.~LogBm,data)
summary(OLS)
ggplot(data, aes(x = LogBm, y = LogQ.)) +
  geom_point(aes(x = LogBm, y = LogQ., shape=Group)) +
  geom_smooth(method="lm" , color="skyblue4", lty=5, se=TRUE) +
  scale_colour_brewer(palette="Set2")

#Fit evolutionary models to Bm and Q. data
fitBM = fitContinuous(treeMorpho,data[c(1,4)],model="BM")
fitOU = fitContinuous(treeMorpho,data[c(1,4)],model="OU")
fitWN = fitContinuous(treeMorpho,data[c(1,4)],model="white")

#Check best fit of different evolutionary models for Q. and Bm using Akaike weights
aic.valsBm<-setNames(c(fitBM$LogBm$opt$aicc,fitOU$LogBm$opt$aicc,fitWN$LogBm$opt$aicc),c("BM","OU","WN"))
aic.valsBm
aic.w(aic.valsBm)
aic.valsQ.<-setNames(c(fitBM$LogQ.$opt$aicc,fitOU$LogQ.$opt$aicc,fitWN$LogQ.$opt$aicc),c("BM","OU","WN"))
aic.valsQ.
aic.w(aic.valsQ.)

#PGLS using OU model
model1<-gls(LogQ.~LogBm,data=data,correlation = corBrownian(value=2.5, phy=treeMorpho), method = "ML")
anova(model1)
summary(model1)
model2<-gls(LogQ.~LogBm*Group,data=data,correlation = corMartins(value=2.5, phy=treeMorpho), method = "ML")
anova(model2)
summary(model2)
model3<-gls(LogQ.~LogBm+Group,data=data,correlation = corMartins(value=2.5, phy=treeMorpho), method = "ML")
anova(model3)
summary(model3)

#Check significance of different models
anova(model1,model2,model3)

#Check for normality of transformed residuals according to Butler et al. (2000)
transfResiduals<-chol(solve(vcv(treeMorpho)))%*%residuals(model3)
hist <- data %>%
  ggplot( aes(x=transfResiduals)) +
  geom_histogram(fill="cyan4", color="lightgray", alpha=0.9) +
  ggtitle("Transformed Residuals") +
  theme_ipsum() +
  theme(plot.title = element_text(size=15))
plot(hist)
lillie.test(transfResiduals)
#ResOutlier <- transfResiduals[-which(rownames(transfResiduals) == "Lobodon_carcinophaga"), ]
#lillie.test(ResOutlier)

#Post-Hoc pairwise tests for Phylogenetic ANCOVA
post.hoc<-glht(model3,linfct=mcp(Group="Tukey"))
summary(post.hoc, test = adjusted("bonferroni")) #P-corrected according to the Bonferroni method

#Plot PGLS regressions
pred = predict(model3)
ggplot(data, aes(x = LogBm, y = LogQ., color = Group) ) +
  geom_point() +
  geom_line(aes(y = pred), linewidth = 1) +
  scale_colour_brewer(palette="Set2")

###############################################
#Estimation of MMR values using Q. and Bm data#
###############################################

#Data preparation
data <- as.data.frame(treedata(treeMorpho,data,sort=T,warnings=T)$data)
data[,1] = as.numeric(data[,1])
data[,4] = as.numeric(data[,4])
data[,6] = as.numeric(data[,6])

#Use only taxa with MMR data
data.MMR = data[!is.na(data$LogMMR) | data$Group == "Xenarthra" | data$Group == "Xenarthra fossil (Glyptodontinae+Folivora)",]
spmatch <- match(treeMorpho$tip.label, row.names(data.MMR))
spmatch
treeMorpho.prun <- drop.tip(treeMorpho,treeMorpho$tip.label[is.na(spmatch)])

#Fit multivariate evolutionary models to data
mvBM = mvBM(tree=treeMorpho.prun, data=data.MMR[c(1,4,6)], error = NULL, model = "BM1")
mvOU = mvOU(tree=treeMorpho.prun, data=data.MMR[c(1,4,6)], error = NULL, model = "OU1", param = list(root=FALSE, sigma = NULL, alpha = NULL, vcv = "fixedRoot", decomp = "spherical", decompSigma = "spherical"), method = "rpf", scale.height = TRUE, optimization = "L-BFGS-B", control = list(maxit = 100000), precalc = NULL,
            diagnostic = TRUE, echo = TRUE)

aicw(c(mvBM$AICc,mvOU$AICc))

#Estimate MMR values 
data.estim = data.MMR
estim = estim(tree=treeMorpho.prun, data=data.estim[c(1,4,6)], object=mvOU, error=NULL, asr=FALSE)
estim$estimates
Lower = estim$estimates - estim$se*1.96
Upper = estim$estimates + estim$se*1.96

#Cross-validate
data.cross = data.MMR
data.cross[,6] = NA
estim.cross = estim(tree=treeMorpho.prun, data=data.cross[c(1,4,6)], object=mvOU, error=NULL, asr=FALSE)
Lower.cross = estim.cross$estimates - estim.cross$se*1.96
Upper.cross = estim.cross$estimates + estim.cross$se*1.96

cross = cbind.data.frame(data.MMR[,6L],estim.cross$estimates[,3L])
cross <- setNames(cross,c("Observed", "Predicted"))
cross <- na.omit(cross)
lm = lm(cross$Predicted~cross$Observed)
summary(lm)
cross.pred = predict(lm)
ggplot(cross, aes(x = Observed, y = Predicted) ) +
  geom_point() +
  geom_line(aes(y = cross.pred), linewidth = 1) +
  scale_colour_brewer(palette="Set2")

#Ancestral state reconstruction of MMR values
taxa = row.names(data.estim)
estimMMR = cbind.data.frame(taxa,estim$estimates[,3],Lower[,3],Upper[,3])

anc = estim(tree=treeMorpho.prun, data=data.estim[c(1,4,6)], object=mvOU, error=NULL, asr=TRUE)
fit <- contMap(treeMorpho.prun,(log10(10^estim$estimates[,3]/(10^estim$estimates[,2])^0.67)),plot=FALSE,method="user",anc.states=(log10(10^anc$estimates[,3]/(10^anc$estimates[,2])^0.67)))
cols <- c("darkorchid4","dodgerblue3","cyan3","lightskyblue","lightcyan","white","lightyellow","lightsalmon","coral","orangered2","firebrick")
fit <- setMap(fit,cols)
plot(fit,fsize=0.7)

#Plot estimated MMR values and confidence intervals
estimMMR[,1] = as.factor(estimMMR[,1])
estimMMR[,1] = factor(estimMMR[,1], levels = row.names(estimMMR))

ggplot(estimMMR) +
  geom_segment( aes(x=taxa, xend=taxa, y=(10^Lower.cross[,3]/(10^estim.cross$estimates[,2])^0.67), yend=(10^Upper.cross[,3]/(10^estim.cross$estimates[,2])^0.67)), color="cadetblue", alpha=0.3, linewidth=2) +
  geom_point( aes(x=taxa, y=(10^Lower.cross[,3]/(10^estim.cross$estimates[,2])^0.67)), color="lightskyblue4", alpha=0.8, size=1) +
  geom_point( aes(x=taxa, y=(10^Upper.cross[,3]/(10^estim.cross$estimates[,2])^0.67)), color="lightskyblue4", alpha=0.8, size=1) +
  geom_point( aes(x=taxa, y=(10^estim.cross$estimates[,3]/(10^estim.cross$estimates[,2])^0.67)), color="firebrick4", alpha=0.5, size=3) +
  geom_point( aes(x=taxa, y=(10^data.MMR[,6]/(10^data.MMR[,4])^0.67)), color="gray0", alpha=1, size=2 ) +
  coord_flip()+
  theme_ipsum() +
  theme(legend.position = "none",) +
  xlab("") +
  scale_y_continuous(trans='log10') +
  ylab("MMR")


#--------------------------------------------#

##############################################
#Duplicating analysis using Qi instead of Q.#
##############################################

################################
#Relationship between Bm and Qi#
################################

#OLS regression
OLS=lm(LogQi~LogBm,data)
summary(OLS)
ggplot(data, aes(x = LogBm, y = LogQi)) +
  geom_point(aes(x = LogBm, y = LogQi, shape=Group)) +
  geom_smooth(method="lm" , color="skyblue4", lty=5, se=TRUE) +
  scale_colour_brewer(palette="Set2")

#Fit evolutionary models to Bm and Qi data
fitBM = fitContinuous(tree,data[c(2,4)],model="BM")
fitOU = fitContinuous(tree,data[c(2,4)],model="OU")
fitWN = fitContinuous(tree,data[c(2,4)],model="white")

#Check best fit of different evolutionary models for Qi and Bm using Akaike weights
aic.valsBm<-setNames(c(fitBM$LogBm$opt$aicc,fitOU$LogBm$opt$aicc,fitWN$LogBm$opt$aicc),c("BM","OU","WN"))
aic.valsBm
aic.w(aic.valsBm)
aic.valsQi<-setNames(c(fitBM$LogQi$opt$aicc,fitOU$LogQi$opt$aicc,fitWN$LogQi$opt$aicc),c("BM","OU","WN"))
aic.valsQi
aic.w(aic.valsQi)

#PGLS using OU model
model1<-gls(LogQi~LogBm,data=data,correlation = corBrownian(value=4.1, phy=tree), method = "ML")
anova(model1)
summary(model1)
model2<-gls(LogQi~LogBm*Group,data=data,correlation = corMartins(value=4.1, phy=tree), method = "ML")
anova(model2)
summary(model2)
model3<-gls(LogQi~LogBm+Group,data=data,correlation = corMartins(value=4.1, phy=tree), method = "ML")
anova(model3)
summary(model3)

#Check significance of different models
anova(model1,model2,model3)

#Check for normality of transformed residuals according to Butler et al. (2000)
transfResiduals<-chol(solve(vcv(tree)))%*%residuals(model3)
hist <- data %>%
  ggplot( aes(x=transfResiduals)) +
  geom_histogram(fill="cyan4", color="lightgray", alpha=0.9) +
  ggtitle("Transformed Residuals") +
  theme_ipsum() +
  theme(plot.title = element_text(size=15))
plot(hist)
lillie.test(transfResiduals)
#Res-Outlier <- transfResiduals[-which(rownames(transfResiduals) == "Lobodon_carcinophaga"), ]
#lillie.test(Res-Outlier)

#Post-Hoc pairwise tests for Phylogenetic ANCOVA
post.hoc<-glht(model3,linfct=mcp(Group="Tukey"))
summary(post.hoc, test = adjusted("bonferroni")) #P-corrected according to the Bonferroni method

#Plot PGLS regressions
pred = predict(model3)
ggplot(data, aes(x = LogBm, y = LogQi, color = Group) ) +
  geom_point() +
  geom_line(aes(y = pred), linewidth = 1) +
  scale_colour_brewer(palette="Set2")

###############################################
#Estimation of MMR values using Qi and Bm data#
###############################################

#Data preparation
data <- as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)
data[,2] = as.numeric(data[,2])
data[,4] = as.numeric(data[,4])
data[,6] = as.numeric(data[,6])

#Use only taxa with MMR data
data.MMR = data[!is.na(data$LogMMR) | data$Group == "Xenarthra" | data$Group == "Xenarthra fossil (Glyptodontinae+Folivora)",]
spmatch <- match(tree$tip.label, row.names(data.MMR))
spmatch
tree.prun <- drop.tip(tree,tree$tip.label[is.na(spmatch)])

#Fit multivariate evolutionary models to data
mvBM = mvBM(tree=tree.prun, data=data.MMR[c(2,4,6)], error = NULL, model = "BM1")
mvOU = mvOU(tree=tree.prun, data=data.MMR[c(2,4,6)], error = NULL, model = "OU1", param = list(root=FALSE, sigma = NULL, alpha = NULL, vcv = "fixedRoot", decomp = "spherical", decompSigma = "spherical"), method = "rpf", scale.height = TRUE, optimization = "L-BFGS-B", control = list(maxit = 100000), precalc = NULL,
            diagnostic = TRUE, echo = TRUE)

aicw(c(mvBM$AICc,mvOU$AICc))

#Estimate MMR values 
data.estim = data.MMR
estim = estim(tree=tree.prun, data=data.estim[c(2,4,6)], object=mvOU, error=NULL, asr=FALSE)
estim$estimates
Lower = estim$estimates - estim$se*1.96
Upper = estim$estimates + estim$se*1.96

#Cross-validate
data.cross = data.MMR
data.cross[,6] = NA
estim.cross = estim(tree=tree.prun, data=data.cross[c(2,4,6)], object=mvOU, error=NULL, asr=FALSE)
Lower.cross = estim.cross$estimates - estim.cross$se*1.96
Upper.cross = estim.cross$estimates + estim.cross$se*1.96

cross = cbind.data.frame(data.MMR[,6L],estim.cross$estimates[,3L])
cross <- setNames(cross,c("Observed", "Predicted"))
cross <- na.omit(cross)
lm = lm(cross$Predicted~cross$Observed)
summary(lm)
cross.pred = predict(lm)
ggplot(cross, aes(x = Observed, y = Predicted) ) +
  geom_point() +
  geom_line(aes(y = cross.pred), linewidth = 1) +
  scale_colour_brewer(palette="Set2")

#Ancestral state reconstruction of MMR values
taxa = row.names(data.estim)
estimMMR = cbind.data.frame(taxa,estim$estimates[,3],Lower[,3],Upper[,3])

anc = estim(tree=tree.prun, data=data.estim[c(2,4,6)], object=mvOU, error=NULL, asr=TRUE)
fit <- contMap(tree.prun,(log10(10^estim$estimates[,3]/(10^estim$estimates[,2])^0.67)),plot=FALSE,method="user",anc.states=(log10(10^anc$estimates[,3]/(10^anc$estimates[,2])^0.67)))
cols <- c("darkorchid4","dodgerblue3","cyan3","lightskyblue","lightcyan","white","lightyellow","lightsalmon","coral","orangered2","firebrick")
fit <- setMap(fit,cols)
plot(fit,fsize=0.7)

#Plot estimated MMR values and confidence intervals
estimMMR[,1] = as.factor(estimMMR[,1])
estimMMR[,1] = factor(estimMMR[,1], levels = row.names(estimMMR))

ggplot(estimMMR) +
  geom_segment( aes(x=taxa, xend=taxa, y=(10^Lower.cross[,3]/(10^estim.cross$estimates[,2])^0.67), yend=(10^Upper.cross[,3]/(10^estim.cross$estimates[,2])^0.67)), color="cadetblue", alpha=0.3, linewidth=2) +
  geom_point( aes(x=taxa, y=(10^Lower.cross[,3]/(10^estim.cross$estimates[,2])^0.67)), color="lightskyblue4", alpha=0.8, size=1) +
  geom_point( aes(x=taxa, y=(10^Upper.cross[,3]/(10^estim.cross$estimates[,2])^0.67)), color="lightskyblue4", alpha=0.8, size=1) +
  geom_point( aes(x=taxa, y=(10^estim.cross$estimates[,3]/(10^estim.cross$estimates[,2])^0.67)), color="firebrick4", alpha=0.5, size=3) +
  geom_point( aes(x=taxa, y=(10^data.MMR[,6]/(10^data.MMR[,4])^0.67)), color="gray0", alpha=1, size=2 ) +
  coord_flip()+
  theme_ipsum() +
  theme(legend.position = "none",) +
  xlab("") +
  scale_y_continuous(trans='log10') +
  ylab("MMR")

