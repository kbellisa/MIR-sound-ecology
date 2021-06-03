##################################################################################
# MIR Study (no removal of outliers)
# USES linear and quadratic discriminant analysis
# Research question? Can we use spectral features to determine soundscape class?
# Kristen Bellisario
# with additional help from Zhao Zhao, Cristan Graupe, and Jack VanSchaik

##################################################################################
# Required Packages
library(rgl)
library(ggplot2)
library(colorspace)
library(vegan) #adonis
library(MASS) #lda
library(gplots) #heatmap

#180 sound files with 23ms and 3s window lengths (note -- label states 1s and was labeled incorrectly)
d.1 <- read.csv("TestFeatures_180.csv", header=T)
rownames(d.1) <- d.1[,1]
d.s.1 <- data.frame(cbind(d.1[,1], scale(d.1[,3:12])))
rownames(d.s.1) <- d.s.1[,1]

#COMPLETE FEATURE GROUPS FOR EACH FRAME - only running LDA on 3s model [optimized sound recording length] // including outliers

model_1.t <- cbind(d.s.1[,1], d.s.1[,c(2,4,6,8,10)])
names(model_1.t) <- c("ID", "Centroid1", "Skew1", "Slope1", "Spread1", "Var1")

#KMEANS PREDICTORS - sorted by soundfile
# each observation has three predictions when overlapping classes
# goal is to find out which class is dominant in single class membership
# note 3e is an internal naming convention for dataset group

p_set.t <- read.csv("test_data_v2.csv", header=T)
pred_p.t <- p_set.t[,c(1,3)]
factors3e.1 <- as.factor(pred_p.t[,2]) #factors 1, 3, 4, 5, 6, 7
full.3e.1.p1 <- cbind(factors3e.1, model_1.t[,c(3,5,6)])


########## NEED SPECIES FACTORS / SPECIES COUNTS
#species factors
factors3e.1

## species counts
#spc.3e.23.p1 <- full.3e.23.p1[,-1]
spc.3e.1.p1 <- full.3e.1.p1[,-1]

########### MODEL: spc.3e 

class.adon <- adonis(spc.3e.1.p1 ~ factors3e.1, method="gower", data=full.3e.1.p1,
                     control=permControl(strata=factors3e.1), permutations=999, by="terms")

###### PERMUTATION TEST

adon.disp <- betadisper(vegdist(spc.3e.1.p1, method="gower"), factors3e.1)
boxplot(adon.disp, col="blue", main="Beta Dispersion sp.3e.3.p1", cex.main =1, cex.axis=.8, cex.lab=0.8)  

###########LDA

disc.class3e.1.p1 <- lda(spc.3e.1.p1, factors3e.1) #raw

###########LDA FOR CONFUSION MATRIX / CLASS ASSESSMENT
#use quadratic discriminant analysis with outliers included
#jackknife cross validation / used for predictive value
disc.class3e.1.p1 <- qda(spc.3e.1.p1, factors3e.1, CV=T) #raw

############CONFUSION MATRIX

assess3 <- table(factors3e.1, disc.class3e.1.p1$class)  #        				
diag(prop.table(assess3,1))
sum(diag(prop.table(assess3))) 

colnames(assess3) <- c("1","3","4","5","6","7")
rownames(assess3) <- c("1","3","4","5","6","7")

heatmap.2(assess3, dendrogram="none", trace="none", key=F, srtCol=0, Rowv=F, 
          Colv=FALSE, cexRow = .7, cexCol = .7, cellnote=assess3, notecol="black", 
          col=topo.colors(50), add=T, main="Results")
