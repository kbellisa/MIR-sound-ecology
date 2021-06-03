# MIR ANALYSIS
# Kristen Bellisario and Jack VanSchaik
# Final: September 13, 2017

#install.packages("alluvial")
library(vegan) #permutational anova
library(MASS) #LDA and QDA
library(colorspace) #new colors
library(caret) #sensitivity
library(dplyr) #alluvial sort
library(alluvial) #plot
library(gplots) #heatmap


# Load Data / Prep ---------------------------------------------------------------


#IMPORT ALL FEATURES FOR WINDOW OBSERVATION LENGTHS
feature <- read.csv("/Users/kristenbellisario/Documents/R_Dir/R_SoundMetricProject/TestingPhase/TestFeatures_180.csv", header=T)

#IMPORT LABELS ASSIGNED VIA SEM
predictor <- read.csv("/Users/kristenbellisario/Documents/R_Dir/R_SoundMetricProject/TestingPhase/Master_Codes2.csv", header=T)
  # pull out full label set that doesn't exclude "outliers"
  cond <- predictor[,8]
  # correction of data for row 155
  cond[155] <- 4

#COMBINE FEATURE AND LABEL IN DATAFRAME 
data1 <- cbind(cond, feature)
co.fl <- as.factor(data1[,1])

# SELECTED FEATURE GROUPS
# FG4 - centroid, skewness
# FG5 - centroid, skewness, variance
# FG11 - centroid, slope, spread
# FG12 - centroid, slope
# FG26 - spread, variance

# FEATURES centroid 6, skew 9, slope 12, spread 15, variance 18
#CREATE LIST OF ALL COMBINATIONS OF FEATURES
L <- list(
  data1[,c(6,9,12,15,18)],
  data1[,c(6,9,12,15)],
  data1[,c(6,9,12)],
  data1[,c(6,9)],
  data1[,c(6,9,15)],
  data1[,c(6,9,18)],
  data1[,c(6,9,12,18)],
  data1[,c(6,9,15,18)],
  data1[,c(6,12,15,18)],
  data1[,c(6,12,18)],
  data1[,c(6,12,15)],
  data1[,c(6,12)],
  data1[,c(6,15,18)],
  data1[,c(6,15)],
  data1[,c(6,18)],
  data1[,c(9,12,15,18)],
  data1[,c(9,12,15)],
  data1[,c(9,12,18)],
  data1[,c(9,12)],
  data1[,c(9,18)],
  data1[,c(9,15,18)],
  data1[,c(9,15)],
  data1[,c(12,15,18)],
  data1[,c(12,15)],
  data1[,c(12,18)],
  data1[,c(15,18)]
)


# Permutational Anova Results ---------------------------------------------


# CHECK PERMUTATIONAL ANOVA RESULTS TO DETERMINE WHICH FEATURE GROUPS ARE SIGNIFICANT WITH BONFERRONI CORRECTION (p=0.001)
perm.results <- lapply(L, function(x) {
  class.adon1 <- adonis(x ~ co.fl, method="euclid", data=data1,
                        control=permControl(strata=col), permutations=9999, by="terms")
  return(class.adon1)
})


# Beta Dispersion between groups ------------------------------------------


# CHECK BETA DISPERSION BETWEEN GROUPS
adon.disp <- lapply(L, function(y) {
  betadisper(vegdist(x, method="euclid"), co.fl)
  bp <- boxplot(adon.disp, col="blue", main="Beta Dispersion", labels=names, cex.main =1, cex.axis=.8, cex.lab=0.8)
  paste( 'SAVE THE IMAGE' )
  return(adon.disp)
})


# LDA with Cross-validation / temp variable -------------------------------------------------------------------------

z <- data1[,c(6,9)]

# RUN LDA WITH CROSS-VALIDATION
lda.results <- lapply(L, function(z) {
  ld.5.3.cv <- lda(z, co.fl, CV=T)
  assess.l <- table(co.fl, ld.5.3.cv$class) # compute accuracy step 1
  diag(prop.table(assess.l,1)) # compute accuracy step 2
  sum(diag(prop.table(assess.l))) # compute accuracy step 3
  pred <- ld.5.3.cv$class # establish predictor class from lda
  test1 <- table(pred, cond) # create table of predictor LDA class and SEM class
  sensitivity(test1, "4") # compute sens of biophony
  specificity(test1, c("6", "7")) #compute specificity of biophony
  sensitivity(test1, "6") # compute sens of geophony
  specificity(test1, c("4", "7")) # compute spec of geophony
  sensitivity(test1, "7") # compute sens of anthrophony
  specificity(test1, c("4", "6"))   #compute spec of anthrophony
})
  return(lda.results)
})


# QDA with Cross-validation -----------------------------------------------


# RUN QDA WITH CROSS-VALIDATION
qda.results <- lapply(L, function(xx) {
  qd.5.3.cv <- qda(xx, co.fl, CV=T)
  assess.q <- table(co.fl, qd.5.3.cv$class) #compute accuracy step 1
  diag(prop.table(assess.q,1)) #compute accuracy step 2
  sum(diag(prop.table(assess.q))) #compute accuracy step 3
  pred <- qd.5.3.cv$class # establish predictor class from lda
  test2 <- table(pred, cond) # create table of predictor LDA class and SEM class
  sensitivity(test2, "4") # compute sens of biophony
  specificity(test2, c("6", "7")) #compute specificity of biophony
  sensitivity(test2, "6") # compute sens of geophony
  specificity(test2, c("4", "7")) # compute spec of geophony
  sensitivity(test2, "7") # compute sens of anthrophony
  specificity(test2, c("4", "6"))   #compute spec of anthrophony
  return(qda.results)
})


# Plot SEM Class Distribution - Alluvial ----------------------------------

# PLOT SEM CLASS DISTRIBUTION USING ALLUVIAL PLOT

raw.feature <- read.csv("/Users/kristenbellisario/Documents/R_Dir/R_SoundMetricProject/TestingPhase/MIR_masterdata_appendix.csv", header=T)
names(raw.feature) <- c("ID", "File", "B", "G", "A", "SEM", "Class", "Time", "Site")
name.label <- ifelse(raw.feature$Class == 4, "biophony", ifelse(raw.feature$Class == 6, "geophony", "anthrophony"))
#carry forward data1 / 

rf<-raw.feature %>% group_by(Class, Time, Site) %>% summarise(numObs=length(Class))

alluvial(
  rf[,1:3],
  freq=rf$numObs,
  col = ifelse(rf$Class == 'biophony', "mediumseagreen", ifelse(rf$Class == 'geophony', "cadetblue1","tomato3")),
  cex=.5,
  cex.axis = .7,
  alpha = 0.5,
  border="gray",
  gap.width=.1,
  blocks=FALSE
)

#gray version
alluvial(
  rf[,1:3],
  freq=rf$numObs,
  col = ifelse(rf$Class == 'biophony', "gray20", ifelse(rf$Class == 'geophony', "gray54","gray75")),
  cex=.5,
  cex.axis = .7,
  alpha = 0.7,
  border="gray",
  gap.width=.1,
  blocks=FALSE
)


# Combos for LDA ----------------------------------------------------------

#6 cen 9 skew 12 sl 15 sp 18 var
#LDA
FG1 <- data1[,c(6,9,12,15,18)] #NO
FG2 <- data1[,c(6,9,12,15)] # NO
FG3 <- data1[,c(6,9,12)] # NO 6 49 18
FG4 <- data1[,c(6,9)] # Y 1 54 18 ---- GEOPHONY
FG5 <- data1[,c(6,9,15)] # 1 53 19
FG6 <- data1[,c(6,9,18)] # NO
FG7 <- data1[,c(6,9,12,18)] #N
FG8 <- data1[,c(6,9,15,18)] #N
FG9 <- data1[,c(6,12,15,18)] #N
FG10 <- data1[,c(6,12,18)] #N
FG11 <- data1[,c(6,12,15)] #Y 0 54 19 GEO  / #Y for ANTHRO
FG12 <- data1[,c(6,12)] #Y 0 53 20
FG13 <- data1[,c(6,15,18)] #NULL
FG14 <- data1[,c(6,15)] #N
FG15 <- data1[,c(6,18)] #N
FG16 <- data1[,c(9,12,15,18)] #N
FG17 <- data1[,c(9,12,15)] #N
FG18 <- data1[,c(9,12,18)] #N
FG19 <- data1[,c(9,12)] #N
FG20 <- data1[,c(9,18)] #N
FG21 <- data1[,c(9,15,18)] #N
FG22 <- data1[,c(9,15)] # 73 for geophony was PERFECT
FG23 <- data1[,c(12,15,18)] #N
FG24 <- data1[,c(12,15)] #1 54 18 SAME
FG25 <- data1[,c(12,18)] #N
FG26 <- data1[,c(15,18)] #N


# Assess using Binary Classification --------------------------------------

# ASSESS USING BINARY CLASSIFICATION
co.fl.new <- read.csv("/Users/kristenbellisario/Documents/R_Dir/R_SoundMetricProject/co.fl.complete.csv")
co.geo <- as.factor (co.fl.new[,3])
co.anthro <- as.factor (co.fl.new[,4])
co.bio <- as.factor (co.fl.new[,5])

lda.geo <- qda(FG13, co.bio, CV=T)
lda.geo.o <- lda(FG4, co.anthro)

assess.q <- table(co.bio, lda.geo$class) #compute accuracy step 1
diag(prop.table(assess.q,1)) #compute accuracy step 2
sum(diag(prop.table(assess.q))) #compute accuracy step 3
# TN (NN) = bottom left //FP (NY) = top right
# FN (YN) = right top //  TP (YY) = left bottom

colnames(assess.q) <- c("1", "4")
rownames(assess.q) <- c("1", "4")

heatmap.2(assess.q, dendrogram="none", trace="none", key=F, srtCol=0, Rowv=F, 
          Colv=FALSE, cexRow = .7, cexCol = .7, cellnote=assess.q, notecol="black", 
          col=colors(distinct=F), add=T, main="BIO: FG13 2")


# Compute ROC Curve -------------------------------------------------------

## TO COMPUTE ROC CURVE
library(pROC)
roc_obj <- roc(as.numeric(co.anthro), as.numeric(lda.geo.predict$class))
a <- auc(roc_obj)
plot(roc_obj, main="FG4", xlab="FPR", ylab="TPR")
text(0.0, .4, paste("auc=", round(a, digits=2)))

#RETEST GEOPHONY REMOVING WEATHER OUTLIERS
#subset of outliers: rain, thunder, etc. // all outliers were put into a class called 1
#Y <- c(-1,-7,-8,-20,-26,-53,-86, -93,-94,-98,-101,-102,-103,-104, -109,-111,-116,-117,-160,-169)

data.no <- data1[c(-1,-7,-8,-20,-26,-53,-86, -93,-94,-98,-101,-102,-103,-104, -109,-111,-116,-117,-160,-169),]
lda.geo2 <- data.frame(co.fl.new[c(-1,-7,-8,-20,-26,-53,-86, -93,-94,-98,-101,-102,-103,-104, -109,-111,-116,-117,-160,-169),3])
lda.geo2.f <- as.factor(unlist(lda.geo2))

FGn1 <- data.no[,c(6,9,12,15,18)] #NO
FGn2 <- data.no[,c(6,9,12,15)] # NO
FGn3 <- data.no[,c(6,9,12)] # NO 6 49 18
FGn4 <- data.no[,c(6,9)] # Y 1 54 18 ---- GEOPHONY
FGn5 <- data.no[,c(6,9,15)] # 1 53 19
FGn6 <- data.no[,c(6,9,18)] # NO
FGn7 <- data.no[,c(6,9,12,18)] #N
FGn8 <- data.no[,c(6,9,15,18)] #N
FGn9 <- data.no[,c(6,12,15,18)] #N
FGn10 <- data.no[,c(6,12,18)] #N
FGn11 <- data.no[,c(6,12,15)] #Y 0 54 19 GEO  / #Y for ANTHRO
FGn12 <- data.no[,c(6,12)] #Y 0 53 20
FGn13 <- data.no[,c(6,15,18)] #NULL
FGn14 <- data.no[,c(6,15)] #N
FGn15 <- data.no[,c(6,18)] #N
FGn16 <- data.no[,c(9,12,15,18)] #N
FGn17 <- data.no[,c(9,12,15)] #N
FGn18 <- data.no[,c(9,12,18)] #N
FGn19 <- data.no[,c(9,12)] #N
FGn20 <- data.no[,c(9,18)] #N
FGn21 <- data.no[,c(9,15,18)] #N
FGn22 <- data.no[,c(9,15)] # 73 for geophony was PERFECT
FGn23 <- data.no[,c(12,15,18)] #N
FGn24 <- data.no[,c(12,15)] #1 54 18 SAME
FGn25 <- data.no[,c(12,18)] #N
FGn26 <- data.no[,c(15,18)] #N

lda.geo <- lda(FGn21, lda.geo2.f, CV=T)
lda.geo <- lda(FGn21, lda.geo2.f, CV=T)

assess.q <- table(lda.geo2.f, lda.geo$class) #compute accuracy step 1
diag(prop.table(assess.q,1)) #compute accuracy step 2 (1 row percentage, 2 column percentage)
sum(diag(prop.table(assess.q))) #compute accuracy step 3

colnames(assess.q) <- c("1","6")
rownames(assess.q) <- c("other","geo")

heatmap.2(assess.q, dendrogram="none", trace="none", key=F, srtCol=0, Rowv=F, 
          Colv=FALSE, cexRow = .7, cexCol = .7, cellnote=assess.q, notecol="black", 
          col=colors(distinct=T), add=T, main="FILTER GEOPHONY: lFGn21")

