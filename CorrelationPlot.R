
library(MASS)

setwd("/Users/kristenbellisario/Documents/R_Dir/SoundMetricProject/frame_threesecond_test")

#################################
# READ, ANALYZE, TRANSFORM DATA #
#################################

data <- read.csv(file="values.csv", header=T)

# data for 23ms to 1 second comparison
centroid <- data [,c(1,2,22)]
skew <- data[,c(1, 3, 23)]
slope <- data[,c(1, 4, 25)] #observation 48278 outlier
spread <- data[c(-10),c(1, 5, 24)] #observation 14483 also outlier
variance <- data[,c(1, 6, 26)]

# data for 23ms to 3 second comparison
centroid <- data[,c(1,2,7)]
skew <- data[,c(1, 3, 8)]
slope <- data[,c(1, 4, 9)] #observation 48278 outlier
spread <- data[c(-10),c(1, 5, 10)] #observation 14483 also outlier
variance <- data[,c(1, 6, 11)]

# data for 23ms to 10 second comparison
centroid <- data[,c(1,2,12)] #
skew <- data[,c(1, 3, 13)]
slope <- data[,c(1, 4, 14)] #observation 48278 outlier
spread <- data[c(-10),c(1, 5, 15)] #observation 14483 also outlier
variance <- data[,c(1, 6, 16)]

#data for 23ms to 1 minute comparison
centroid <- data[,c(1,2,17)] #
skew <- data[,c(1, 3, 18)]
slope <- data[,c(1, 4, 20)] #observation 48278 outlier
spread <- data[c(-10),c(1, 5, 19)] #observation 14483 also outlier
variance <- data[,c(1, 6, 21)]


################## CORRELATION FUNCTION #########################
# function
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.6/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r+.3)
}

#removed observation 48278
pairs(centroid[,2:3], lower.panel=panel.smooth, upper.panel=panel.cor) #0.98 / 10sec 1.0
pairs(skew[,2:3], lower.panel=panel.smooth, upper.panel=panel.cor) #0.076 - one outlier / 0.99 / 10 sec 0.92
pairs(slope[,2:3], lower.panel=panel.smooth, upper.panel=panel.cor) #0.98 / 0.20
pairs(spread[,2:3], lower.panel=panel.smooth, upper.panel=panel.cor) #0.61 ******* improved to 0.70 when removed outlier 14483 / 10 sec 0.36
pairs(variance[,2:3], lower.panel=panel.smooth, upper.panel=panel.cor) #0.92 / 10sec 0.99

#PLOT CORRELATION VALUES [centroid, skew, slope, spread, variance]
r.23_1 <- data.frame(t(c(1.0, 0.96, 0.35, .99, .39)))
r.23_3 <- data.frame(t(c(0.98, 0.99, 0.98, 0.92, 0.61)))
r.23_10 <- data.frame(t(c(1.0, 0.92, 0.20, 0.99, 0.36)))
r.23_60 <- data.frame(t(c(0.98, 0.97, 0.34, 0.98, 0.34)))
n <- c("Centroid", "Skew", "Slope", "Variance", "Spread")

M <- as.matrix(rbind(r.23_1, r.23_3, r.23_10, r.23_60))
names(M) <- n
bp <- barplot(M,names=n, beside=T, main="Comparison of Window Lengths")
legend("topright", pch=15, col=c("black", "darkgrey","grey", "lightgrey"), c("1s", "3s", "10s", "60s"))

############ FEATURE GROUP SEM VS CLUSTER EVALUATION BARPLOT
models <- c("3a", "3e", "4c", "c1", "c5")
percent <- c(0.577, 0.714, 0.627, 0.542, 0.495)

barplot(percent, names=models, beside=T, xlab = "Feature Group", ylab = "Relevance Degree")

