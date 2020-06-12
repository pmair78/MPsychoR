## ------------------------------ Chapter 8: Gifi Methods -------------------------------

## ----- princals
library("MPsychoR")
data("ASTI")
st <- ASTI[ ,c(2,4,7,13,16,24,25)]      
pg <- ASTI[ ,c(11,14,15,17,18,23)]  
stpg <- data.frame(st = st, pg = pg)
pcafit <- prcomp(stpg, scale = TRUE)

library("Gifi")
knotslin <- knotsGifi(stpg, type = "E")
prlin <- princals(stpg, knots = knotslin, degrees = 1)
prlin

(pcafit$sdev^2)[1:2]   
prlin$evals[1:2]
head(round(pcafit$rotation[,1:2], 3), 3)                   
head(round(prlin$loadings, 3), 3)

apply(pcafit$rotation[,1:2], 2, function(pc) sum(pc^2))
apply(prlin$loadings, 2, function(pc) sum(pc^2))

prloads1 <- apply(prlin$loadings, 2, function(pc){
  pc/sqrt(sum(pc^2))
})
head(round(prloads1, 3), 3)

prord <- princals(stpg)
prord

plot(prord, plot.type = "transplot", var.subset = c(1:2, 8:9), lwd = 2)

head(round(prord$objectscores, 3), 3)
head(stpg[,1:5])
head(round(prord$scoremat[,1:5], 3))

plot(prord, main = "ASTI Loadings Plot")
plot(prord, plot.type = "screeplot")


si <- rowSums(ASTI[ ,c(10, 19, 20, 21)])    
pm <- rowSums(ASTI[ ,c(1, 5, 9, 22)])    
na <- rowSums(ASTI[ ,c(3, 6, 8, 12)])     
asti2 <- data.frame(stpg, si, pm, na)

knotsord <- knotsGifi(asti2[,1:13], type = "D")
knotslin <- knotsGifi(asti2[,14:16], type = "E")
knotslist <- c(knotsord, knotslin)

prordlin <- princals(asti2, knots = knotslist, degrees = 1, ndim = 3)
colvec <- c(rep("gray", 13), rep("coral", 3))
plot(prordlin, col.loadings = colvec, plot.dim = c(1, 2))


## ------ homals 
library("MPsychoR")
library("Gifi")
data("WilPat")
WP6 <- WilPat[,c(32, 38, 41, 44, 45, 46, 47)]

homwp <- homals(WP6)
homwp

colvec <- c(rep("gray", 6), "coral4")
plot(homwp, col.points = colvec, main = "Wilson-Patterson Joint Plot")

WPmixed <- WilPat[,c(32, 38, 41, 44, 45, 46, 47:51)]
WPmixed$LibCons <- cut(WPmixed$LibCons, breaks = c(0,2,4,6,9), labels = 1:4)
WPmixed <- na.omit(WPmixed)
itknots <- knotsGifi(WPmixed[,1:6], "D")      ## item knots (data)
cknots <- knotsGifi(WPmixed[,7], "D")         ## country knots (data)
lcknots <- knotsGifi(WPmixed[,8], "D")        ## lib-cons knots (data)
lrknots <- knotsGifi(WPmixed[,9], "Q", n = 2) ## left-right (terciles)
genknots <- knotsGifi(WPmixed[,10], "D")      ## gender knots (data)
ageknots <- knotsGifi(WPmixed[,11], "E")      ## age knots (empty)
knotlist <- c(itknots, cknots, lcknots, lrknots, genknots, ageknots)

ordvec <- c(rep(FALSE, 6), FALSE, TRUE, TRUE, FALSE, TRUE)
degvec <- c(rep(-1, 7), 1, 2, -1, 1)

hommix <- homals(WPmixed, knots = knotlist, ordinal = ordvec, degrees = degvec)
plot(hommix, "transplot", var.subset = 6:11)

prinwp1 <- princals(WP6, ordinal = FALSE)

library("colorspace")
op <- par(mfrow = c(2,1))
plot(prinwp1, main = "Nominal Princals Loadings")
class(prinwp1) <- "homals"
for (i in 1:ncol(WP6)) {
  quants <- prinwp1$quantifications[[i]]
  ind <- apply(quants, 1, function(xx) all(xx == 0))
  quants <- quants[!ind, ]
  rownames(quants) <- rownames(homwp$quantifications[[i]])
  prinwp1$quantifications[[i]] <- quants
}
colvec <- rainbow_hcl(ncol(WP6), 80)
plot(prinwp1, col.points = colvec, main = "Nominal Princals Joint Plot")
for (i in 1:ncol(WP6)) {
  xy <- prinwp1$loadings[i,]
  fit <- lm(xy[2] ~ -1 + xy[1])
  abline(fit, col = colvec[i])
}
par(op)

prinwp2 <- princals(WP6, ordinal = FALSE, copies = 2)
prinwp2

copvec <- c(rep(2, 7), 1, 1, 2, 1)
prinmix <- princals(WPmixed, knots = knotlist, ordinal = ordvec, degrees = degvec, copies = copvec)

op <- par(mfrow = c(3,2), mar = c(4,4,3,2))
## marijuana
plotvars <- as.matrix(prinmix$datanum[, 6])
xlabels <- as.data.frame(prinmix$data[, 6])
x1 <- as.matrix(prinmix$datanum[, 6])
y1 <- as.matrix(prinmix$transform[, c(11, 12)])
xy <- cbind(x1, y1)
ord <- order(xy[, 1])
sfun0 <- stepfun(xy[ord, 1][-1], xy[ord, 2], f = 0)
plot(sfun0, xlab = "original", ylab = "transformed", main = "LegalizedMarijuana", xaxt = "n", col = "black", 
     do.points = FALSE, verticals = FALSE, ylim = range(xy[,2:3]))
axis(1, labels = unique(xlabels[,1]), at = unique(x1))
sfun1 <- stepfun(xy[ord, 1][-1], xy[ord, 3], f = 0)
plot(sfun1, add = TRUE, col = "red", do.points = FALSE, verticals = FALSE)
## country
plotvars <- as.matrix(prinmix$datanum[, 7])
xlabels <- as.data.frame(prinmix$data[, 7])
x1 <- as.matrix(prinmix$datanum[, 7])
y1 <- as.matrix(prinmix$transform[, 13])
xy <- cbind(x1, y1)
ord <- order(xy[, 1])
sfun0 <- stepfun(xy[ord, 1][-1], xy[ord, 2], f = 0)
plot(sfun0, xlab = "original", ylab = "transformed", main = "Country", xaxt = "n", col = "black", 
     do.points = FALSE, verticals = FALSE)
axis(1, labels = unique(xlabels[,1]), at = unique(x1))
## libcons
plotvars <- as.matrix(prinmix$datanum[, 8])
xlabels <- as.data.frame(prinmix$data[, 8])
x1 <- as.matrix(prinmix$datanum[, 8])
y1 <- as.matrix(prinmix$transform[, 14])
xy <- cbind(x1, y1)
ord <- order(xy[, 1])
sfun0 <- stepfun(xy[ord, 1][-1], xy[ord, 2], f = 0)
plot(sfun0, xlab = "original", ylab = "transformed", main = "LibCons", xaxt = "n", col = "black", 
     do.points = FALSE, verticals = FALSE)
axis(1, labels = unique(xlabels[,1]), at = unique(x1))
## left-right
plotvars <- as.matrix(prinmix$datanum[, 9])
xlabels <- as.data.frame(prinmix$data[, 9])
x1 <- as.matrix(prinmix$datanum[, 9])
y1 <- as.matrix(prinmix$transform[, 15])
xy <- cbind(x1, y1)
ord <- order(xy[, 1])
plot(xy[ord, 1], xy[ord, 2], type = "l", xlab = "original", ylab = "transformed", main = "LeftRight", 
     xaxt = "n", col = "black")
axis(1, labels = unique(xlabels[, 1]), at = unique(x1))
## Gender
plotvars <- as.matrix(prinmix$datanum[, 10])
xlabels <- as.data.frame(prinmix$data[, 10])
x1 <- as.matrix(prinmix$datanum[, 10])
y1 <- as.matrix(prinmix$transform[, 16])
xy <- cbind(x1, y1)
ord <- order(xy[, 1])
sfun0 <- stepfun(xy[ord, 1][-1], xy[ord, 2], f = 0)
plot(sfun0, xlab = "original", ylab = "transformed", main = "Gender", xaxt = "n", col = "black", 
     do.points = FALSE, verticals = FALSE)
axis(1, labels = unique(xlabels[,1]), at = unique(x1))
## Age
plotvars <- as.matrix(prinmix$datanum[, 11])
xlabels <- as.data.frame(prinmix$data[, 11])
x1 <- as.matrix(prinmix$datanum[, 11])
y1 <- as.matrix(prinmix$transform[, 17])
xy <- cbind(x1, y1)
ord <- order(xy[, 1])
plot(xy[ord, 1], xy[ord, 2], type = "l", xlab = "original", ylab = "transformed", main = "Age", xaxt = "n", 
     col = "black")
axis(1, labels = unique(xlabels[, 1]), at = unique(x1))
par(op)


loads <- prinmix$loadings[c(14, 15, 17), ]/10
rownames(loads) <-  c("LibCons", "LeftRight", "Age")
plot(loads, type = "p", pch = ".", xlab = "Dimension 2", ylab = "Dimension 2", main = "Princals Joint Plot", 
     col = "gray", asp = 1, ylim = c(-0.04, 0.06), xlim = c(-0.07, 0.05))
abline(h = 0, v = 0, col = "gray", lty = 2)
arrows(0, 0, loads[1, 1], loads[1, 2], length = 0.08, col = 1)
arrows(0, 0, loads[2, 1], loads[2, 2], length = 0.08, col = 1)
arrows(0, 0, loads[3, 1], loads[3, 2], length = 0.08, col = 1)
text(loads, labels = rownames(loads), pos = c(3,1,2), cex = 0.8, col = 1)
quants <- prinmix$quantifications[c(1:7, 10)]
for (i in 1:6) rownames(quants[[i]]) <- c(0,1,2)
rownames(quants[[7]]) <- c("India", "Hungary")
rownames(quants[[8]]) <- c("Female", "Male")
colvec <- rainbow_hcl(length(quants), 80)
for (i in 1:length(quants)) {
  points(quants[[i]], pch = 20, col = colvec[i], cex = 0.7)
  labs <- paste(names(quants)[i], rownames(quants[[i]]), sep = ".")
  text(quants[[i]], labels = labs, pch = 20, col = colvec[i], cex = 0.7, pos = 3)
}
lcquants <- prinmix$quantifications[[8]]
rownames(lcquants) <- 1:4
points(lcquants, pch = 20, col = "gray", cex = 0.7)
labs <- paste("LibCons", rownames(lcquants), sep = ".")
text(lcquants, labels = labs, pch = 20, col = "gray", cex = 0.7, pos = c(4,2,3,4))

newdat <- prord$scoremat
newR <- prord$rhat

## ----- lineals
library("MPsychoR")
library("aspect")
data("BSSS")
linbs <- lineals(BSSS)
linbs

library("lavaan")
BSSSnew <- linbs$scoremat
BSSS.model <- 'ES =~ Explore + Trip
               BS =~ Restless + Friends
               TAS =~ Frightning + Bungee
               DIS =~ Party + Illegal
               Risk =~ ES + BS + TAS + DIS'
cfabs <- cfa(BSSS.model, data = BSSSnew, estimator = "WLS", optim.method = "BFGS", check.gradient = FALSE)
round(fitMeasures(cfabs)["rmsea"], 3)
