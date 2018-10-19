## ----------------------- Chapter 10: Biplots ----------------------------

## ----- variable and subject space
library("MPsychoR")
Enorm <- function(x) sqrt(sum(x^2))     ## computes Euclidean norm

data("BrainIQ")   
x <- as.vector(scale(BrainIQ$VIQ, scale = FALSE))   ## center x
y <- as.vector(scale(BrainIQ$PIQ, scale = FALSE))   ## center y
xn <- Enorm(x)      ## length x
yn <- Enorm(y)      ## length y
theta <- acos(crossprod(x, y)/(Enorm(x)*Enorm(y)))  ## angle (in radians)

op <- par(mfrow = c(1,2))
plot(BrainIQ$VIQ, BrainIQ$PIQ, pch = 19, xlab = "VIQ", ylab = "PIQ", main = "IQ Variable Space")
plot(0:160, 0:160, ylim = c(0, 110), type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, 
     main = "IQ subject space")  ## create empty plot panel
arrows(0, 0, xn, 0, length = 0.10, lwd = 2)              ## x vector 
text(xn, 0, "VIQ", pos = 4, offset = 0.4)
y1 <- yn * cos(theta)                           ## x-coordinate y vector (orthogonal y projection on x)
y2 <- yn * sin(theta)                           ## y-coordinate y vector (orthogonal y projection on vector orthogonal to x)
arrows(0, 0, y1, y2, length = 0.10, lwd = 2)             ## y vector 
text(y1, y2, "PIQ", pos = 4, offset = 0.4)
par(op)

## ----- regression biplots
library("MPsychoR")
data("BrainIQ")
BrainIQ <- na.omit(BrainIQ[,-1])  ## we omit NAs and gender
rownames(BrainIQ) <- 1:nrow(BrainIQ)  ## relabel persons
BrainIQ1 <- as.data.frame(scale(BrainIQ)) ## standardize 

regfit <- lm(cbind(VIQ, PIQ, Weight) ~ -1 + Height + MRI_Count, data = BrainIQ1)
colnames(regfit$coef) <- c("VIQ", "PIQ", "Weight")
round(regfit$coef, 3)   ## vector coordinates
R2vec <- sapply(summary(regfit), `[[`, "r.squared")

plot(BrainIQ1$Height, BrainIQ1$MRI_Count, pch = 20, col = "gray", xlab = "Body Height (standardized)", 
     ylab = "Brain Size (standardized)", main = "Regression Biplot", asp = 1)
abline(h = 0, v = 0, lty = 2, col = "darkgray")
arrows(0, 0, regfit$coef[1,1], regfit$coef[2,1], length = 0.1, col = "brown")
text(regfit$coef[1,1], regfit$coef[2,1], labels = "VIQ", col = "brown", pos = 2, cex = 0.9)
arrows(0, 0, regfit$coef[1,2], regfit$coef[2,2], length = 0.1, col = "brown")
text(regfit$coef[1,2], regfit$coef[2,2], labels = "PIQ", col = "brown", pos = 3, cex = 0.9)
arrows(0, 0, regfit$coef[1,3], regfit$coef[2,3], length = 0.1, col = "brown")
text(regfit$coef[1,3], regfit$coef[2,3], labels = "Weight", col = "brown", pos = 4, cex = 0.9)

round(cor(BrainIQ[, 2:4]), 3)

library("calibrate")
plot(BrainIQ1$Height, BrainIQ1$MRI_Count, pch = 20, cex = 0.8, xlab = "Height", ylab = "MRI", 
     col = "darkblue", asp = 1, main = "Orthogonal Projections")
text(BrainIQ1$Height, BrainIQ1$MRI_Count, labels = 1:nrow(BrainIQ1), cex = 0.7, pos = 3, col = "darkblue")
abline(h = 0, v = 0, lty = 2, col = "darkgray")
calibrate.Z <- calibrate(regfit$coef[,1], BrainIQ1$VIQ, seq(-2,2, by = 0.5), 
                         cbind(BrainIQ1$Height, BrainIQ1$MRI_Count), dp = TRUE, axiscol = "brown",
                         axislab = "VIQ", labpos = 3, verb = FALSE)

VIQcal <- calibrate.Z$yt*sd(BrainIQ$VIQ) + mean(BrainIQ$VIQ)

summary(BrainIQ[, 2])   ## basic location measures
BrainIQ[c(20, 11), -c(1, 4)]  ## biplot projections
round(R2vec, 3)

## ------ principal component biplot
library("MPsychoR")
data("BrainIQ")
BrainIQ1 <- na.omit(BrainIQ[,-1])
head(BrainIQ1, 3)

pca_biq1 <- prcomp(BrainIQ1)
pca_biq2 <- prcomp(BrainIQ1, scale = TRUE)

op <- par(mfrow = c(1,2))
biplot(pca_biq1, pc.biplot = TRUE, cex = c(0.6, 0.8), col = c("gray", "coral1"), arrow.len = 0.05, 
       main = "Biplot (Unstandardized)", xlim = c(-4, 4), asp = 1, cex.axis = 0.8)
biplot(pca_biq2, pc.biplot = TRUE, cex = c(0.6, 0.8), col = c("gray", "coral1"), arrow.len = 0.05, 
       main = "Biplot (Standardized)", asp = 1, cex.axis = 0.8)
layout(1)
par(op)

data("yaass")
pca_yaass1 <- prcomp(yaass[,1:5])
pca_yaass2 <- prcomp(yaass[,1:5], scale = TRUE)

op <- par(mfrow = c(2,2), mar = c(4, 3, 4, 1))
biplot(pca_yaass1, pc.biplot = TRUE, cex = c(0.6, 0.8), col = c("gray", "coral1"), arrow.len = 0.05, 
       main = expression(paste("Biplot (Unstandardized, ", alpha, "=0)")), scale = 0, cex.axis = 0.8)
abline(h = 0, v = 0, col = "gray", lty = 2)
biplot(pca_yaass1, pc.biplot = TRUE, cex = c(0.6, 0.8), col = c("gray", "coral1"), arrow.len = 0.05, 
       main = expression(paste("Biplot (Unstandardized, ", alpha, "=1)")), scale = 1, cex.axis = 0.8)
abline(h = 0, v = 0, col = "gray", lty = 2)
biplot(pca_yaass2, pc.biplot = TRUE, cex = c(0.6, 0.8), col = c("gray", "coral1"), arrow.len = 0.05, 
       main = expression(paste("Biplot (Standardized, ", alpha, "=0)")), scale = 0, cex.axis = 0.8)
abline(h = 0, v = 0, col = "gray", lty = 2)
biplot(pca_yaass2, pc.biplot = TRUE, cex = c(0.6, 0.8), col = c("gray", "coral1"), arrow.len = 0.05, 
       main = expression(paste("Biplot (Standardized, ", alpha, "=1)")), scale = 1, cex.axis = 0.8)
abline(h = 0, v = 0, col = "gray", lty = 2)
par(op)

round(cor(yaass[,1:5]), 2)
round(apply(yaass[,1:5], 2, sd), 3)

X <- pca_yaass2$x[, 1:2]  ## extract PC scores
Y <- scale(yaass[,1:5])   ## standardize variables
fitlms <- lm(Y ~ -1 + X)  ## fit regressions
round(coef(fitlms), 3)
round(t(pca_yaass2$rotation[,1:2]), 3)

plot(X[,1], X[,2], pch = 20, xlab = "PC1", ylab = "PC2", col = "darkblue", asp = 1, main = "Biplot Axis", 
     xlim = c(-3.2, 3.2))
text(X[,1], X[,2], labels = rownames(X), cex = 0.7, pos = 3, col = "darkblue")
abline(h = 0, v = 0, lty = 2, col = "gray")
calAE <- calibrate(fitlms$coef[,"AE"], Y[,"AE"], tm = seq(-2, 2, by = 0.5), Fr = X, dp = TRUE, 
                   axiscol = "brown", axislab = "AE", labpos = 3, verb = FALSE)

R2vec <- sapply(summary(fitlms), "[[", "r.squared")
round(R2vec[2], 3)

library("bpca")
resbi <- bpca(yaass[, 1:5], scale = FALSE, method = c("gh"))  
colvec <- c("cadetblue", "chartreuse4")[unclass(yaass$Group)]
plot(resbi, main = "YAASS Biplot", obj.color = colvec, var.factor = 3, obj.cex = 0.8, asp = 1, 
     xlim = c(-2, 2), ylim = c(-2, 2))
legend("topleft", legend = c("high risk psychosis", "healthy controls"), pch = 20, col = unique(colvec))

library("Gifi")
ABC6 <- ABC[,6:11] 
fitabc <- princals(ABC6)

plot(fitabc, plot.type = "biplot", main = "Princals Biplot", expand = 0.7, cex.scores = 0.6, col.scores = "gray")
abline(h = 0, v = 0, lty = 2, col = "gray")

## ----- MDS biplot
library("MPsychoR")
library("smacof")
data("NeuralActivity")
delta <- Reduce("+", NeuralActivity)/length(NeuralActivity)
fit_neural <- mds(delta, type = "interval")
fit_neural

data("NeuralScales")
mdsbi <- biplotmds(fit_neural, NeuralScales)
plot(mdsbi, main = "Neural Activity MDS Biplot", col = "gray", label.conf = list(col = "gray"), 
     vec.conf = list(col = "brown", cex = 1))

X <- fit_neural$conf
Y <- scale(NeuralScales, scale = TRUE)
plot(X, pch = 20, cex = 0.6, xlab = "Dimension 1", ylab = "Dimension 2", col = "darkblue", asp = 1, 
     main = "Biplot MDS Emotion Axis")
text(X, labels = rownames(X), cex = 0.7, pos = 3, col = "darkblue")
abline(h = 0, v = 0, lty = 2, col = "gray")
calEm <- calibrate(mdsbi$coef[,"Emotion"], Y[,"Emotion"], tm = seq(-2, 1.5, by = 0.5), Fr = X, 
                   dp = TRUE, axiscol = "brown", axislab = "Emotion", labpos = 3, verb = FALSE)



## ----- correspondence analysis biplots
library("anacor")
superfan <- as.table(matrix(c(9, 12, 8, 1, 13, 1, 6, 20, 15, 4, 23, 18), ncol = 3))
attr(superfan, "dimnames") <- list(c("Slayer", "Iron Maiden", "Metallica", "Judas Priest"), 
                                   c("Horst", "Helga", "Klaus"))
fit_fans <- anacor(superfan, scaling = c("Benzecri", "standard"))
plot(fit_fans, main = "Asymmetric Superfan CA Map")

Srole <- as.table(matrix(c(64, 94, 58, 46, 
                           57, 94, 54, 40,
                           57, 105, 65, 60,
                           72, 141, 77, 94,
                           36, 97, 54, 78,
                           21, 71, 54, 71), nrow = 4))
attr(Srole, "dimnames") <- list(mhealth = c("well", "mild", "moderate", "impaired"), ses = LETTERS[1:6])
Srole

fit_ses <- anacor(Srole, scaling = c("standard", "Benzecri"))
plot(fit_ses, arrows = c(T, F), main = "Asymmetric CA Map")

