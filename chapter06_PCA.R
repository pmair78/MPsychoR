## ------------------------ Chapter 6: principal component analysis -------------
library("mvtnorm")
sigma <- matrix(c(2, 0.8, 0, 0.8, 0.5, 0, 0, 0, 0.1), ncol = 3)
set.seed(123)
baguette <- rmvnorm(1000, mean = c(0,0,0), sigma = sigma)
svdb <- svd(baguette)
D <- diag(svdb$d)   ## singular values
round(D, 3)
V <- svdb$v         ## right singular vectors
round(V, 3)

v1 <- V[,1]
v2 <- V[,2]
round(v1 %*% v2, 7)  ## inner product
sqrt(sum(v1^2))      ## length 
sqrt(sum(v2^2))

library("rgl")
plot3d(baguette, col = "gray", xlim = c(-4, 4), ylim = c(-4,4), zlim = c(-4,4), aspect = 1, size = 2)
arrow3d(c(0,0,0), V[,1], col = "red")
arrow3d(c(0,0,0), V[,2], col = "red")
arrow3d(c(0,0,0), V[,3], col = "red")

U <- svdb$u
baguette1 <- U %*% D %*% t(V)
identical(round(baguette, 7), round(baguette1, 7))
covb <- cov(baguette)
eigend <- eigen(covb)
L <- diag(eigend$values)    ## eigenvalues
round(L, 3)
U <- eigend$vectors         ## eigenvectors
round(U, 3)

n <- nrow(baguette)
l <- (diag(D)/sqrt(n-1))^2  ## from singular values
round(l, 3)                 ## eigenvalues

## ----- principal component analysis
library("MPsychoR")
data("BrainIQ")
X <- BrainIQ[, c("VIQ", "PIQ", "MRI_Count")]
head(X, 4)

PCAfit <- princomp(X, cor = TRUE)
PCAfit
round(unclass(PCAfit$loadings), 3)  ## weights (loadings)
evIQ <- eigen(cor(X))               ## eigenvalue decomposition
sqrt(evIQ$values)                   ## standard deviations
round(evIQ$vectors, 3)              ## weights (loadings)
summary(PCAfit)
evIQ$values/sum(evIQ$values)

PCAfit2 <- prcomp(X, scale = TRUE)
print(PCAfit2, digits = 3)

n <- nrow(X)
svdIQ <- svd(scale(X)/sqrt(n-1))      ## SVD
round(svdIQ$d, 3)         ## singular values
round(svdIQ$v, 3)         ## right singular vectors
head(PCAfit2$x, 4)
head(svdIQ$u %*% diag(svdIQ$d)*sqrt(n-1), 4)


data("Privacy")
pcaPriv <- prcomp(Privacy, scale = TRUE)
screeplot(pcaPriv, type = "lines", main = "Privacy Scree Plot")
abline(h = 1, col = "gray", lty = 2)
round(pcaPriv$rotation[,1:3], 3)

library("psych")
pcaPriv1 <- principal(cor(Privacy), 3, rotate = "none")
pcaPriv1$loadings

pcaPriv2 <- principal(cor(Privacy), 3, rotate = "varimax")
pcaPriv2$loadings

pcaPriv3 <- principal(cor(Privacy), 3, rotate = "promax")
pcaPriv3$loadings
round(pcaPriv3$Phi, 3)

biplot(pcaPriv, col = c("gray", "black"), cex = c(0.6, 0.9))
abline(h = 0, v = 0, lty = 2)

## ----- some PCA variants
set.seed(123)
pcarob <- princomp(covmat = MASS::cov.rob(Privacy), cor = TRUE)

library("elasticnet")
spcaPriv <- spca(scale(Privacy), K = 3, sparse = "varnum", para = c(3, 4, 3))
spcaPriv


## ----- three-way principal component analysis
data("SDOwave")
SDOar <- array(unlist(SDOwave), dim = c(612, 4, 5))
dnames <- dimnames(SDOar) <- list(1:612, paste0("SDO", 1:4), 1996:2000)
dim(SDOar)

library("multiway")
SDOarc <- ncenter(SDOar, mode = 1)
dimnames(SDOarc) <- dnames

set.seed(66)
sdopara <- parafac(SDOarc, nfac = 2)
str(sdopara[1:4])
dimnames(sdopara$B) <- list(dnames[[2]], paste0("Comp.", 1:2))
round(sdopara$B, 3)  ## loadings

op <- par(mfrow = c(1,2))
plot(sdopara$A, xlab = "Dimension 1", ylab = "Dimension 2", main = "SDO Parafac Biplot", col = "gray", 
     cex = 0.5, asp = 1, ylim = c(-4, 2))
abline(h = 0, v = 0, lty = 2, col = "gray")
arrows(0, 0, sdopara$B[,1], sdopara$B[,2], length = 0.10)
text(sdopara$B, labels = dnames[[2]], cex = 0.8, pos = 3)
plot(1:5, sdopara$C[,1], type = "l", xaxt = "n", ylab = "Component Scores", main = "Parafac Occasion Components",
     xlab = "Time", col = "cadetblue", lwd = 2, ylim = c(-0.5, -0.8))
axis(1, 1:5, labels = 1996:2000)
lines(1:5, sdopara$C[, 2], type = "l", col = "coral", lwd = 2)
legend("bottomright", legend = c("Component 1", "Component 2"), lty = 1, col = c("cadetblue", "coral"))
par(op)

Xhat <- fitted(sdopara)  ## structural image
SSE <- sum((SDOar - Xhat)^2)/prod(dim(SDOar))
round(SSE, 3)

set.seed(111)
ndims <- c(3, 2, 2)  ## set P, Q, R
sdotuck <- tucker(SDOar, nfac = ndims)
dimnames(sdotuck$B) <- list(dnames[[2]], paste0("Comp.", 1:2))
round(sdotuck$B, 3)

op <- par(mfrow = c(1,2))
plot(sdotuck$A[,1:2], xlab = "Dimension 1", ylab = "Dimension 2", main = "SDO Tucker Biplot", col = "gray", cex = 0.5, asp = 1)
abline(h = 0, v = 0, lty = 2, col = "gray")
arrows(0, 0, sdotuck$B[,1]/10, sdotuck$B[,2]/10, length = 0.10)
text(sdotuck$B/10, labels = dnames[[2]], cex = 0.8, pos = 3)
plot(1:5, sdotuck$C[,1], type = "l", xaxt = "n", ylab = "Component Scores", main = "Tucker Occasion Components",
     xlab = "Time", col = "cadetblue", lwd = 2, ylim = c(-0.5, 1))
axis(1, 1:5, labels = 1996:2000)
lines(1:5, sdotuck$C[, 2], type = "l", col = "coral", lwd = 2)
legend("topright", legend = c("Component 1", "Component 2"), lty = 1, col = c("cadetblue", "coral"))
par(op)

round(sdotuck$G, 2)

## ----- independent component analysis
library("eegkit")
data("storcap")
elecvec <- c("FP1","FP2", "F3", "F4", "FC5", "FC6", "T7", "T8", "FC1", "FC2", "C3", "C4", 
             "P7", "P8", "P3", "P4", "CP5", "CP6", "CP1", "CP2", "PO7", "PO8", "PO3", "PO4", "O1", "O2")
eegcap(elecvec, type = "2d", col.point = "gray", col.label = "black", cex.label = 0.8, cex.point = 4)
data(eeghead)
shade3d(eeghead)
eeghead$material$color <- rep(1,length(eeghead$material$color))
wire3d(eeghead)
eegcap(elecvec, type = "3d", col.point = "coral4", cex.point = 0.3, col.label = "coral4", 
       head = FALSE, add = TRUE)

eegcap(elecvec, type = "2d", col.point = "gray", col.label = "black", cex.label = 0.8, cex.point = 4)

library("plyr")
PO34 <- subset(storcap, channel == "PO3_4")    ## select PO3/4
PO34agg <- ddply(PO34, .(time, cond), summarize, mean = mean(voltage), sd = sd(voltage))

library("colorspace")
ylims <- range(PO34agg$mean)
cols <- rainbow_hcl(4, 80)
eegtime(PO34agg$time[PO34agg$cond == "SS1_Ips"], xlab = "Time", voltage = PO34agg$mean[PO34agg$cond == "SS1_Ips"], 
        main = "Electrode PO3/4", vcol = cols[1], ylim = ylims)
eegtime(PO34agg$time[PO34agg$cond == "SS1_Contra"], voltage = PO34agg$mean[PO34agg$cond == "SS1_Contra"], 
        add = TRUE, vcol = cols[2])
eegtime(PO34agg$time[PO34agg$cond == "SS3_Ips"], voltage = PO34agg$mean[PO34agg$cond == "SS3_Ips"], 
        add = TRUE, vcol = cols[3])
eegtime(PO34agg$time[PO34agg$cond == "SS3_Contra"], voltage = PO34agg$mean[PO34agg$cond == "SS3_Contra"], 
        add = TRUE, vcol = cols[4])
abline(v = c(300, 1200), lty = 2, col = "gray")
legend("bottomright", c("Set Size 1 (Ipsilateral)", "Set Size 1 (Contralateral)", "Set Size 3 (Ipsilateral)", 
                        "Set Size 3 (Contralateral)"), lty = 1, col = cols, bty = "n")

storcap1 <- storcap
SS1_CDA <- storcap1[storcap1$cond == "SS1_Contra",]$voltage - storcap1[storcap1$cond == "SS1_Ips",]$voltage
SS3_CDA <- storcap1[storcap1$cond == "SS3_Contra",]$voltage - storcap1[storcap1$cond == "SS3_Ips",]$voltage
n <- length(c(SS1_CDA, SS3_CDA))/2
storcap1 <- storcap1[1:(2*n),]
storcap1$voltage <- c(SS1_CDA, SS3_CDA)
storcap1$cond <- factor(rep(c("SS1_CDA", "SS3_CDA"), each = n))

storcap700 <- subset(storcap1, time == 700)
CDA700agg <- ddply(storcap700, .(channel, cond), summarize, mean = mean(voltage)) 
CDA700agg1 <- CDA700agg2 <- CDA700agg
even <- seq(2, nrow(CDA700agg), 2)    ## even row selector
odd <- seq(1, nrow(CDA700agg), 2)     ## odd row selector
CDA700agg1[even, ] <- CDA700agg1[odd, ]  
CDA700agg2[odd, ] <- CDA700agg2[even, ]  
data(eegcoord)        ## coordinate template from eegkit
coords <- eegcoord[elecvec, 1:3]      ## 3D coordinates

eegspace(coords, CDA700agg1[,3], main = "Set Size 1 CDA", vlim = range(CDA700agg[,3]), colorlab = "")
eegspace(coords, CDA700agg2[,3], main = "Set Size 3 CDA", vlim = range(CDA700agg[,3]), colorlab = "")

CDA1 <- subset(storcap1, cond == "SS1_CDA")
CDA1agg <- daply(CDA1, .(channel, time), function(x) mean(x$voltage)) 
tempICA1 <- eegica(CDA1agg, nc = 4, type = "time")   
CDA3 <- subset(storcap1, cond == "SS3_CDA")
CDA3agg <- daply(CDA3, .(channel, time), function(x) mean(x$voltage)) 
tempICA3 <- eegica(CDA3agg, nc = 4, type = "time")   
round(tempICA1$vafs, 3)
round(tempICA3$vafs, 3)

tvec <- unique(storcap1$time)
mixmat1 <- tempICA1$M[rep(1:nrow(tempICA1$M), each = 2), ]  
mixmat3 <- tempICA3$M[rep(1:nrow(tempICA3$M), each = 2), ] 
vlims <- range(c(mixmat1[,1], mixmat1[,2], mixmat3[,1]))
op <- par(mfrow = c(3,2))
eegtime(tvec, tempICA1$S[, 1], main = "Component 1 (Set Size 1 CDA)", vcol = 1)
eegspace(coords[, 1:2], mixmat1[, 1], vlim = vlims)
eegtime(tvec, tempICA1$S[, 2], main = "Component 2 (Set Size 1 CDA)", vcol = 1)
eegspace(coords[, 1:2], mixmat1[, 2], vlim = vlims)
eegtime(tvec, tempICA3$S[, 1], main = "Component 1 (Set Size 3 CDA)", vcol = 1)
eegspace(coords[, 1:2], mixmat3[, 1], vlim = vlims)
par(op)

