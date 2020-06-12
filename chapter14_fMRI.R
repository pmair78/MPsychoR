## ------------------ Chapter 14: Analysis of fMRI Data ------------------

## ----- fMRI data manipulation
library("fmri")
scanS1R1 <- read.NIFTI("s01_r01_FUNC.nii") 
summary(scanS1R1)
imageS1R1 <- fmri::extractData(scanS1R1)

imdim <- dim(imageS1R1)
rm(scanS1R1)
imageS1R1[50, 50, 10, 100]

mask <- read.ANALYZE("globalmask.img")
mask <- fmri::extractData(mask)[,,,1]
imageS1R1r <- read.NIFTI("wFUN1.nii")
imageS1R1 <- fmri::extractData(imageS1R1r)
dim(imageS1R1)

imageS1R1T1 <- imageS1R1[,,,1]

library("AnalyzeFMRI")
sigma <- 2
unsmooth <- GaussSmoothArray(imageS1R1T1, sigma = diag(0, 3), mask = mask)   
smooth <- GaussSmoothArray(imageS1R1T1, sigma = diag(sigma^2, 3), mask = mask)  

library("ggBrain")   ## devtools::install_github('aaronjfisher/ggBrain',build_vignettes=TRUE)
plotsmooth <- ggBrain(brains = list(unsmooth, smooth), mask = mask, mar = c(3,3), mar_ind = c(50,50), 
                      brain_ind = c(1,2), col_ind = c("Unsmoothed", "Smoothed"), type = 'signed')
plotsmooth

library("RNiftyReg")  
template <- RNiftyRef::readNifti(system.file("MNI152_T1_2mm_brain.nii.gz", package = "brainR"))
dim(template)

imageS1R1a <- RNiftyReg::readNifti("wFUN1.nii")[,,,10]
dim(imageS1R1a)

template1 <- niftyreg(source = template, target = imageS1R1a, scope = "affine")$image 
template1[is.na(template1)] <- 0
dim(template1)

library("brainR")
contour3d(template1, level = 3500, alpha = 0.1, draw = TRUE)
contour3d(imageS1R1a, level = c(1000, 2000), add = TRUE, alpha = c(0.8, 0.9), 
          color = c("yellow", "red"), mask = mask)

voxelts <- imageS1R1[50, 50, 10,]
plot(1:length(voxelts), voxelts, type = "l", main = "Voxel Time Series", ylab = "BOLD Intensity", xlab = "Time")

dim(imageS1R1)                      ## 4D represenation
nt <- dim(imageS1R1)[4]             ## number of scans
imFlat <- fourDto2D(imageS1R1, nt)  ## flatten 4D input data
dim(imFlat)                         ## all voxels
maskind <- which(mask > 0)          
imFlatm <- imFlat[, maskind]        
dim(imFlatm)                        ## non-background voxels

im4D <- twoDto4D(imFlat, dim = dim(imageS1R1))
dim(im4D)

library("brainGraph")
head(aal90)

ROI <- factor(aal90$name)    
mnixyz <- as.data.frame(aal90[,2:4, with = FALSE])

dims <- dim(imageS1R1)[1:3]  ## spatial dimensions
Cmat <- cbind(expand.grid(z = 1:dims[1], y = 1:dims[2], x = 1:dims[3])[,3:1], 1)

Tmat <- rbind(imageS1R1r$header$srowx, imageS1R1r$header$srowy, imageS1R1r$header$srowz, c(0, 0, 0, 1))
Tmat

mni <- Tmat %*% t(Cmat)
mni <- (t(mni))[,-4]

maskTF <- as.logical(mask)  ## convert into boolean mask
mni <- mni[maskTF, ]        ## apply mask 

library("class")
fitknn <- class::knn(mnixyz, mni, ROI, k = 3)
length(fitknn)

imageS1R1 <- GaussSmoothArray(imageS1R1, sigma = diag(sigma^2, 3), mask = mask)
imgMat <- fourDto2D(imageS1R1, dim(imageS1R1)[4])
imgMat <- imgMat[,maskTF]        ## apply mask
ROIMat <- t(apply(imgMat, 1, function(vox) tapply(vox, fitknn, mean)))
dim(ROIMat)



## ----- linear modeling of fMRI data
library("MPsychoR")
data("NeuralScanner")
onsetsS1R1 <- NeuralScanner$TRIAL_START[NeuralScanner[,2] == 1] 
durationS1R1 <- NeuralScanner$RT[NeuralScanner[,2] == 1]  + 1 
nt <- 162

library("fmri")
xt <- fmri.stimulus(scans = nt, onsets = onsetsS1R1, dur = durationS1R1, TR = 2.5, times = TRUE)
plot(1:nt, as.vector(xt), type = "l", main = "Expected BOLD Time Series", ylab = "Expected BOLD Response", 
     xlab = "Time")

X <- fmri.design(xt, order = 2)
head(X)
dim(X)

library("fields")
col5 <- colorRampPalette(c('cadetblue4', 'white', 'coral4'))
collev <- 21
maxval <- max(abs(X))
colbreak <- seq(-maxval, maxval, length.out = collev + 1)
dev.new()
op <- par(mar = c(5,5,5,7))
image(t(X), axes = FALSE, main = "Simple fMRI Design Matrix", xlab = "Contrasts", 
      col = col5(n = collev), breaks = colbreak)
axis(1, at = c(0, 0.33, 0.66, 1), labels = c("E(BOLD)", "Intercept", "Linear", "Quadratic"))
box()
image.plot(t(X), legend.only = TRUE , col = col5(n = collev), breaks = colbreak)
par(op)

imageS1R1 <- read.NIFTI("wFUN1.nii")
dim(fmri::extractData(imageS1R1))
spmS1R1 <- fmri.lm(imageS1R1, X, verbose = TRUE)
dim(spmS1R1$beta) 

spmsmoothS1R1 <- fmri.smooth(spmS1R1)

library("MPsychoR")
data("NeuralScanner")
nruns <- max(NeuralScanner[,2])   ## number of runs (16)
ncond <- max(NeuralScanner[,3])   ## number of conditions (60)
nmeas <- nrow(NeuralScanner)      ## number of measurents (960)
nt <- 162                         ## number of scans per run
durRun <- 400                     ## 1 run is 400sec

add <- rep(seq(0, 6000, by = durRun), each = ncond)
TrStart <- NeuralScanner$TRIAL_START + add              

onsetsMat <- do.call(rbind, split(TrStart, NeuralScanner$STATE_INDEX))    
durationMat <- do.call(rbind, split(NeuralScanner$RT + 1, NeuralScanner$STATE_INDEX))  

library("fmri")
Xcond <- matrix(0, nt*nruns, ncond)
for (i in 1:ncond) {
  Xcond[,i] <- fmri.stimulus(scans = nt*nruns, onsets = onsetsMat[i,], dur = durationMat[i,], 
                             TR = 2.5, times = TRUE)
}

ind <- rep(1, nt)           ## interept
lin <- (1:nt)/nt            ## linear trend
quad <- (0.5 - lin)^2       ## quadratic trend
Xtrend1 <- cbind(ind, lin, quad)

blow <- diag(1, nruns)
Xtrend <- blow %x% Xtrend1     ## expand to 16 runs

library("magic")
data("NeuralHM")
Xhm <- NeuralHM[[1]]
for (i in 2:length(NeuralHM)) {
  Xhm <- adiag(Xhm, NeuralHM[[i]])  
}
X <- cbind(Xcond, Xtrend, Xhm)

col5 <- colorRampPalette(c('deepskyblue4', 'lightgray', 'darkred'))  
collev <- 41
maxval <- max(abs(X))
colbreak <- seq(-maxval, maxval, length.out = collev + 1)
op <- par(mar=c(5,5,5,7))
image(t(X)[, ncol(t(X)):1], axes = FALSE, main = "Design Matrix", col = col5(n = collev), breaks = colbreak)
axis(1, at = c(0.15, 0.4, 0.75), labels = c("Conditions", "Trend", "Head Motion"), tick = FALSE)
axis(1, at = c(0, ncol(Xcond)/ncol(X), (ncol(Xcond) + ncol(Xtrend))/ncol(X), 1), tick = TRUE, labels = FALSE)
box()
image.plot(t(X), legend.only=TRUE , col = col5(n = collev), breaks = colbreak)
par(op)

nruns <- 3
X3 <- X[1:(nt*nruns), ]
X3 <- X3[, colSums(X3) != 0]

mask <- read.ANALYZE("globalmask.img")
mask <- fmri::extractData(mask)[,,,1]

scanfiles <- paste0("wFUN", 1:nruns, ".nii")
imageS1 <- array(0, dim = c(dim(mask), nt*nruns))
sigma <- 2
for (i in 1:nruns) {
  scanRi <- read.NIFTI(scanfiles[i])
  scanRi <- fmri::extractData(scanRi)
  scanRi <- GaussSmoothArray(scanRi, sigma = diag(sigma^2,3), mask = mask)
  start4d <- (i-1)*nt + 1
  end4d <- i*nt
  imageS1[,,,start4d:end4d] <- scanRi
}

f.write.nifti(imageS1, "imageS1", nii = TRUE)
imageS1 <- read.NIFTI("imageS1")

maskTF <- twoDto4D(as.logical(mask), dim = dim(mask))
imageS1$mask <- maskTF

spm3 <- fmri.lm(imageS1, X3, actype = "noac")
str(spm3$beta) 
spm3_2D <- t(fourDto2D(spm3$beta, dim(spm3$beta)[4]))[,1:60]
dim(spm3_2D)


## ----- multiple comparisons in fMRI
mask <- read.ANALYZE("globalmask.img")
mask <- fmri::extractData(mask)[,,,1]
imgdim <- dim(mask)
n <- 20
sigma <- 2
files <- paste0("con_s", sprintf("%02d", 1:20), "_0002")
spmC1 <- numeric()
for (i in 1:n) {
  spmi <- read.ANALYZE(files[i])
  spmi <- extractData(spmi)[,,,1]
  spmiSmooth <- GaussSmoothArray(spmi, sigma = diag(sigma^2,3), mask = mask)
  spmiSmooth <- spmiSmooth[mask != 0]
  spmC1 <- rbind(spmC1, spmiSmooth)
}

library("genefilter")     ## BiocManager::install("genefilter")
fitT <- colttests(spmC1, gl(1, nrow(spmC1)))
sum(fitT$p.value <= 0.05)

nvox <- length(fitT$statistic)
alphanew <- 0.05/nvox
sum(fitT$p.value <= alphanew)

pFDR <- p.adjust(fitT$p.value, method = "fdr")
sum(pFDR <= 0.05)

sigma <- 2
covmat <- diag(sigma^2,3)
trf <- Threshold.RF(p.val = 0.05, sigma = covmat, num.vox = nvox, type = "t", df = nvox-1)
trf
sum(fitT$statistic >= trf)

set.seed(123)
nperm <- 1000
maxtvec <- numeric(nperm)
for (i in 1:nperm) {
   print(i)
   signperm <- sample(c(-1, 1), n, replace = TRUE)
   permmat <- spmC1*signperm
   maxtvec[i] <- max(colttests(permmat, gl(1, nrow(permmat)))$statistic)
}

tcrit <- quantile(maxtvec, probs = 0.95)
tcrit
sum(fitT$statistic >= tcrit)

pthresh <- 0.001
p01obs <- array(ifelse(fitT$p.value <= pthresh, 1, 0), imgdim)

library("mmand")
kwidth <- 5
kernel <- shapeKernel(kwidth, dim = 3)
compsobs <- components(p01obs, kernel)
clustsize <- table(c(compsobs), useNA = "no")

set.seed(111)
nperm <- 1000
maxS <- numeric(nperm)
for (i in 1:nperm) {
  print(i)
  signperm <- sample(c(-1, 1), n, replace = TRUE)
  permmat <- spmC1*signperm
  pvec <- colttests(permmat, gl(1, nrow(permmat)))$p.value
  p01perm <- array(ifelse(pvec <= pthresh, 1, 0), imgdim)
  comps <- components(p01perm, kernel)
  maxS[i] <- suppressWarnings(max(table(c(comps), useNA = "no")))
}

Scrit <- quantile(maxS, probs = 0.95)
Scrit
sum(clustsize >= Scrit)

op <- par(mfrow = c(2,1))
hist(maxtvec, freq = FALSE, breaks = 20, xlab = "max-T values", main = "Max-T Permutation Distribution")
abline(v = tcrit, lty = 2)
hist(maxS, freq = FALSE, breaks = 20, xlab = "max-S values", main = "Max-S Permutation Distribution")
abline(v = Scrit, lty = 2)
par(op)


## ----- independent component analysis in fMRI
library("AnalyzeFMRI")
f.write.nifti(imageS1, "imageS1", nii = TRUE)

set.seed(123)
f.icast.fmri("imageS1.nii", "globalmask.img", is.spatial = TRUE)

tsICAs <- as.matrix(read.table("imageS1-ICAs-time-series.dat", header = TRUE, row.names = 1))
str(tsICAs)

library("ggplot2")
library("cowplot")
imgICAS1 <- read.NIFTI("imageS1_ICAs.nii")
imgICAS1 <- extractData(imgICAS1)[,,,1]
tsdf <- data.frame(time = 1:486, IC1 = tsICAs[,1])
bp <- ggBrain(brains = imgICAS1, mask = mask, mar = 3, mar_ind = 40, type = 'signed') + theme_gray()
tsp <- ggplot(tsdf, aes(time, IC1)) + geom_line() + theme_gray()
plot_grid(bp, tsp, labels=c("Activation", "Time Series"), ncol = 1, nrow = 2)

## ----- representational similarity analysis
library("smacof")
R <- cor(spm3_2D)
RDM <- sim2diss(R)

library("MPsychoR")
data("NeuralActivity")
Delta <- Reduce("+", NeuralActivity)/n

RSAfit <- mds(Delta, ndim = 3, type = "interval")
RSAfit

library("plot3D")
col5 <- colorRampPalette(c('cadetblue', 'coral'))
text3D(RSAfit$conf[,1], RSAfit$conf[,2], RSAfit$conf[,3], labels = row.names(RSAfit$conf),
       bty = "g", ticktype = "detailed", d = 2, adj = 0.5, font = 2, cex = 0.7, xlab = "D1", ylab = "D2", zlab = "D3",
       col = col5(21), theta = 60, phi = 20, colvar = RSAfit$spp, main = "3D RSA Representation", clab = "SPP")


## ----- functional connectivity analysis
library("brainGraph")
library("brainwaver")
data(brain)
colnames(brain) <- aal90$name
dim(brain)

PreCG.L <- brain[,1]  ## seed time series

library("psych")
ROIs <- brain[,-1]    ## remaining ROI time series
corvec <- cor(PreCG.L, ROIs)  ## correlation
zvec <- fisherz(corvec)       ## z-transformation

SCAdf <- data.frame(x = aal90$x, y = aal90$y, names = aal90$name, z = c(NA, zvec))
sca <- ggplot(SCAdf, aes(x = x, y = y, label = names, color = z)) + 
  scale_color_gradient(low = "white", high = "cadetblue", na.value = "black")
sca + geom_text(size = 3) + ggtitle("Seed-Based Correlations") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_blank(),  axis.ticks.y = element_blank())

wavecor <- const.cor.list(as.matrix(brain), n.levels = 6)

R <- 0.4
pvalues <- p.value.compute(wavecor[[4]], proc.length = nrow(brain), sup = R, num.levels = 4)
FDRthresh <- compute.FDR(pvalues, q = 0.05)
p01vec <- ifelse(pvalues <= FDRthresh, 1, 0)  
p01 <- diag(0, ncol(brain))
p01[lower.tri(p01)] <- p01vec
p01 <- p01 + t(p01)    
diag(p01) <- 1

op <- par(mfrow = c(1,2), mar = c(2, 3, 4, 2))
image(wavecor[[4]][nrow(wavecor[[4]]):1, ], main = "ROI Wavelet Correlation", axes = FALSE)
box()
image(p01[nrow(p01):1,], main = "Significance Map", col = c(0, 1), axes = FALSE)
box()
par(op)

library("eigenmodel")
diag(p01) <- NA
colnames(p01) <- rownames(p01) <- colnames(brain)
brainNet <- eigenmodel_mcmc(p01, R = 2, S = 1000, burn = 200, seed = 123)
evecs <- eigen(brainNet$ULU_postmean)$vec[, 1:2]

plot(evecs, type = "n", xlab = "Dimension 1", ylab = "Dimension 2", main = "Brain Network")
addlines(evecs, p01 > 0, col = "gray")  
text(evecs, labels = colnames(p01), cex = 0.7)
