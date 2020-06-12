## -------------------- Chapter 9: multidimensional scaling -------------------

## ----- exploratory MDS
library("MPsychoR")
library("smacof")
data(Wenchuan)
Wdelta <- dist(t(Wenchuan))     ## Euclidean distances                   
fit.wenchuan1 <- mds(Wdelta, type = "ordinal")  ## MDS fit  
fit.wenchuan1
plot(fit.wenchuan1, main = "Wenchuan MDS")

set.seed(123)
rsvec <- randomstress(n = attr(Wdelta, "Size"), ndim = 2, nrep = 500, type = "ordinal")
mean(rsvec)
mean(rsvec) - 2*sd(rsvec)

set.seed(123)
permmds <- permtest(fit.wenchuan1, data = Wenchuan, method.dat = "euclidean", nrep = 500, verbose = FALSE)
permmds

n <- attr(Wdelta, "Size")
svec <- NULL
for (i in 1:(n-1)) {
  svec[i] <- mds(Wdelta, ndim = i, type = "ordinal")$stress
}

plot(1:(n-1), svec, type = "b", main = "MDS Scree Plot", pch = 20, xlab = "Number of Dimensions", ylab = "Stress")

set.seed(123)
fit.wenchuan <- NULL  ## 100 random starts
for(i in 1:100) {
  fit.wenchuan[[i]] <- mds(Wdelta, type = "ordinal", init = "random") 
}
## extract the best solution
ind <- which.min(sapply(fit.wenchuan, function(x) x$stress))
fit.wenchuan2 <- fit.wenchuan[[ind]]
fit.wenchuan2$stress        ## lowest stress (random start)
fit.wenchuan1$stress        ## stress (classical scaling start)

fit.wenchuan3 <- mds(Wdelta, type = "interval")
fit.wenchuan3

op <- par(mfrow = c(1,2))
plot(fit.wenchuan2, plot.type = "Shepard", main = "Shepard Diagram (Ordinal MDS)")
plot(fit.wenchuan3, plot.type = "Shepard", main = "Shepard Diagram (Interval MDS)")
par(op)

plot(fit.wenchuan2, plot.type = "stressplot", main = "Wenchuan Stress-per-Point")

library("colorspace")
set.seed(123)
bootWen <- bootmds(fit.wenchuan2, data = Wenchuan, method.dat = "euclidean", nrep = 100)
cols <- rainbow_hcl(17, l = 60)
plot(bootWen, col = cols)

colpal <- c(rainbow_hcl(3, c = 100))
pal <- palette(colpal)
memb <- c(1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3) 
plot(fit.wenchuan2, main = "Wenchuan Configuration", label.conf = list(label = FALSE), col = memb)
legend("topright", legend = c("intrusive recollection", "avoidance/numbing", "arousal"), col = colpal, pch = 20)
abline(-0.05, 0.2, col = "lightgray")
lines(c(0.1, -1),c(-0.03, 1), col = "lightgray")
palette(pal)
colpal <- c(rainbow_hcl(3, c = 100, l = 30))
pal <- palette(colpal)
text(fit.wenchuan2$conf[-7,], labels = rownames(fit.wenchuan2$conf)[-7], col = 
memb[-7], pos = 3, cex = 0.8)
text(fit.wenchuan2$conf[7,1], fit.wenchuan2$conf[7,2], labels = 
rownames(fit.wenchuan2$conf)[7], col = memb[7], pos = 1, cex = 0.8)
palette(pal)

## ----- confirmatory MDS
library("MPsychoR")
library("smacof")
ocpD <- sim2diss(OCP[,1:54])
fit.ocp1 <- mds(ocpD, type = "ordinal")

library("colorspace")
library("gplots")
cols <- as.numeric(OCP[,55])   ## colors
colpal <- c(rainbow_hcl(4, c = 100, l = 40), col2hex("lightgray"))
pal <- palette(colpal)
op <- par(mar = c(5.1, 4.1, 4.1, 10.1), xpd = TRUE)
plot(fit.ocp1, col = cols, label.conf = list(col = cols), main = "OCP Configuration (Unrestricted)")
legend("right", legend = levels(OCP[,55]), inset=c(-0.48,0), col = colpal, pch = 19, title = "Classes")
par(op)

Z <- OCP[,56:57]
head(Z)
fit.ocp2 <- smacofConstraint(ocpD, type = "ordinal", constraint = "diagonal", init = fit.ocp1$conf,
                             external = Z, constraint.type = "ordinal")

op <- par(mar=c(5.1, 4.1, 4.1, 10.1), xpd = TRUE)
plot(fit.ocp2, col = cols, label.conf = list(col = cols), main = "OCP Configurations (Regional Restrictions)")
legend("right", legend = levels(OCP[,55]), inset=c(-0.48,0), col = colpal, pch = 19, title = "Facets")
par(op)
palette(pal)

fit.ocp3 <- smacofSphere(ocpD, penalty = 1000, type = "ordinal")

colpal <- c(rainbow_hcl(4, c = 100, l = 40), col2hex("lightgray"))
pal <- palette(colpal)
op <- par(mar=c(5.1, 4.1, 4.1, 10.1), xpd = TRUE)
plot(fit.ocp3, col = cols, label.conf = list(col = cols))
legend("right", legend = levels(OCP[,55]), inset=c(-0.48,0), col = colpal, pch = 19, title = "Facets")
par(op)
palette(pal)

## ----- unfolding
library("MPsychoR")
library("smacof")
fitSchwartz <- unfolding(indvalues, type = "interval")

plot(fitSchwartz, label.conf.rows = list(label = FALSE))
circ <- fitCircle(fitSchwartz$conf.col[,1], fitSchwartz$conf.col[,2])
draw.circle(circ$cx, circ$cy, circ$radius, border = "gray", lty = 2)

## ---- MDS extensions and related models
library("MPsychoR")
data("Pashkam")
fitcolor <- mds(Pashkam$color, type = "interval") ## color task
fitshape <- mds(Pashkam$shape, type = "interval") ## shape task

fitproc <- Procrustes(X = fitcolor$conf, Y = fitshape$conf)
fitproc

op <- par(mfrow = c(1,2))
plot(fitcolor, col = "cadetblue", label.conf = list(col = "cadetblue"), main = "Separate MDS Configurations", 
     ylim = c(-1, 1))
points(fitshape$conf, col = "coral1", pch = 20)
text(fitshape$conf, labels = rownames(fitshape$conf), col = "coral1", cex = 0.8, pos = 3)
legend("bottomright", col = c("cadetblue", "coral"), legend = c("Color Task", "Shape Task"), pch = 19)
plot(fitproc, main = "Procrustes Configuration", ylim = c(-1, 1))
legend("bottomright", col = c("cadetblue", "coral"), legend = c("Color Task", "Shape Task"), pch = 19)
par(op)


data("NeuralActivity")
fitNeuro <- indscal(NeuralActivity[1:10], type = "interval", itmax = 5000)

library("plotfunctions")
data("NeuralScales")
cols <- cut(NeuralScales$Social, 5, labels = FALSE)
colpal <- rev(sequential_hcl(5))
pal <- palette(colpal)
plot(fitNeuro, col = cols, label.conf = list(col = cols), main = "Neural Activity Space")
gradientLegend(valRange = c(1,5), color = colpal, n.seg = 1:5, dec = 0, inside = TRUE, pos.num = 4)
palette(pal)

