## -------------------- Chapter 13: Modeling Trajectories and Time Series

## ----- hidden Markov models
library("markovchain")
mcStats <- new("markovchain", 
  state = c("bored", "horrified", "sleeping", "somewhat awake"), 
  transitionMatrix = matrix(c(0.60, 0.05, 0.25, 0.10,
                              0.15, 0.30, 0.10, 0.45, 
                              0.10, 0.10, 0.50, 0.30,  
                              0.30, 0.10, 0.50, 0.10), 
                              byrow = TRUE, nrow = 4), 
  name = "StatsClass")
print(mcStats)
A <- mcStats@transitionMatrix
rowSums(A)

library("qgraph")
qgraph(A, edge.labels = TRUE, edge.color = "black")

pi <- c(0, 1, 0, 0)
pi * mcStats^2  ## 10 minutes later
pi * mcStats^3  ## 15 minutes later

p <- steadyStates(mcStats)       ## limiting distribution
p
p %*% A


library("MPsychoR")
library("depmixS4")
data("iatfaces")
p1dat <- subset(iatfaces, id == 1)

plot(1:80, log(p1dat$latency), type = "l", ylab = "Latency (log)", xlab = "Trial", main = "IAT Latency (Person 1)")
abline(v = 40, col = "gray", lty = 2)

set.seed(123)
p1obj1 <- depmix(log(latency) ~ 1, data = p1dat, nstates = 1)
p1fit1 <- fit(p1obj1)
p1obj2 <- depmix(log(latency) ~ 1, data = p1dat, nstates = 2)
p1fit2 <- fit(p1obj2)
p1obj3 <- depmix(log(latency) ~ 1, data = p1dat, nstates = 3)
p1fit3 <- fit(p1obj3)
p1obj4 <- depmix(log(latency) ~ 1, data = p1dat, nstates = 4)
p1fit4 <- fit(p1obj4)
c(BIC(p1fit1), BIC(p1fit2), BIC(p1fit3), BIC(p1fit4))

summary(p1fit2)
round(posterior(p1fit2)[35:45, ], 3)
table(state = p1fit2@posterior$state, block = p1dat$block)

set.seed(123)
p2dat <- subset(iatfaces, id == 2)
p2obj <- depmix(log(latency) ~ 1, data = p2dat, nstates = 2)
p2fit <- fit(p2obj)
p3dat <- subset(iatfaces, id == 3)
p3obj <- depmix(log(latency) ~ 1, data = p3dat, nstates = 2)
p3fit <- fit(p3obj)
p4dat <- subset(iatfaces, id == 4)
p4obj <- depmix(log(latency) ~ 1, data = p4dat, nstates = 2)
p4fit <- fit(p4obj)

table(state = p2fit@posterior$state, block = p2dat$block)
table(state = p3fit@posterior$state, block = p3dat$block)
table(state = p4fit@posterior$state, block = p4dat$block)

p1obj2 <- depmix(log(latency) ~ block, data = p1dat, ns = 1)
p1fit2a <- fit(p1obj2, verbose = FALSE)
summary(p1fit2a, which = "response")
lm(log(latency) ~ block, data = p1dat)

set.seed(123)
p1obj3 <- depmix(log(latency) ~ block, data = p1dat, ns = 2) 
p1fit2b <- fit(p1obj3, verbose = FALSE)
summary(p1fit2b)
c(BIC(p1fit2a), BIC(p1fit2b))

set.seed(123)
p1obj4 <- depmix(log(latency) ~ 1, data = p1dat, nstates = 2, transition = ~ block)
p1fit2c <- fit(p1obj4,emcontrol = em.control(maxit = 5000))
summary(p1fit2c, which = "transition")

eta1 <- matrix(getpars(p1fit2c)[3:6], 2, byrow = TRUE)
eta1          
eta2 <- matrix(getpars(p1fit2c)[7:10], 2, byrow = TRUE)
eta2          

exps1b1 <- exp(eta1[1,])
a11c <- 1/sum(exps1b1)
a11c
exps1b2 <- exp(eta1[2,])
a11ic <- 1/sum(exps1b2)
a11ic
exps2b1 <- exp(eta2[1,])
a22c <- 1-1/sum(exps2b1)
a22c
exps2b2 <- exp(eta2[2,])
a22ic <- 1-1/sum(exps2b2)
a22ic

Acong <- round(matrix(c(a11c, 1-a11c, 1-a22c, a22c), 2, byrow = TRUE), 5)
Aicong <- round(matrix(c(a11ic, 1-a11ic, 1-a22ic, a22ic), 2, byrow = TRUE), 5)
dimnames(Acong) <- dimnames(Aicong) <- list(c("fromS1", "fromS2"), c("toS1", "toS2"))
Acong           ## congruent condition
Aicong          ## incongruent condition


## ----- time series analysis
library("MPsychoR")
data("ageiat")
yts <- ts(ageiat, start = c(2007, 1), frequency = 12)   
plot(yts, ylab = "d-measure", main = "Age IAT Time Series")

library("forecast")
tr <- time(yts) - 2007
fitlm <- tslm(yts ~ tr)
fitlm
Acf(residuals(fitlm), main = "ACF of Residuals")

library("lmtest")
dwtest(ageiat ~ tr)

library("strucchange")
bp <- breakpoints(yts ~ tr)
bp
round(coef(bp), 3)

op <- par(mfrow = c(1,2))
plot(yts, col = "gray", main = "Structural Change Regression", ylab = "d-measure")
lines(window(fitted(bp), start = 2007, end = breakdates(bp)), col = "cadetblue", lwd = 1)
lines(window(fitted(bp), start = c(2010, 7), end = c(2015, 12)), col = "cadetblue", lwd = 1)
lines(confint(bp))
Acf(residuals(bp), main = "ACF Residuals")
par(op)

Box.test(residuals(bp), lag = 2)

ytslog <- log(yts)
op <- par(mfrow = c(2,1))
plot(yts, ylab = "d-measure", main = "IAT Time Series (Levels)")
plot(ytslog, ylab = "d-measure (log)", main = "IAT Time Series (Logs)")
par(op)

tsdec <- decompose(yts)
plot(tsdec)

## ------------------------------------------------------------------------
sadj <- yts - tsdec$seasonal

yts1 <- diff(yts, difference = 1)
yts2 <- diff(yts, difference = 2)
op <- par(mfrow = c(1,2))
plot(yts1, ylab = "difference d-measure", main = "First Difference Time Series")
plot(yts2, ylab = "difference d-measure", main = "Second Difference Time Series")
par(op)

library("tseries")
adf.test(yts)
adf.test(yts1)

op <- par(mfrow = c(1,2))
Acf(yts1, main = "Differenced IAT ACF")
Pacf(yts1, main = "Differenced IAT Partial ACF")
par(op)

iatima <- Arima(yts, order = c(0, 1, 2))
iatima

iatauto <- auto.arima(yts)
iatauto

BIC(iatima, iatauto)[,2]
AIC(iatima, iatauto)[,2]

tsdiag(iatima)
plot(iatauto$x, ylab = "D-measure", main = "IAT Time Series")
lines(fitted(iatauto), col = "salmon", lwd = 2)

cor(fitted(iatauto), yts)^2
plot(forecast(iatauto, h = 48), ylab = "d-measure", main = "Age IAT Forecasts")


aca <- c(rep(0, 39), rep(1, 69)) 
preaca <- window(yts, end = c(2010, 3)) ## pre-event series
tspre <- Arima(preaca, order = c(0,1,1), include.drift = TRUE)
preds <- forecast(tspre, h = 69)    ## post-event predictions

plot(yts, col = "gray", main = "Age IAT Time Series", ylab = "d-measure")
abline(v = time(yts)[40])
lines(preds$mean, col = "cadetblue", lwd = 1.5)
lines(preds$lower[,1], col = "cadetblue", lty = 2)
lines(preds$upper[,1], col = "cadetblue", lty = 2)
text(2010.3, 1.01, "Obamacare passed", pos = 4)

tsreg <- Arima(yts, order = c(0, 1, 1), include.drift = TRUE, xreg = aca)
print(coeftest(tsreg), 3, signif.legend = FALSE)

## ----- functional data analysis
library("MPsychoR")
library("fda.usc")
data("tension")
tension1 <- as.matrix(tension[,1:800])  ## tension time series
cond <- tension$cond                    ## condition
ftension <- fdata(tension1, argvals = seq(1, 80, length.out = 800), 
                  names = list(main = "Music tension", xlab = "Time (sec)", ylab = "Tension"))

ftensionNP <- optim.np(ftension)
plot(ftension, main = "Original Data")
plot(ftensionNP$fdata.est, main = "Smooth Data")

ftension1 <- ftensionNP$fdata.est
deriv1 <- fdata.deriv(ftension1)
op <- par(mfrow = c(2,1), mar = c(4,4,3,1) + 0.1)
plot(ftension1[7,], main = "Smooth Tension (Person 7)", ylab = "Functional Tension Values")
plot(deriv1[7,], main = "", ylab = "Velocity Tension Values")
par(op)

fsplit <- split(ftension1, cond)
Amean <- func.mean(fsplit$Auditory)
Vmean <- func.mean(fsplit$Visual)
AVmean <- func.mean(fsplit$AuditoryVisual)
Amedian <- func.med.FM(fsplit$Auditory)
Vmedian <- func.med.FM(fsplit$Visual)
AVmedian <- func.med.FM(fsplit$AuditoryVisual)

dev.new()
op <- par(mfrow = c(2,1))
plot(Amean, lwd = 2, main = "Mean Tension Trajectories", ylim = c(-2.5, 2.5))
lines(Vmean, col = "coral", lwd = 2)
lines(AVmean, col = "cadetblue", lwd = 2)
legend("bottomright", col = c(1, "coral", "cadetblue"), lty = 1, legend = c("Auditory", "Visual", "Auditory & Visual"))
plot(Amedian, lwd = 2, main = "Median Tension Trajectories", ylim = c(-2.5, 2.5))
lines(Vmedian, col = "coral", lwd = 2)
lines(AVmedian, col = "cadetblue", lwd = 2)
legend("bottomright", col = c(1, "coral", "cadetblue"), lty = 1, legend = c("Auditory", "Visual", "Auditory & Visual"))
par(op)

fsplit$AuditoryVisual$names$main <- NULL
set.seed(123)
AVboot <- fdata.bootstrap(fsplit$AuditoryVisual, draw = TRUE)
title("Bootstrap Mean (Auditory-Visual)")

set.seed(123)
tension1way <- fanova.onefactor(ftension1, cond, nboot = 50)
tension1way$pvalue

fdat <- as.data.frame(ftension1$data)
gdat <- as.data.frame(cond)
ctrAudio <- contr.treatment(3)
set.seed(222)
fitrpm <- fanova.RPm(fdat, ~ cond, gdat, RP = 30, contrast = list(cond = ctrAudio))
summary(fitrpm)

library("refund")
ftension2 <- fdata2fd(ftension1)
X <- model.matrix(~ cond)  ## auditory as baseline
tenreg <- fosr(fdobj = ftension2, X = X, lambda = 100)

refund:::plot.fosr(tenreg, titles = c("Intercept", "Auditory vs. Visual", "Auditory vs. Visual-Auditory"))

fpca <- fdata2pc(ftension1, ncomp = 2)
summary(fpca, biplot = FALSE)  
pcscores <- fpca$x[,1:2]
op <- par(mfrow = c(1,2))
plot(fpca$rotation[1,]$argvals, fpca$rotation$data[1,], type = "l", ylim = c(-0.3, 0.3), xlab = "Time", 
     ylab = "Loadings", main = "Functional PCA Loadings", lwd = 2)
lines(fpca$rotation$argvals, fpca$rotation$data[2,], type = "l", col = "gray", lwd = 2)
legend("topright", legend = c("PC1", "PC2"), lty = 1, col = c(1, "gray"))
plot(pcscores, type = "n", asp = 1, main = "Functional PC Scores")
cols <- c("black", "coral", "cadetblue")[as.numeric(cond)]
text(pcscores, col = cols)
legend("topright", legend = c("Auditory", "Visual", "Auditory & Visual"), text.col = c("black", "coral", "cadetblue"))
par(op)

fmean <- func.mean(ftension1)
pc1plus <- fmean$data[1,] + 3*fpca$rotation$data[1,]
pc1minus <- fmean$data[1,] - 3*fpca$rotation$data[1,]
pc2plus <- fmean$data[1,] + 3*fpca$rotation$data[2,]
pc2minus <- fmean$data[1,] - 3*fpca$rotation$data[2,]

op <- par(mfrow = c(1,2))
plot(fmean, lwd = 2, main = "Mean (PC1)", ylim = c(-2, 2))
lines(ftension1$argvals, pc1plus, col = "salmon")
lines(ftension1$argvals, pc1minus, col = "cadetblue")
legend("bottomright", legend = c("mean + PC1", "mean - PC1"), lty = 1, col = c("salmon", "cadetblue"))
plot(fmean, lwd = 2, main = "Mean (PC2)", ylim = c(-2, 2))
lines(ftension1$argvals, pc2plus, col = "salmon")
lines(ftension1$argvals, pc2minus, col = "cadetblue")
legend("bottomright", legend = c("mean + PC2", "mean - PC2"), lty = 1, col = c("salmon", "cadetblue"))
par(op)

