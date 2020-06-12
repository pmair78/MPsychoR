## ------------------- Chapter 4: Item Response Theory
library("MPsychoR")
library("mirt")
data("zareki")
zarsub <- zareki[, grep("subtr", colnames(zareki))]

## ----- dimensionality assessment
library("Gifi")
prinzar <- princals(zarsub)
plot(prinzar, main = "Zareki Loadings")

library("psych")
nfactors(zarsub, n = 4, cor = "tet")

fitifa1 <- mirt(zarsub, 1, verbose = FALSE)
fitifa2 <- mirt(zarsub, 2, verbose = FALSE, TOL = 0.001)
anova(fitifa1, fitifa2, verbose = FALSE)

## ----- unidimensional dichotomous IRT models
library("eRm")
fitrasch1 <- RM(zarsub)
fitrasch1
round(fitrasch1$betapar, 3)
round(sort(-fitrasch1$betapar), 3)

timecat <- factor(zareki$time <= median(zareki$time), labels = c("fast", "slow"))
fitLR <- LRtest(fitrasch1, timecat)
fitLR
Waldtest(fitrasch1, timecat)
plotGOF(fitLR, ctrline = list(col = "gray"), conf = list())

fitrasch2 <- RM(zarsub[, -5])
LRtest(fitrasch2, timecat)

set.seed(123)
T1 <- NPtest(as.matrix(zarsub[, -5]), n = 1000, method = "T1")
T1
T11 <- NPtest(as.matrix(zarsub[, -5]), n = 1000, method = "T11")
T11

round(sort(-fitrasch2$betapar), 2)

plotjointICC(fitrasch2, xlab = "Subtraction Trait", main = "ICCs Subtraction Items")

zarppar <- person.parameter(fitrasch2)

zareki$theta <- zarppar$theta.table[,1]
summary(aov(theta ~ class, data = zareki))

library("ltm")
data("RWDQ")
fit2pl1 <- ltm(RWDQ ~ z1)
head(coef(fit2pl1))

RWDQ1 <- RWDQ[,-1]
fit2pl2 <- ltm(RWDQ1 ~ z1)
head(coef(fit2pl2))

item.fit(fit2pl2)

library("colorspace")
cols <- rainbow_hcl(8, c = 80)
plot(fit2pl2, item = 1:5, legend = TRUE, col = cols, lwd = 2)
round(coef(fit2pl2)[1:5, 2], 3)

ppars <- ltm::factor.scores(fit2pl2, resp.patterns = RWDQ1)$score.dat[, "z1"]

data("Wilmer")
VPMT <- Wilmer[,3:27]
fit3pl <- tpm(VPMT)
round(head(coef(fit3pl)), 3)
cols <- rainbow_hcl(6, c = 80)
plot(fit3pl, item = 1:6, legend = TRUE, col = cols, lwd = 2)

## ----- unidimensional polytomous IRT models
data("CEAQ")
itceaq <- CEAQ[,1:16] - 1

fitrsm <- RSM(itceaq)
ppar <- person.parameter(fitrsm)
ifit0 <- eRm::itemfit(ppar)
ifit0

ind <- match("ceaq10", colnames(itceaq))
itceaq1 <- itceaq[,-ind]
fitrsm1 <- RSM(itceaq1)
ppar1 <- person.parameter(fitrsm1)
ifit1 <- eRm::itemfit(ppar1)

ind <- match("ceaq15", colnames(itceaq1))
itceaq2 <- itceaq1[, -ind]
fitrsm2 <- RSM(itceaq2)
ppar2 <- person.parameter(fitrsm2)
ifit2 <- eRm::itemfit(ppar2)

library("mice")
set.seed(222)
imp <- mice(CEAQ)
gradevec <- complete(imp)$grade
levels(gradevec) <- c("grade56","grade56","grade78","grade78")
LRtest(fitrsm2, gradevec)

thpar <- thresholds(fitrsm2)
thpar

th1 <- thpar$threshpar[1:2]
plotICC(fitrsm2, item.subset = 1, xlab = "Empathy", col = rainbow_hcl(3, 80), lwd = 2)
abline(v = th1, col = "gray", lty = 2)

plotPImap(fitrsm2, latdim = "Empathy", main = "Person-Item Map CEAQ")

data("ASTI")
PGitems <- ASTI[ ,c(11,14,15,17,18,23)]   ## extract PG items
fitpcm <- PCM(PGitems)
thresholds(fitpcm)

plotPImap(fitpcm, latdim = "Presence/Growth", main = "Person-Item Map ASTI")

data("ASTI")
STitems <- ASTI[ ,c(2,4,7,13,16,24,25)]       ## ST items
stpcm <- gpcm(STitems, constraint = "rasch")  ## PCM
stgpcm <- gpcm(STitems)                       ## GPCM
anova(stpcm, stgpcm)                          ## LR-test

library("colorspace")
cols <- rainbow_hcl(4, 80)
op <- par(mfrow = c(2,2))
plot(stpcm, items = 1, main = c("ICC PCM (Item 1)"), col = cols)
plot(stpcm, items = 2, main = c("ICC PCM (Item 2)"), col = cols)
plot(stgpcm, items = 1, main = c("ICC GPCM (Item 1)"), col = cols)
plot(stgpcm, items = 2, main = c("ICC GPCM (Item 2)"), col = cols)
par(op)

fitgrm <- grm(STitems)
ppargrm <- ltm::factor.scores(fitgrm)

op <- par(mfrow = c(1,2))
plot(fitgrm, type = "OCCu", items = 1, main = "OCC Item 1")
abline(h = 0.5, lwd = 0.8, col = "gray")
matlines(rbind(coef(fitgrm)[[1]][1:3], coef(fitgrm)[[1]][1:3]), t(matrix(rep(c(0, 0.5), 3), ncol = 2, byrow = TRUE)), 
         col = "gray", lty = 2)
text(coef(fitgrm)[[1]][1:3], 0.05, labels = c(expression(beta[11]), expression(beta[12]), expression(beta[13])), 
     pos = 1, col = "gray", lty = 2)
plot(fitgrm, items = 1, main = "ICC Item 1")
par(op)

library("Gifi")
data("WilPat")
wpit15 <- WilPat[,1:15]
wpiprin <- princals(wpit15, ordinal = FALSE)
elim <- c("Nationalism", "Patriotism", "ChurchAuthority", "Obedience")
ind <- match(elim, colnames(wpit15))
wpitnew <- wpit15[, -ind]
wpihom <- homals(wpitnew)

plot(wpiprin, main = "Loadings WP Conservatism")
plot(wpihom, main = "Categories WP Conservatism")

library("mirt")
nrmwp <- mirt(wpitnew, 1, itemtype = "nominal")
ip1 <- itemplot(nrmwp, 1, main = colnames(wpitnew)[1], auto.key = list(text = c("disapprove", "approve", "don't know"), cex = 0.7))
ip2 <- itemplot(nrmwp, 3, main = colnames(wpitnew)[3], auto.key = list(text = c("disapprove", "approve", "don't know"), cex = 0.7))
ip3 <- itemplot(nrmwp, 4, main = colnames(wpitnew)[4], auto.key = list(text = c("disapprove", "approve", "don't know"), cex = 0.7))
ip4 <- itemplot(nrmwp, 10, main = colnames(wpitnew)[10], auto.key = list(text = c("disapprove", "approve", "don't know"), cex = 0.7))
print(ip1, split = c(1, 1, 2, 2), more = TRUE)
print(ip2, split = c(2, 1, 2, 2), more = TRUE)
print(ip3, split = c(1, 2, 2, 2), more = TRUE)
print(ip4, split = c(2, 2, 2, 2), more = FALSE)

M2(nrmwp)

## ---- item and test information
plot(nrmwp, type = "infotrace", main = "Item Information")
plot(nrmwp, type = "info")

## ----- IRT sample size determination
library("SimDesign")
m <- 20
n <- c(50, 75, 100, 150, 200, 300)
design <- as.data.frame(n)
set.seed(222)
poppars <- rbind(alpha = round(rlnorm(m, 0, 0.25), 2), d = round(rnorm(m), 2))

irtGenerate <- function(condition, fixed_objects = FALSE) {
  n <- condition$n
  a <- fixed_objects['alpha', ]
  d <- fixed_objects['d', ]
  dat <- simdata(a, d, n, itemtype = '2PL')
  return(dat)
}

irtAnalyze <- function(condition, dat, fixed_objects = NULL) {
    mod <- mirt(dat, 1, itemtype = '2PL', verbose = FALSE)
    simpars <- coef(mod, simplify = TRUE, digits = Inf)$items
    irtpars <- c(a = simpars[,1], d = simpars[,2])
    return(irtpars)
}

irtSummarize <- function(condition, results, fixed_objects = NULL) {
  apop <- fixed_objects['alpha', ]
  dpop <- fixed_objects['d', ]
  simrmse <- RMSE(results, c(apop, dpop))
  out <- c(RMSE = simrmse)
  return(out)
}

set.seed(222)
simres <- runSimulation(design, replications = 100, parallel = TRUE,  generate = irtGenerate, analyse = irtAnalyze, 
                        summarise = irtSummarize, packages = c('mirt'), fixed_objects = poppars)
simres

colind <- grep(".a.", colnames(simres))
sima <- as.data.frame(simres[, colind])
nvec <- as.numeric(levels(simres$n))
matplot(nvec, log(sima), type = "l", col = 1, lty = 1, ylab = "log(RMSE)", xlab = "sample size", 
        main = "2-PL Monte Carlo", xaxt = "n")
axis(1, at = nvec)

meanRMSE <- rowMeans(sima)
names(meanRMSE) <- n
round(meanRMSE, 2)

## ----- differential item functioning
library("lordif")
library("MPsychoR")
data("YouthDep")
cdi <- YouthDep[,1:26]        ## extract CDI items
cdiDIF <- lordif(cdi, YouthDep$race, criterion = "Chisqr")
cdiDIF$stats[1:3, 1:5]

plot(cdiDIF, labels = c("White", "Black", "Asian", "Latino"))
head(cdiDIF$ipar.sparse, 10)
ppar <- cdiDIF$calib.sparse$theta

library("psychotree")
library("psychotools")
data("MathExam14W")
itmath <- as.list.data.frame(MathExam14W$solved)
covars <- MathExam14W[,3:9]
mex <- data.frame(solved = itemresp(itmath), covars)
mex <- subset(mex, nsolved > 0 & nsolved < 13)
mex$tests <- ordered(mex$tests)
mex$nsolved <- ordered(mex$nsolved)
mex$attempt <- ordered(mex$attempt)

set.seed(1)
mrt <- raschtree(solved ~ group + tests + nsolved + gender + attempt + study + semester, 
                 data = mex, vcov = "info", minsize = 50, ordinal = "l2", nrep = 1e5)
plot(mrt)
round(itempar(mrt)[,1:4], 2)


## ----- multidimensional IRT models
library("MPsychoR")
library("ltm")
data("RWDQ")
RWDQ1 <- RWDQ[,-1]  ## eliminate first item (misfit)
irtpar <- ltm(RWDQ1 ~ z1)
fapar <- ltm(RWDQ1 ~ z1, IRT.param = FALSE)
round(head(cbind(coef(irtpar), coef(fapar))), 3)

irtppar <- factor.scores(irtpar)$score.dat$z1
fappar <- factor.scores(fapar)$score.dat$z1
identical(irtppar, fappar)

library("MPsychoR")
library("Gifi")
data("zareki")
itzareki <- zareki[, 1:16]
przar <- princals(itzareki)
op <- par(mfrow = c(1,2))
plot(przar)
plot(przar, "screeplot")
par(op)

zar1d <- mirt(itzareki, 1, itemtype = "2PL")
zar2d <- mirt(itzareki, 2, itemtype = "2PL")

anova(zar1d, zar2d)
M2(zar2d)

ifit2D2pl <- mirt::itemfit(zar2d)
ifit2D2pl[ifit2D2pl[, 4] < 0.05, ]  ## misfitting items

summary(zar2d, rotate = "varimax")
summary(zar2d, rotate = "oblimin")

itemplot(zar2d, 3, main = "ICS addit3", rot = list(xaxis = -70, yaxis = 50, zaxis = 10))
head(MDIFF(zar2d))
head(fscores(zar2d))

class2 <- zareki$class
levels(class2) <- c("second", "thirdfourth", "thirdfourth")
modMG <- multipleGroup(itzareki, model = 2, group = class2, SE = TRUE, verbose = FALSE)
astiDIF <- DIF(modMG, c('a1', 'd'), Wald = TRUE, p.adjust = 'fdr')
round(astiDIF$adj_pvals[astiDIF$adj_pvals < 0.05], 4)

data("ASTI")
itasti <- ASTI[, 1:25]
modASTI <- mirt.model('
  si = 10,19,20,21
  pm = 1,5,9,22
  na = 3,6,8,12
  st = 2,4,7,13,16,24,25
  pg = 11,14,15,17,18,23
  COV = si*pm*na*st*pg
')
asti5d <- mirt(itasti, model = modASTI, itemtype = 'graded', method = 'MHRM', SE.type = 'MHRM', verbose = FALSE)

astisum <- summary(asti5d, verbose = FALSE)
round(astisum$fcor, 3)
round(astisum$rotF["ASTI18",], 4)
M2(asti5d, QMC = TRUE)


## ----- longitudinal IRT models
library("MPsychoR")
data("SDOwave")
SDO3 <- SDOwave[,c(1:12)]
SDO3 <- sapply(SDO3, function(bin) ifelse(bin == 1, 0, 1))

library("eRm")
sdolltm1 <- LLTM(SDO3, mpoints = 3)
sdolltm1$W
summary(sdolltm1)

W0 <- sdolltm1$W[,-c(4:5)]
sdolltm0 <- LLTM(SDO3, W0)
anova(sdolltm0, sdolltm1)

group <- rep(c(1, 2), each = nrow(SDO3)/2)  ## fake covariate
sdolltm2 <- LLTM(SDO3, mpoints = 3, group = group)
anova(sdolltm2, sdolltm1)

sdollra <- LLRA(SDO3, mpoints = 3)
summary(sdollra)

SDO2 <- SDOwave[,1:8]
SDO2 <- sapply(SDO2, function(co) cut(co, c(0,1,2,3,4,7), labels = 1:5))
class(SDO2) <- "numeric"

iloads <- rep(1:4, 2)
ttmodel <- mirt.model('
  T1996 = 1-4
  T1997 = 5-8
  COV = T1996*T1997, T1997*T1997
  MEAN = T1997
  CONSTRAIN = (1, 5, d1), (2, 6, d1), (3, 7, d1), (4, 8, d1),
              (1, 5, d2), (2, 6, d2), (3, 7, d2), (4, 8, d2),
              (1, 5, d3), (2, 6, d3), (3, 7, d3), (4, 8, d3),
              (1, 5, d4), (2, 6, d4), (3, 7, d4), (4, 8, d4)')
fitSDO2 <- bfactor(SDO2, iloads, ttmodel, SE = TRUE)
round(coef(fitSDO2)$GroupPars[,2], 4)

itloads <- rep(1:4, 2)
modgr <- mirt.model('
    Intercept = 1-8
    Slope = 1-8
    COV = Intercept*Slope, Intercept*Intercept, Slope*Slope
    MEAN = Intercept, Slope
    START = (1-8, a1, 1), (1-4, a2, 0), (5-8, a2, 1)
    FIXED = (1-8, a1), (1-4, a2), (5-8, a2)
    CONSTRAIN = (1, 5, d1), (2, 6, d1), (3, 7, d1), (4, 8, d1),
                (1, 5, d2), (2, 6, d2), (3, 7, d2), (4, 8, d2),
                (1, 5, d3), (2, 6, d3), (3, 7, d3), (4, 8, d3),
                (1, 5, d4), (2, 6, d4), (3, 7, d4), (4, 8, d4)')
fitGIRT <- bfactor(SDO2, itloads, modgr, SE = TRUE)
coef(fitGIRT)$GroupPars[, c(2,8)]

library("MPsychoR")
data("RWDQ")
RWDQ1 <- RWDQ[, 2:9]  ## select 8 items
freq2pl <- ltm(RWDQ1 ~ z1)
intstart <- -coef(freq2pl)[,1]
discstart <- coef(freq2pl)[,2]

## ----- Bayesian IRT
library("MCMCpack")
chainWDQ1 <- MCMCirt1d(RWDQ1, burnin = 5000, mcmc = 50000, seed = 111, AB0 = 0.15, store.item = TRUE, 
                       store.ability = FALSE, verbose = TRUE, alpha.start = intstart, beta.start = discstart)

chainWDQ2 <- MCMCirt1d(RWDQ1, burnin = 5000, mcmc = 50000, seed = 222, AB0 = 0.15, store.item = TRUE, 
                       store.ability = FALSE, verbose = TRUE, alpha.start = intstart, beta.start = discstart)

chainWDQ3 <- MCMCirt1d(RWDQ1, burnin = 5000, mcmc = 50000, seed = 333, AB0 = 0.15, store.item = TRUE,
                       store.ability = FALSE, verbose = TRUE, alpha.start = intstart, beta.start = discstart)
WDQlist <- mcmc.list(chainWDQ1, chainWDQ2, chainWDQ3)

plot(WDQlist, auto.layout = FALSE, ask = FALSE)

data("HRB")
HRB1 <- HRB[rowSums(HRB) > 0, ]
rownames(HRB1) <- 1:nrow(HRB1)
time <- rep(1:5, each = 4)
fitdyn <- MCMCdynamicIRT1d(HRB1, item.time.map = time, mcmc = 20000, burnin = 5000, seed = 111, 
                           store.ability = TRUE, store.item = FALSE, verbose = TRUE)
dynsum <- summary(fitdyn)

nt <- 5                            ## number of time points
postmean <- dynsum$statistics[,1]  ## posterior means
pertraj <- t(matrix(postmean[1:(nrow(HRB1)*nt)], nrow = nt))
colnames(pertraj) <- paste0("T", 1:5)
round(head(pertraj), 3)

matplot(t(pertraj), type = "l", lty = 1, cex = 0.8, col = adjustcolor(1, alpha.f = 0.3), 
        ylab = "person parameter", xlab = "time points", main = "Individual Trajectories")
