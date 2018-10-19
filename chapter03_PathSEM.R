## ------------------- Chapter 3: Path Analysis and Structural Equation Models

## ----- multivariate regression
library("MPsychoR")
data("Bergh")
fitmvreg <- lm(cbind(EP, DP) ~ A1 + A2 + O1 + O2, data = Bergh)

library("car")
Manova(fitmvreg)

library("lavaan")
mvreg.model <- '
  EP ~ b11*A1 + b12*A2 + b13*O1 + b14*O2
  DP ~ b21*A1 + b22*A2 + b23*O1 + b24*O2'
fitmvreg2 <- sem(mvreg.model, data = Bergh)

library("semPlot")
semPaths(fitmvreg2, what = "est", edge.label.cex = 1, layout = "tree", residuals = FALSE, edge.color = 1, 
         esize = 1, rotation = 3, sizeMan = 8, asize = 2.5, fade = FALSE, optimizeLatRes = TRUE)
parameterEstimates(fitmvreg2)[c(1:8, 11),]

## ----- moderator and mediator models
library("MPsychoR")
data("Paskvan")
wintense.c <- scale(Paskvan$wintense, scale = FALSE)  ## center
fit.YX <- lm(cogapp ~ wintense.c, data = Paskvan)  ## Y on X
round(summary(fit.YX)$coefficients, 4)
pclimate.c <- scale(Paskvan$pclimate, scale = FALSE)  ## center
fit.YZ <- lm(cogapp ~ pclimate.c, data = Paskvan)  ## Y on Z
round(summary(fit.YZ)$coefficients, 4)

library("QuantPsyc")
fit.mod <- moderate.lm(x = wintense, z = pclimate, y = cogapp, data = Paskvan)
round(summary(fit.mod)$coefficients, 4)

fit.ss <- sim.slopes(fit.mod, Paskvan$pclimate)
round(fit.ss, 4)

library("mediation")
fit.MX <- lm(cogapp ~ wintense, data = Paskvan) 
fit.YXM <- lm(emotion ~ wintense + cogapp, data = Paskvan) 

set.seed(123)
fitmed <- mediation::mediate(fit.MX, fit.YXM, treat = "wintense", mediator = "cogapp", 
                             sims = 999, boot = TRUE, boot.ci.type = "bca")
summary(fitmed)

library("lavaan")
med.model <- '
  emotion ~ c*wintense + b*cogapp
  cogapp ~ a*wintense
  ind := a*b
  tot := ind+c
  prop := ind/tot'
set.seed(123)
fitmedsem <- lavaan::sem(med.model, Paskvan, se = "bootstrap", bootstrap = 999)
parameterEstimates(fitmedsem, zstat = FALSE, pvalue = FALSE, boot.ci.type = "bca.simple")[c(7,1,8,9),]

quantile(Paskvan$pclimate)
medmod.model <- '
## set of regressions
cogapp ~ a1*wintense + a2*pclimate + a3*wintense:pclimate  
emotion ~ c*wintense + b*cogapp

## conditional indirect effects 
cie.q1 := (a1 + a3*2)*b      ## first quartile
cie.q2 := (a1 + a3*3)*b      ## median
cie.q3 := (a1 + a3*3.5)*b    ## third quartile
'
set.seed(123)
fitmedmod <- lavaan::sem(medmod.model, data = Paskvan, se = "bootstrap", bootstrap = 999)

semPaths(fitmedmod, layout = "spring", asize = 2.5, sizeMan = 10, residuals = FALSE, nCharNodes = 7, 
         edge.label.cex = 1)

parameterEstimates(fitmedmod, zstat = FALSE, pvalue = FALSE, boot.ci.type = "bca.simple")[c(3, 4, 14:16),]



## ----- structural equation models
library("MPsychoR")
library("lavaan")
data("Bergh")
Bergh.model <- 'GP =~ EP + HP + DP + SP
                Agree =~ A1 + A2 + A3
                Open =~ O1 + O2 + O3
                GP ~ Agree + Open'
fitGP <- sem(Bergh.model, data = Bergh, estimator = "MLR")

semPaths(fitGP, what = "std", edge.label.cex = 0.7, esize = 1, intercepts = FALSE, rotation = 4, 
         edge.color = 1, asize = 2.5, sizeMan = 5, mar = c(1, 1.5, 1.5, 3), fade = FALSE)

summary(fitGP, standardized = TRUE, fit.measures = TRUE)

fit.free <- sem(Bergh.model, group = "gender", group.equal = c("intercepts"),
                group.partial = c("DP~1", "HP~1", "SP~1"), data = Bergh, estimator = "MLR")

fit.load <- sem(Bergh.model, group = "gender", group.equal = c("loadings", "intercepts"),
                group.partial = c("GP=~SP", "DP~1", "HP~1", "SP~1"), data = Bergh, estimator = "MLR")

fit.prestrict <- sem(Bergh.model, group = "gender", group.equal = c("intercepts", "regressions"),
                     group.partial = c("DP~1", "HP~1", "SP~1"), data = Bergh, estimator = "MLR")

anova(fit.free, fit.prestrict)

library("nonnest2")
fit.load1 <- update(fit.load, estimator = "ML")
fit.prestrict1 <- update(fit.prestrict, estimator = "ML")
compIC <- icci(fit.load1, fit.prestrict1)
compIC

vuongtest(fit.load1, fit.prestrict1)

## ----- latent growth models
library("lavaan")
library("aspect")
data("duncan")
model_shape <- ' 
   inter =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
   shape =~ 0*CIG_T1 + 1*CIG_T2 + CIG_T3 + CIG_T4'
fitCig1 <- growth(model_shape, data = duncan, estimator = "WLS")

semPaths(fitCig1, what = "std", edge.label.cex = 0.7, esize = 1, edge.color = 1, sizeMan = 6, asize = 2.5, 
         intercepts = FALSE, rotation = 4, mar = c(3, 5, 3.5, 5), fade = FALSE)

summary(fitCig1, header = FALSE)

model_lin <- '
  inter =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
  linear =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4'
fitCig2 <- growth(model_lin, data = duncan, estimator = "WLS")
parameterEstimates(fitCig2)[21, ]
round(fitMeasures(fitCig2)[c("rmsea", "cfi", "srmr")], 3)

model_quad <- '
  inter =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
  linear =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4
  quad =~ 0*CIG_T1 + 1*CIG_T2 + 4*CIG_T3 + 9*CIG_T4'
fitCig3 <- growth(model_quad, data = duncan, estimator = "WLS")
parameterEstimates(fitCig3)[28:29,]
round(fitMeasures(fitCig3)[c("rmsea", "cfi", "srmr")], 3)

anova(fitCig1, fitCig2)

model_ord <- '
  inter =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
  linear =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4
  CIG_T1 | 0*t1 + t2 + t3 + t4
  CIG_T2 | 0*t1 + t2 + t3 + t4'
fitCigord <- growth(model_ord, data = duncan, ordered = names(duncan)[5:8])

model_pc <- '
  cint =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
  clin =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4
  pint =~ 1*POT_T1 + 1*POT_T2 + 1*POT_T3 + 1*POT_T4
  plin =~ 0*POT_T1 + 1*POT_T2 + 2*POT_T3 + 3*POT_T4

  ## correlated errors
  CIG_T1 ~~ CIG_T2; CIG_T2 ~~ CIG_T3; CIG_T3 ~~ CIG_T4
  POT_T1 ~~ POT_T2; POT_T2 ~~ POT_T3; POT_T3 ~~ POT_T4

  ## fix error variances
  CIG_T1 ~~ rc*CIG_T1
  CIG_T2 ~~ rc*CIG_T2
  CIG_T3 ~~ rc*CIG_T3
  CIG_T4 ~~ rc*CIG_T4
  POT_T1 ~~ rp*POT_T1
  POT_T2 ~~ rp*POT_T2
  POT_T3 ~~ rp*POT_T3
  POT_T4 ~~ rp*POT_T4'
fitPC1 <- growth(model_pc, data = duncan, estimator = "WLS")
round(fitMeasures(fitPC1)[c("rmsea", "cfi", "srmr")], 3)

semPaths(fitPC1, what = "est", edge.label.cex = 0.7, esize = 1, edge.color = 1, sizeMan = 6, asize = 2.5, 
         intercepts = FALSE, rotation = 4, mar = c(3, 5, 3.5, 5), fade = FALSE)
inspect(fitPC1, "std")$psi
duncan$ALCavg <- rowMeans(duncan[, 9:12])

model_pca1 <- '
  cint =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
  clin =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4
  pint =~ 1*POT_T1 + 1*POT_T2 + 1*POT_T3 + 1*POT_T4
  plin =~ 0*POT_T1 + 1*POT_T2 + 2*POT_T3 + 3*POT_T4

  ## effects of alcohol on marijuana
  pint ~ ALCavg
  plin ~ ALCavg'
fitPCA1 <- growth(model_pca1, data = duncan, estimator = "MLR")
round(fitMeasures(fitPCA1)[c("rmsea", "cfi", "srmr")], 3)
parameterEstimates(fitPCA1)[17:18,]

semPaths(fitPCA1, what = "std", edge.label.cex = 0.7, edge.color = 1, esize = 1, sizeMan = 6, asize = 2.5,
         intercepts = FALSE, rotation = 4, thresholdColor = "red", mar = c(3, 5, 3.5, 5), fade = FALSE)

model_pca2 <- '
  cint =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
  clin =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4
  pint =~ 1*POT_T1 + 1*POT_T2 + 1*POT_T3 + 1*POT_T4
  plin =~ 0*POT_T1 + 1*POT_T2 + 2*POT_T3 + 3*POT_T4

  ## effects of alcohol on marijuana
  POT_T1 ~ ALC_T1
  POT_T2 ~ ALC_T2
  POT_T3 ~ ALC_T3
  POT_T4 ~ ALC_T4'
fitPCA2 <- growth(model_pca2, data = duncan, estimator = "MLR")
round(fitMeasures(fitPCA2)[c("rmsea", "cfi", "srmr")], 3)

semPaths(fitPCA2, what = "std", edge.label.cex = 0.7, edge.color = 1, esize = 1, sizeMan = 6, asize = 2.5,
         intercepts = FALSE, rotation = 4, thresholdColor = "red", mar = c(3, 5, 3.5, 5), fade = FALSE)
parameterEstimates(fitPCA2)[17:20,]