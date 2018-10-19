## ----------------------------------- Chapter 1: Classical Test Theory --------------------------
library("MPsychoR")
library("psych")
data("Rmotivation")

## ----- reliability
ind <- grep("hyb", colnames(Rmotivation))  
HybMotivation <- na.omit(Rmotivation[,ind]) ## item selection     
k <- ncol(HybMotivation)                    ## number of items
vcmat <- cov(HybMotivation)
sigma2_Xi <- tr(vcmat)
sigma2_X <- sum(vcmat)

cronalpha <- k/(k-1)*(1-sigma2_Xi/sigma2_X)
round(cronalpha, 2)
sqrt(sigma2_X)*sqrt(1-cronalpha)

alpha.hyb <- psych::alpha(HybMotivation)
round(alpha.hyb$total[1], 2)          ## Cronbach's alpha

res.glb <- glb(HybMotivation)
res.omega <- omega(HybMotivation, plot = FALSE)

## ----- generalizability theory
library("reshape2")
Hyb1 <- data.frame(HybMotivation, person = 1:nrow(HybMotivation))
Hyblong <- melt(Hyb1, id.vars = c("person"), variable.name = "item")
Hyblong$person <- as.factor(Hyblong$person)
summary(aov(value ~ person + item, data = Hyblong))
round((0.85-0.15)/0.85, 2)

icchyb <- ICC(HybMotivation)
sqrt((0.85-0.15)/19)
sqrt((31.88-0.15)/777)

library("lme4")
VarCorr(lmer(value ~ (1|person) + (1|item), data = Hyblong))

library("gtheory")
gfit <- gstudy(data = Hyblong, formula = value ~ (1|person) + (1|item))
dfit <- dstudy(gfit, colname.objects = "person", colname.scores = "value", data = Hyblong)
round(dfit$generalizability, 3)

data("Lakes")
phydat <- subset(Lakes, subtest == "physical")
phydat$item <- droplevels(phydat$item)
head(phydat)

formula <- score ~ (1|personID) + (1|raterID) + (1|item) + (1|personID:raterID) + (1|personID:item) + (1|raterID:item)
gfit <- gstudy(formula = formula, data = phydat)
gfit

dfit <- dstudy(gfit, colname.objects = "personID", colname.scores = "score", data = phydat)
dfit$components
dfit$var.error.abs
dfit$sem.abs
dfit$var.error.rel
dfit$sem.rel
dfit$dependability
dfit$generalizability