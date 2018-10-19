## ------------------------- Chapter 5: Preference Modeling ----------------

## ----- models for paired comparisons
library("MPsychoR")
data("bandpref")
bandpref

library("BradleyTerry2")
bandsBT <- BTm(cbind(Win1, Win2), Band1, Band2, data = bandpref)
bandsAbil <- BTabilities(bandsBT)
round(sort(bandsAbil[,1]), 3)
alphas <- exp(bandsAbil[,1])
sort(alphas)
alphas["Rush"]/alphas["Death"]

library("psychotree")
data("Topmodel2007")
topmtree <- bttree(preference ~ age + gender + q1 + q2 + q3, data = Topmodel2007)
plot(topmtree)
pcmat <- as.matrix(Topmodel2007[[1]])
pcvec <- as.vector(pcmat)
pcvec[pcvec == -1] <- 0

library("stringr")
modnames <- t(str_split(colnames(pcmat), ":", simplify = TRUE))
preds <- Topmodel2007[2:6]
ind <- sapply(preds, is.factor)
preds[ind] <- sapply(preds[ind], function(f) c(0:(length(levels(f))-1))[f])
preds <- scale(preds)
rownames(preds) <- paste0("P", 1:nrow(preds))
head(preds)

library("BTLLasso")
sid <- rep(rownames(preds), ncol(pcmat))
mfirst <- rep(modnames[1,], each = nrow(preds))
msecond <- rep(modnames[2,], each = nrow(preds))
BTresp <- response.BTLLasso(response = pcvec, subject = sid, first.object = mfirst, second.object = msecond)

set.seed(123)
lambda <- exp(seq(log(10), log(1), length = 20))-1
modLasso <- cv.BTLLasso(Y = BTresp, X = preds, lambda = lambda, folds = 5, trace = FALSE)
modLasso
plot(modLasso, plots_per_page = 5, which = 2:6, ask_new = FALSE)

optind <- which.min(modLasso$criterion)
lambda2 <- lambda[(optind-3):(optind+3)]
set.seed(111)
bootLasso <- boot.BTLLasso(modLasso, lambda = lambda2, cores = 4, B = 50, trace = FALSE, trace.cv = FALSE)
plot(bootLasso, plots_per_page = 6, ask_new = FALSE)

## ----- log-linear models for preference
library("prefmod")
music5 <- music[,c("mood", "regg", "rap", "hvym", "conr", "sex")]
head(music5)

pattmus <- pattL.fit(music5, nitems = 5, formel = ~sex)
pattmus
round(patt.worth(pattmus, outmat = "lambda"), 3)

pworth <- patt.worth(pattmus)
colnames(pworth) <- c("male", "female")
plot(pworth, main = "Musical Preferences")

library("MPsychoR")
data("learnemo")
head(learnemo)
pattemo <- pattPC.fit(learnemo, nitems = 5, formel = ~ sex, 
                      obj.names = c("enjoyment", "pride", "anger", "anxiety", "boredom"))
pattemo

pwemo <- patt.worth(pattemo)
colnames(pwemo) <- c("male", "female")
plot(pwemo, main = "Achievement Emotions")


carconf1 <- carconf[, c(1:6, 8)]
head(carconf1)
pattcar <- pattR.fit(carconf1, nitems = 6, formel = ~ age)
pattcar
pwcar <- patt.worth(pattcar)
colnames(pwcar) <- c("17-29 years", "30-49 years", "50+ years")
plot(pwcar, main = "Car Ratings")

