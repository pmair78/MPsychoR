## ------------------- Chapter 12: Parametric Cluster Analysis and Mixture Regression ---------------

## ----- model-based clustering approaches
library("MPsychoR")
library("mclust")
data("Rmotivation2")
Rwd <- Rmotivation2[, c(4:6)]   ## work design variables

set.seed(123)
clustbic <- mclustBIC(Rwd, G = 2:10)
clustbic
clusfit <- Mclust(Rwd, x = clustbic)
summary(clusfit, parameters = TRUE)

head(round(clusfit$z, 3))   ## soft
cluster <- as.factor(clusfit$classification)  ## hard

clustred <- MclustDR(clusfit)
plot(clustred, what = "boundaries", ngrid = 200)
plot(clustred, what = "density", dimens = 1)

library("lattice")
densityplot(~npkgs|cluster, data = Rmotivation2, col = "gray", index.cond = list(c(3, 4, 1, 2)))
histogram(~lists|cluster, data = Rmotivation2, col = "gray", index.cond = list(c(3, 4, 1, 2)))


library("MPsychoR")
library("poLCA")
data("AvalanchePrep")
formula <- cbind(info, discuss, gear, decision) ~ 1
set.seed(1)
fitlca1 <- poLCA(formula, data = AvalanchePrep, nclass = 1, nrep = 3)
fitlca2 <- poLCA(formula, data = AvalanchePrep, nclass = 2, nrep = 3)
fitlca3 <- poLCA(formula, data = AvalanchePrep, nclass = 3, nrep = 3, maxiter = 2000)
fitlca4 <- poLCA(formula, data = AvalanchePrep, nclass = 4, nrep = 3, maxiter = 5000)
c(fitlca1$bic, fitlca2$bic, fitlca3$bic, fitlca4$bic)
fitlca3

library("scatterplot3d")
poLCA:::plot.poLCA(fitlca3)

clusmemb <- as.factor(fitlca3$predclass)


library("flexmix")
data("zareki")
set.seed(123)
zarflex <- stepFlexmix(~ 1, data = zareki, k = 1:4, nrep = 3, 
                       model = list(FLXMRmultinom(addit7 ~ .), 
                                    FLXMRmultinom(addit8 ~ .),
                                    FLXMRmultinom(subtr3 ~ .),
                                    FLXMRmultinom(subtr7 ~ .),
                                    FLXMRglm(time ~ ., family = "gaussian")))

zarflex2 <- getModel(zarflex, "BIC")
cluster <- zarflex2@cluster
table(cluster)

catpars <- sapply(parameters(zarflex2)[1:4], exp)
colnames(catpars) <- c("addit7", "subtr3", "addit8", "subtr7")
catpars                     ## categorical variables
parameters(zarflex2)[[5]]   ## time variable

library("vcd")
A <- mosaic(addit7 ~ cluster, data = zareki, return_grob = TRUE)
B <- mosaic(addit8 ~ cluster, data = zareki, return_grob = TRUE)
C <- mosaic(subtr3 ~ cluster, data = zareki, return_grob = TRUE)
D <- mosaic(subtr7 ~ cluster, data = zareki, return_grob = TRUE)
mplot(A, B, C, D)

zarflexc <- flexmix(~ 1, data = zareki, cluster = posterior(zarflex2), concomitant = FLXPmultinom(~ class),
                    model = list(FLXMRmultinom(addit7 ~ .),
                                 FLXMRmultinom(subtr3 ~ .),
                                 FLXMRmultinom(addit8 ~ .),
                                 FLXMRmultinom(subtr7 ~ .),
                                 FLXMRglm(time ~ ., family = "gaussian")))
zarflexc@prior   ## weights with concomitant
zarflex2@prior   ## weights without concomitant
exp(zarflexc@concomitant@coef)

clusterc <- zarflexc@cluster
table(zareki$class, clusterc)  ## with concomitant variable
table(zareki$class, cluster)   ## without concomitant variable


## ----- mixture regression models
set.seed(123)
X <- rnorm(100)
Y1 <- 2*X + rnorm(100, sd = 0.5)
Y2 <- (-2)*X + rnorm(100, sd = 0.5)
toydat <- data.frame(X = c(X, X), Y = c(Y1, Y2))
lm(Y ~ X, data = toydat)

library("flexmix")
toymix <- flexmix(Y ~ X, k = 2, data = toydat)
parameters(toymix)

op <- par(mfrow = c(1,2))
with(toydat, plot(X, Y, main = "Scatterplot", cex = 0.8))
abline(lm(Y ~ X, data = toydat))
cols <- c("cadetblue", "salmon")
cmemb <- toymix@cluster
with(toydat, plot(X, Y, col = cols[cmemb], cex = 0.8, main = "Mixture Regression Fit"))
ind <- order(toydat$X)
pred1 <- predict(toymix)$Comp.1[ind]
lines(toydat$X[ind], pred1, col = "cadetblue", lwd = 2)
pred2 <- predict(toymix)$Comp.2[ind]
lines(toydat$X[ind], pred2, col = "coral", lwd = 2)
par(op)

library("MPsychoR")
data("KoreanSpeech")
set.seed(123)
koreamix <- flexmix(frequency ~ attitude|subject, k = 2, data = na.omit(KoreanSpeech))
table(koreamix@cluster)
parameters(koreamix) 

table(Cluster = clusters(koreamix), Gender = na.omit(KoreanSpeech)$gender)


library("gamair")
library("lattice")
data("brain")
brain <- brain[brain$medFPQ > 5e-3,] ## exclude outliers
trellis.par.set(regions = list(col = colorRampPalette(c('cadetblue4', 'white', 'coral4'))))
levelplot(log(medFPQ) ~ Y*X, data = brain)     

set.seed(123)
fitgammix <- flexmix(log(medFPQ) ~ s(Y, X, k = 30), model = FLXMRmgcv(), k = 3, 
                     data = brain, control = list(tolerance = 10^-3))
table(clusters(fitgammix))

xvals <- seq(40, 90, by = 0.5)    ## X
yvals <- seq(10, 60, by = 0.5)    ## Y
predvals <- predict(fitgammix, newdata = data.frame(X = xvals, Y = yvals))
color_transparent <- adjustcolor(clusters(fitgammix), alpha.f = 0.3) 
op <- par(mfrow = c(1,2))
with(brain, plot(log(medFPQ) ~ X, col = color_transparent , main = "Fit for X Coordinates", ylab = "Median FPQ [log]"))  ## Plot for X
lines(xvals, predvals[[1]], col = 1, lwd = 2)
lines(xvals, predvals[[2]], col = 2, lwd = 2)
lines(xvals, predvals[[3]], col = 3, lwd = 2)
legend("bottomright", col = 1:3, legend = paste("Cluster", 1:3), lty = 1)
with(brain, plot(log(medFPQ) ~ Y, col = color_transparent , main = "Fit for Y Coordinates", ylab = "Median FPQ [log]"))  ## Plot for Y
lines(yvals, predvals[[1]], col = 1, lwd = 2)
lines(yvals, predvals[[2]], col = 2, lwd = 2)
lines(yvals, predvals[[3]], col = 3, lwd = 2)
legend("bottomright", col = 1:3, legend = paste("Cluster", 1:3), lty = 1)
par(op)

## ----- dirichlet-based clustering 
library("profdpm")
koreamix2 <- profLinear(frequency ~ attitude, group = subject, data = na.omit(KoreanSpeech))
koreamix2$m           ## intercept and slope parameters
table(Cluster = koreamix2$clust, Gender = na.omit(KoreanSpeech)$gender)


library("MPsychoR")
library("randomNames")
gopraw <- readLines(paste0(path.package("MPsychoR"), "/GOPstatements.txt"))
set.seed(123)
rnames <- randomNames(length(gopraw), which.names = "first", ethnicity = 5, sample.with.replacement = FALSE)
names(gopraw) <- rnames

library("tm")
library("tidyr")
myStopwords <- c("beleive", "shld", "-", "wenot", "etc", "im", "conservatismthe", "fatherthe", 
                 "conservativebelieve", "governmentprolife2nd", "amendmentand", "valuessmall", "ive", 
                 "familyrepublican", "—government", "1st", "dont", "get", "given", "people", "better", 
                 "system", "always", "enough", "yet", "hand")
gopcorp <- Corpus(VectorSource(gopraw)) %>% 
  tm_map(content_transformer(tolower)) %>% 
  tm_map(removePunctuation) %>%
  tm_map(removeWords, c(stopwords("english"), myStopwords)) %>%
  tm_map(removeNumbers) %>%
  tm_map(stripWhitespace)

library("wordcloud")
set.seed(1)
op <- par(mar = c(5,4,4,2), oma = c(15,0,0,4))
wordcloud(gopcorp, colors = brewer.pal(8, "Dark2"), min.freq = 3, random.order = FALSE)
par(op)

DTmat <- DocumentTermMatrix(gopcorp)
dim(DTmat)

library("slam")
tfidf <- tapply(DTmat$v/row_sums(DTmat)[DTmat$i], DTmat$j, mean) * log2(nDocs(DTmat)/col_sums(DTmat > 0))
DTmat2 <- DTmat[, tfidf >= median(tfidf)]

ind <- which(rowSums(as.matrix(DTmat2)) > 0)
DTmat2 <- DTmat2[ind, ]
dim(DTmat2)

library("ldatuning")
Ktopics <- FindTopicsNumber(DTmat2, topics = 2:30, 
                            metrics = c("Arun2010", "CaoJuan2009", "Griffiths2004", "Deveaud2014"))
FindTopicsNumber_plot(Ktopics)

library("topicmodels")
K <- 6
goplda <- LDA(DTmat2, k = K, method = "Gibbs", control = list(seed = 123, iter = 50000, burnin = 1000))

postprob <- posterior(goplda)
pterms <- as.data.frame(t(postprob$terms))
round(head(pterms, 10), 4)

library("plyr")
library("tidyr")
w2 <- pterms %>% mutate(word = rownames(pterms)) %>% gather(topic, weight, -word)
n <- 50
pal <- rep(brewer.pal(9, "Greys"), each = ceiling(n / 9))[n:1]
op <- par(mfrow = c(3,2), mar = c(3,0,2,0))
for (i in 1:K) {
  w3 <- w2 %>% dplyr::filter(topic == i) %>% arrange(desc(weight))
  with(w3[1:n, ], wordcloud(word, freq = weight, scale = c(2, 0.5), random.order = FALSE, ordered.colors = TRUE, 
                            colors = pal))
  title(paste("GOP Topic", i))
}
par(op)

terms(goplda, 5)
topics(goplda)[1:6]