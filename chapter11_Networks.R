## -------------------- Chapter 11: Networks -------------------

## ----- network basics: relational data structures
library("igraph")
album_df <- matrix(c("Helga", "Gertrud", "Horst", "Klaus", "Horst", "Helga", "Gertrud", "Horst", "Klaus", "Horst", 
                     "Horst", "Helga", "Helga", "Klaus", "Klaus", "Klaus", "Gertrud", "Gertrud"), ncol = 2)
album_el <- graph.edgelist(album_df, directed = TRUE)
E(album_el)$weight = c(2, 3, 3, 4, 5, 1, 7, 7, 2)

album_ad <- get.adjacency(album_el, sparse = FALSE, attr = "weight")
album_ad

set.seed(123)
plot(album_el, vertex.size = 0, edge.arrow.size = 0.5, vertex.label.dist = 0.8, vertex.color = "black", 
     vertex.label.color = "black", vertex.label.cex = 1.5)  

## ----- correlation networks
library("MPsychoR")
data(Rogers)
cormat <- cor(Rogers)

library("qgraph")
cornet <- qgraph(cormat, layout = "spring", minimum = 0.2, graph = "cor", 
                 groups = list(Depression = 1:16, OCD = 17:26), color = c("white", "gray"), 
                 labels = colnames(Rogers), title = "Depression/OCD Correlation Network")

centralityPlot(cornet)

qgraph(cormat, layout = "spring", sampleSize = nrow(Rogers), graph = "glasso", 
       groups = list(Depression = 1:16, OCD = 17:26), color = c("white", "gray"), 
       labels = colnames(Rogers), title = "Depression/OCD Graphical Lasso")

## ----- latent network models
library("eigenmodel")
diag(cormat) <- NA   ## NA diagonals required
fitEM <- eigenmodel_mcmc(cormat, R = 2, S = 1000, burn = 200, seed = 123)  
evals <- colMeans(fitEM$L_postsamp)
evals

evecs <- eigen(fitEM$ULU_postmean)$vec[, 1:2]
cols <- c("coral", "cadetblue")
plot(evecs, type = "n", xlab = "Dimension 1", ylab = "Dimension 2", xlim = c(-0.30, 0), 
     main = "Depression/OCD Eigenmodel Network")
corthresh <- 0.2                     ## correlation threshold 
addlines(evecs, abs(cormat) > corthresh, col = "gray")
ind <- c(rep(1, 16), rep(2, 10))
text(evecs, labels = rownames(cormat), col = cols[ind], cex = 0.8)
legend("topright", legend = c("Depression", "OCD"), col = cols, pch = 19)

thresh <- 0.2
cormat01 <- ifelse(abs(cormat) > thresh, 1, 0)   
library("network")
cornet <- network(cormat01, matrix.type = "adjacency", directed = FALSE)

library("latentnet")
set.seed(111)
fitLN1 <- ergmm(cornet ~ euclidean(d = 2, G = 1))
summary(fitLN1)$bic$Z
fitLN2 <- ergmm(cornet ~ euclidean(d = 2, G = 2))
summary(fitLN2)$bic$Z
fitLN3 <- ergmm(cornet ~ euclidean(d = 2, G = 3))
summary(fitLN3)$bic$Z

plot(fitLN2, main = "Latent Class Network (2 Classes)", cluster.col = c("coral", "cadetblue"), labels = TRUE, 
     label.col = "darkgray", label.cex = 0.8, what = "pmean", xlab = "Dimension 1", ylab = "Dimension 2", 
     xlim = c(-10.5, 6.7))

clusmemb2 <- fitLN2$mkl$mbc$Z.pZK
dimnames(clusmemb2) <- list(colnames(cormat01), paste("Cluster", 1:2))
clusmemb2[c("comptime", "suicide", "weightgain"), ]

plot(fitLN3, main = "Latent Class Network", cluster.col = c("coral", "cadetblue", "darkgoldenrod"), 
     pie = TRUE, labels = TRUE, label.cex = 0.8, label.col = "darkgray", vertex.cex = 1.5, what = "pmean", 
     xlim = c(-7, 10), 
     xlab = "Dimension 1", ylab = "Dimension 2")

plot(fitLN3, main = "Latent Class Network", cluster.col = c("coral", "cadetblue", "darkgoldenrod"), what = "density")

## ----- Bayesian networks
library("bnlearn")
Rogers2 <- as.data.frame(apply(Rogers, 2, as.numeric))
set.seed(123)
fitBN <- hc(Rogers2, restart = 10, perturb = 100)

estrength <- arc.strength(fitBN, Rogers2, "bic-g")
head(estrength[order(estrength[,3]), ], 5)

strength.plot(fitBN, estrength, main = "Bayesian Network Depression/OCD", shape = "ellipse")

set.seed(123)
bootnet <- boot.strength(Rogers2, R = 500, algorithm = "hc")
head(bootnet)

avgnet <- averaged.network(bootnet, threshold = 0.85)
estrength <- arc.strength(avgnet, Rogers2, "bic-g")
strength.plot(avgnet, estrength, shape = "ellipse")

subedge <- head(bootnet[bootnet$strength > 0.95, ])
subedge

boottab <- bootnet[bootnet$strength > 0.85 & bootnet$direction > 0.5, ]  
astr <- boottab   
astr$strength <- astr$direction  
strength.plot(avgnet, astr, shape = "ellipse")
