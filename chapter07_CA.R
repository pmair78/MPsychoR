## ----------------------- Chapter 7: Correspondence Analysis 

## ----- simple correspondence analysis
superfan <- as.table(matrix(c(9, 12, 8, 1, 13, 1, 6, 20, 15, 4, 23, 18), ncol = 3))
attr(superfan, "dimnames") <- list(Band = c("Slayer", "Iron Maiden", "Metallica", "Judas Priest"), 
                                   Fan = c("Horst", "Helga", "Klaus"))
superfan

fit_chisq <- chisq.test(superfan)
fit_chisq
S <- fit_chisq$residuals
round(S, 3)

library("vcd")
mosaic(superfan, shade = TRUE)
P <- prop.table(superfan)
round(P, 3)         ## table with relative frequencies

r_mass <- margin.table(P, 1)
round(r_mass, 3)             ## row masses
c_mass <- margin.table(P, 2)
round(c_mass, 3)             ## column masses

r_profile <- prop.table(P, 1)
round(r_profile, 3) ## conditional relative frequencies rows
c_profile <- prop.table(P, 2)
round(c_profile, 3) ## conditional relative frequencies columns

ar_profile <- t(r_profile) %*% r_mass ## average row profile
round(as.vector(ar_profile), 3)
round(as.vector(c_mass), 3)           ## column masses
ac_profile <- c_profile %*% c_mass    ## average column profile
round(as.vector(ac_profile), 3)
round(as.vector(r_mass), 3)           ## row masses 

library("plot3D")
tc <- r_profile
scatter3D(x = tc[,1], y = tc[,2], z = tc[,3], xlab = "Horst", ylab = "Helga", zlab = "Klaus", colkey = FALSE, 
          col = 1, pch = 20, xlim = c(0,1), ylim = c(0,1), zlim = c(0,1), ticktype = "simple", type = "h", 
          phi = 40, theta = 50, main = "Row Profiles", bty = "g")
points3D(x = c(0,0,1), y = c(0,1,0), z = c(1,0,0), col = "red", add = TRUE, pch = 20)
lines3D(x = c(0,0,1,0), y = c(0,1,0,0), z = c(1,0,0,1), col = "red", add = TRUE)
text3D(x = tc[,1], y = tc[,2], z = tc[,3], labels = rownames(tc), pos = 3, add = TRUE, cex = 0.8, adj = -0.1) 

library("ggtern")
tf <- as.data.frame.matrix(superfan)
c_mass <- as.vector(c_mass)
lines <- data.frame(x = c(c_mass[1], 1-c_mass[3], 0),
                    y = c(1-c_mass[1], 0, c_mass[2]),
                    z = c(0, c_mass[3], 1-c_mass[2]),
                    xend = c(c_mass[1], c_mass[1], c_mass[1]),
                    yend = c(c_mass[2], c_mass[2], c_mass[2]),
                    zend = c(c_mass[3], c_mass[3], c_mass[3]), row.names = NULL)
gt <- ggtern(data = tf, aes(Horst, Helga, Klaus))
gt + geom_point() + theme_rgbw() + geom_text(label = rownames(tf), vjust = -0.5) +
geom_point(aes(x = c_mass[1], y = c_mass[2], z = c_mass[3]), colour = "red", size = 4) +
geom_segment(data = lines, aes(x = x, y = y, z = z, xend = xend, yend = yend, zend = zend), 
             color = 'red', size = 0.5) +
labs(title = "Ternary Plot")

## ------------------------------------------------------------------------
sqrt(sum((r_profile["Slayer",] - r_profile["Iron Maiden",])^2/ar_profile)) 
sqrt(sum((r_profile["Slayer",] - r_profile["Judas Priest",])^2/ar_profile))  
sqrt(sum((r_profile["Slayer",] - ar_profile)^2/ar_profile))  
sqrt(sum((r_profile["Iron Maiden",] - ar_profile)^2/ar_profile))  

library("anacor")
ca_fans <- anacor(superfan, ellipse = TRUE)
ca_fans
plot(ca_fans, main = "Symmetric CA Map")


library("MPsychoR")
data("HarvardPsych")
dim(HarvardPsych)  ## researchers in rows, words in columns
fit_HP <- anacor(HarvardPsych)

plot(fit_HP, main = "Harvard Psychology Faculty", asp = NULL, xlim = c(-4, 1))
plot(fit_HP, main = "Harvard Psychology Faculty (Zoom)", asp = NULL, xlim = c(0, 1), ylim = c(-1, 1))


## ------ multiple correspondence analysis
library("MPsychoR")
data("YouthDep")
cdisub <- YouthDep[, c("CDI15r", "CDI21r", "race")]

B <- burtTable(cdisub)
dim(B)   

library("ca")
fit_mca <- mjca(cdisub, lambda = "Burt")
plot(fit_mca, xlim = c(-0.5, 0.5))

## ----- configural frequency analysis
library("cfa")
data("HarvardPsych")
configs <- expand.grid(dimnames(HarvardPsych))
counts <- as.vector(HarvardPsych)
fit.cfa <- cfa(configs, counts, binom.test = TRUE, sorton = "n")
types <- fit.cfa$table[fit.cfa$table$sig.bin == TRUE, 1:3]
head(types, 10)

countdf <- as.data.frame(table(cdisub))
fit.cdi <- cfa(countdf[,1:3], countdf[,4])
fit.cdi$table[fit.cdi$table$sig.chisq == TRUE, 1:3]  #chi2-test
fit.cdi$table[fit.cdi$table$sig.z == TRUE, 1:3] ## z-test