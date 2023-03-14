## loading dataset and reading R sources
load("../data/guide.rda")
source("../code/india_GEE.R")
library(gee)

## fitted model
fm <- gee(bothered ~ gender + age + dayacc + severe + toilet, id = practID, data = guide,
  family = binomial("logit"), corstr = "exchangeable", scale.fix = TRUE, scale.value = 1.)
# next five lines are required for our influence measures...
x <- guide[,c(5,7,9,10,11)]
x <- cbind(1, x)
colnames(x) <- fm$xnames
x$gender1 <- as.numeric(guide$gender) - 1
x <- as.matrix(x)
nobs <- nrow(x)

## removing pacients 8, 44, 64, 86, 88 and 122
rm.08 <- rm.44 <- rm.64 <- rm.86 <- rm.88 <- rm.122 <- rep(TRUE, nobs)
rm.08[8] <- rm.44[44] <- rm.64[64] <- rm.86[86] <- rm.88[88] <- rm.122[122] <- FALSE
fm.08 <- gee(bothered ~ gender + age + dayacc + severe + toilet, id = practID, subset = rm.08,
  data = guide, family = binomial("logit"), corstr = "exchangeable", scale.fix = TRUE, scale.value = 1.)
fm.44 <- gee(bothered ~ gender + age + dayacc + severe + toilet, id = practID, subset = rm.44,
  data = guide, family = binomial("logit"), corstr = "exchangeable", scale.fix = TRUE, scale.value = 1.)
fm.64 <- gee(bothered ~ gender + age + dayacc + severe + toilet, id = practID, subset = rm.64,
    data = guide, family = binomial("logit"), corstr = "exchangeable", scale.fix = TRUE, scale.value = 1.)
fm.86 <- gee(bothered ~ gender + age + dayacc + severe + toilet, id = practID, subset = rm.86,
    data = guide, family = binomial("logit"), corstr = "exchangeable", scale.fix = TRUE, scale.value = 1.)
fm.88 <- gee(bothered ~ gender + age + dayacc + severe + toilet, id = practID, subset = rm.88,
    data = guide, family = binomial("logit"), corstr = "exchangeable", scale.fix = TRUE, scale.value = 1.)
fm.122 <- gee(bothered ~ gender + age + dayacc + severe + toilet, id = practID, subset = rm.122,
    data = guide, family = binomial("logit"), corstr = "exchangeable", scale.fix = TRUE, scale.value = 1.)

## Table 1 (slightly recrafted)
tab1 <- matrix(0, nrow = 14, ncol = 6)
tab1[1,] <- fm$coef
tab1[2,] <- sqrt(diag(fm$robust))
tab1[3,] <- fm.08$coef
tab1[4,] <- sqrt(diag(fm.08$robust))
tab1[5,] <- fm.44$coef
tab1[6,] <- sqrt(diag(fm.44$robust))
tab1[7,] <- fm.64$coef
tab1[8,] <- sqrt(diag(fm.64$robust))
tab1[9,] <- fm.86$coef
tab1[10,] <- sqrt(diag(fm.86$robust))
tab1[11,] <- fm.88$coef
tab1[12,] <- sqrt(diag(fm.88$robust))
tab1[13,] <- fm.122$coef
tab1[14,] <- sqrt(diag(fm.122$robust))
colnames(tab1) <- fm$xnames
rownames(tab1) <- c("all", "SE", "8", "SE", "44", "SE", "64", "SE", "86", "SE", "88", "SE", "122", "SE")
tab1
#     (Intercept)  gender     age dayacc severe toilet
# all     -3.0543 -0.7453 -0.6756 0.3918 0.8124 0.1078
# SE       0.9590  0.6003  0.5606 0.0932 0.3590 0.0989
# 8       -2.6561 -1.0772 -0.9126 0.4579 0.6461 0.1433
# SE       0.8229  0.5603  0.5771 0.0985 0.3303 0.1166
# 44      -3.3655 -0.7602 -0.7821 0.3912 0.7225 0.2090
# SE       1.0284  0.6373  0.5941 0.1013 0.3509 0.1004
# 64      -3.4054 -0.5006 -0.5461 0.4084 0.8235 0.1044
# SE       1.0178  0.5755  0.5882 0.0984 0.3571 0.0963
# 86      -3.0876 -0.7354 -0.6707 0.3949 0.8150 0.1083
# SE       0.9585  0.6022  0.5621 0.0992 0.3602 0.0992
# 88      -3.5041 -0.7666 -0.7423 0.4266 0.9560 0.1120
# SE       0.9920  0.6252  0.5667 0.0976 0.3675 0.1059
# 122     -3.1799 -0.7477 -0.6983 0.3755 0.9444 0.0965
# SE       1.0181  0.5810  0.5498 0.0908 0.3946 0.0974

alpha <- rep(0, 7)
alpha[1] <- fm$work[2,1]
alpha[2] <- fm.08$work[2,1]
alpha[3] <- fm.44$work[2,1]
alpha[4] <- fm.64$work[2,1]
alpha[5] <- fm.86$work[2,1]
alpha[6] <- fm.88$work[2,1]
alpha[7] <- fm.122$work[2,1]
names(alpha) <- c("all", "8", "44", "64", "86", "88", "122")
alpha
#    all      8     44     64     86     88    122
# 0.0932 0.0863 0.1020 0.1031 0.0785 0.0314 0.1035

vol1 <- rep(0, 7) # volume of 'robust' estimator of cov(coef)
vol1[1] <- det(fm$robust)
vol1[2] <- det(fm.08$robust)
vol1[3] <- det(fm.44$robust)
vol1[4] <- det(fm.64$robust)
vol1[5] <- det(fm.86$robust)
vol1[6] <- det(fm.88$robust)
vol1[7] <- det(fm.122$robust)
names(vol1) <- c("all", "8", "44", "64", "86", "88", "122")
10^7 * vol1
#    all      8     44     64     86     88    122 
# 0.3005 0.3487 0.3124 0.4477 0.3519 0.3916 0.3214 

vol2 <- rep(0, 7) # volume of 'naive' estimator of cov(coef)
vol2[1] <- det(fm$naive)
vol2[2] <- det(fm.08$naive)
vol2[3] <- det(fm.44$naive)
vol2[4] <- det(fm.64$naive)
vol2[5] <- det(fm.86$naive)
vol2[6] <- det(fm.88$naive)
vol2[7] <- det(fm.122$naive)
names(vol2) <- c("all", "8", "44", "64", "86", "88", "122")
10^7 * vol2
#    all      8     44     64     86     88    122 
# 0.4884 0.8620 0.8250 0.5724 0.5323 0.7317 0.5571 

## Percentage change in parameter estimates
change <- matrix(0, nrow = 9, ncol = 6)
change[1:6,1] <- 100 * (tab1[3,] -  tab1[1,]) / tab1[1,]
change[1:6,2] <- 100 * (tab1[5,] -  tab1[1,]) / tab1[1,]
change[1:6,3] <- 100 * (tab1[7,] -  tab1[1,]) / tab1[1,]
change[1:6,4] <- 100 * (tab1[9,] -  tab1[1,]) / tab1[1,]
change[1:6,5] <- 100 * (tab1[11,] -  tab1[1,]) / tab1[1,]
change[1:6,6] <- 100 * (tab1[13,] -  tab1[1,]) / tab1[1,]
change[7,] <- 100 * (alpha[-1] -  alpha[1]) / alpha[1]
change[8,] <- 100 * (vol1[-1] -  vol1[1]) / vol1[1]
change[9,] <- 100 * (vol2[-1] -  vol2[1]) / vol2[1]
colnames(change) <- c("8", "44", "64", "86", "88", "122")
rownames(change) <- c(fm$xnames, "alpha", "robust", "naive")
change
#                  8     44     64     86     88    122
# (Intercept) -13.04  10.18  11.50   1.09  14.73   4.11
# gender       44.54   2.00 -32.83  -1.33   2.86   0.33
# age          35.07  15.76 -19.17  -0.73   9.86   3.35
# dayacc       16.86  -0.16   4.24   0.78   8.87  -4.17
# severe      -20.47 -11.07   1.36   0.31  17.67  16.25
# toilet       32.94  93.88  -3.11   0.46   3.90 -10.49
# alpha        -7.40   9.54  10.67 -15.74 -66.31  11.10
# robust       16.06   3.96  49.02  17.11  30.32   6.96
# naive        76.50  68.92  17.21   9.00  49.83  14.07

## Cluster-level gradient distances
bd <- BD.distance(fm, x) # WARNING! lot of messages are displayed
one <- BD.1step(fm, x)
idx <- 1:38
which <- c(3,4,12,25,38)
lab <- c(27,41,107,156,235)
ticklab <- c(8,24,27,41,45,55,56,60,65,89,102,107,108,111,113,118,124,125,127,130,132,137,146,153,156,182,185,195,201,206,207,208,211,216,220,228,232,235)

## Fig 1.a: plot of gradient distances
par(pty = "s")
plot(bd, ylim = c(0,.8), xlab = "Cluster", ylab = "Gradient distances", axes = FALSE, lwd = 2)
box()
axis(2)
axis(1, at = 1:38, lab = as.character(ticklab))
text(idx[which], bd[which], labels = as.character(lab), pos = 3)

## Fig 1.b: plot of gradient distances (one-step approximation)
par(pty = "s")
plot(one, ylim = c(0,.8), xlab = "Cluster", ylab = "Gradient distances, one-step approximation", axes = FALSE, lwd = 2)
box()
axis(2)
axis(1, at = 1:38, lab = as.character(ticklab))
text(idx[which], one[which], labels = as.character(lab), pos = 3)

## Observation-level influence measures
ID <- fm$id
len <- table(ID)
o <- BD.msom(fm, x)
cutoff <- qchisq(0.975, df = 1)

# Fig 2.a: plot of Cook's distances
par(pty = "s")
plot(rep(1:38, len), o$cooks, xlab = "Cluster", ylab = "Cook's distances", ylim = c(0, 1.7), axes = FALSE, lwd = 2)
box()
axis(2)
axis(1, at = 1:38, lab = as.character(ticklab))
text(c(3,12,18), o$cooks[c(8,44,64)], lab = c("8", "44", "64"), pos = 3)

# Fig 2.b: plot of Venezuela's distances
par(pty = "s")
plot(rep(1:38, len), o$venezuelas, xlab = "Cluster", ylab = "Venezuela's distances", ylim = c(0, 1.7), axes = FALSE, lwd = 2)
box()
axis(2)
axis(1, at = 1:38, lab = as.character(ticklab))
text(c(3,12), o$venezuelas[c(8,44)], lab = c("8", "44"), pos = 3)

# Fig 3.a: plot of gradient-type statistics
par(pty = "s")
plot(rep(1:38, len), o$BF, xlab = "Cluster", ylab = "Gradient-type statistics", ylim = c(0,10), axes = FALSE, lwd = 2)
box()
axis(2)
axis(1, at = 1:38, lab = as.character(ticklab))
abline(h = cutoff, col = "gray55", lwd = 2, lty = 2)
text(c(25,27), o$BF[c(86,98)], lab = c("86", "98"), pos = 3)

# Fig 3.b: plot of score-type statistics
par(pty = "s")
plot(rep(1:38, len), o$score, xlab = "Cluster", ylab = "score-type statistics", ylim = c(0,10), axes = FALSE, lwd = 2)
box()
axis(2)
axis(1, at = 1:38, lab = as.character(ticklab))
abline(h = cutoff, col = "gray55", lwd = 2, lty = 2)
text(c(25,27), o$score[c(86,98)], lab = c("86", "98"), pos = 3)
