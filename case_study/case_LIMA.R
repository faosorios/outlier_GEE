## loading dataset and reading R sources
load("../data/Lima.rda")
source("../code/influenceGEE_norm.R")
library(gee)

# Figure E1
library(lattice)
bwtheme <- standard.theme("pdf", color = FALSE)
xyplot(log(error) ~ block|Subject, data = Lima, par.settings = bwtheme, layout = c(10,4), type = "b", lwd = 2, cex = 1)

## fitted model
fm <- gee(log(error) ~ block, id = Subject, data = Lima, family = gaussian, na.action = na.omit, corstr = "AR-M", Mv = 1)
# next four lines are required for our influence measures...
x <- cbind(1, Lima$block)
colnames(x) <- fm$xnames
x <- as.matrix(x)
nobs <- nrow(x)

rm.11 <- rm.51 <- rm.55 <- rm.331 <- rm.391 <- rep(TRUE, nobs)
rm.11[1] <- rm.51[33] <- rm.55[37] <- rm.331[257] <- rm.391[305] <- FALSE
fm11 <- gee(log(error) ~ block, id = Subject, data = Lima, subset = rm.11, family = gaussian, na.action = na.omit, corstr = "AR-M", Mv = 1)
fm51 <- gee(log(error) ~ block, id = Subject, data = Lima, subset = rm.51, family = gaussian, na.action = na.omit, corstr = "AR-M", Mv = 1)
fm55 <- gee(log(error) ~ block, id = Subject, data = Lima, subset = rm.55, family = gaussian, na.action = na.omit, corstr = "AR-M", Mv = 1)
fm331 <- gee(log(error) ~ block, id = Subject, data = Lima, subset = rm.331, family = gaussian, na.action = na.omit, corstr = "AR-M", Mv = 1)
fm391 <- gee(log(error) ~ block, id = Subject, data = Lima, subset = rm.391, family = gaussian, na.action = na.omit, corstr = "AR-M", Mv = 1)

## Table E1 (slightly recrafted)
tab1 <- matrix(0, nrow = 18, ncol = 2)
tab1[1,] <- fm$coef
tab1[2,] <- sqrt(diag(fm$robust))
tab1[3,] <- 2 * pnorm(abs(summary(fm)$coef[,5]), lower.tail = FALSE)
tab1[4,] <- fm11$coef
tab1[5,] <- sqrt(diag(fm11$robust))
tab1[6,] <- 2 * pnorm(abs(summary(fm11)$coef[,5]), lower.tail = FALSE)
tab1[7,] <- fm51$coef
tab1[8,] <- sqrt(diag(fm51$robust))
tab1[9,] <- 2 * pnorm(abs(summary(fm51)$coef[,5]), lower.tail = FALSE)
tab1[10,] <- fm55$coef
tab1[11,] <- sqrt(diag(fm55$robust))
tab1[12,] <- 2 * pnorm(abs(summary(fm55)$coef[,5]), lower.tail = FALSE)
tab1[13,] <- fm331$coef
tab1[14,] <- sqrt(diag(fm331$robust))
tab1[15,] <- 2 * pnorm(abs(summary(fm331)$coef[,5]), lower.tail = FALSE)
tab1[16,] <- fm391$coef
tab1[17,] <- sqrt(diag(fm391$robust))
tab1[18,] <- 2 * pnorm(abs(summary(fm391)$coef[,5]), lower.tail = FALSE)
colnames(tab1) <- fm$xnames
rownames(tab1) <- c("all", "SE", "p-value", "1,1", "SE", "p-value", "5,1", "SE", "p-value", 
                    "5,5", "SE", "p-value", "33,1", "SE", "p-value", "39,1", "SE", "p-value")
tab1
#        (Intercept)   block
#all          3.8500 -0.0513
#SE           0.0665  0.0101
#p-value      0.0000  0.0000
#1,1          3.8299 -0.0481
#SE           0.0644  0.0095
#p-value      0.0000  0.0000
#5,1          3.8317 -0.0483
#SE           0.0694  0.0104
#p-value      0.0000  0.0000
#5,5          3.8519 -0.0511
#SE           0.0667  0.0101
#p-value      0.0000  0.0000
#33,1         3.8335 -0.0487
#SE           0.0613  0.0094
#p-value      0.0000  0.0000
#39,1         3.8349 -0.0489
#SE           0.0640  0.0099
#p-value      0.0000  0.0000

# row: alpha
alpha <- rep(0, 6)
alpha[1] <- fm$work[2,1]
alpha[2] <- fm11$work[2,1]
alpha[3] <- fm51$work[2,1]
alpha[4] <- fm55$work[2,1]
alpha[5] <- fm331$work[2,1]
alpha[6] <- fm391$work[2,1]
alpha
# 0.5309 0.5400 0.5440 0.5342 0.5269 0.5310
100 * (alpha[-1] -  alpha[1]) / alpha[1]
# 1.7075 2.4567 0.6089 -0.7620 0.0045

# row: scale
phi <- rep(0, 6)
phi[1] <- fm$scale
phi[2] <- fm11$scale
phi[3] <- fm51$scale
phi[4] <- fm55$scale
phi[5] <- fm331$scale
phi[6] <- fm391$scale
phi
# 0.1730 0.1687 0.1720 0.1707 0.1681 0.1698
100 * (phi[-1] -  phi[1]) / phi[1]
# -2.4652 -0.5829 -1.3253 -2.8430 -1.8434

## Observation-level influence measures
z <- cooks.GEEnorm(fm)

# Fig E2.a: plot of Cook's distances
obs <- (1:nobs)[z$cooks > 0.06]
par(pty = "s")
plot(rep(1:40, z$len), z$cooks, xlab = "Cluster", ylab = "Cook's distances", ylim = c(0, 0.12), lwd = 2, cex.lab = 1.3)
text(rep(1:40, z$len)[obs], z$cooks[obs], lab = c("1,1", "5,1", "33,1", "39,1"), pos = 3)

# Fig E2.b: plot of Venezuela's distances
obs <- (1:nobs)[z$venez > 0.05]
par(pty = "s")
plot(rep(1:40, z$len), z$venez, xlab = "Cluster", ylab = "Venezuela's distances", ylim = c(0, 0.12), lwd = 2, cex.lab = 1.3)
text(rep(1:40, z$len)[obs], z$venez[obs], lab = c("1,1", "5,1", "33,1"), pos = 3)

# Fig E3:
o <- MSOM.GEEnorm(fm, x)
cutoff <- qchisq(0.975, df = 1)
par(pty = "s")
plot(rep(1:40, z$len), o$gradient, xlab = "Cluster", ylab = "gradient-type statistics", ylim = c(1, 2), lwd = 2, cex.lab = 1.3)
