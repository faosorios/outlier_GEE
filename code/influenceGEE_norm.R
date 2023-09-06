
cooks.GEEnorm <- function(object) 
{ ## Cook's distances for single observations in GEE

  # covariance estimates
  naive <- object$naive
  robust <- object$robust

  # extracting info
  nobs <- length(object$y)
  p <- ncol(naive)
  ID <- object$id
  len <- table(ID)

  # 
  aug <- rep(0, nobs)
  aug <- as.data.frame(aug)

  # result containers
  cf <- matrix(0, nrow = nobs, ncol = p)
  cooks <- rep(0, nobs)
  venezuelas <- rep(0, nobs)
  gammas <- rep(0, nobs)

  for (i in 1:nobs) {
    aug[i,] <- 1
    db <- cbind(Lima, aug) # USING specific dataset from Appendix E
    faug <- gee(log(error) ~ block + aug, id = Subject, data = db, family = gaussian, corstr = "AR-M", Mv = 1)
    # saving coefficient without ij-observation and gamma estimates
    cf[i,] <- faug$coef[1:p]
    gammas[i] <- faug$coef[p + 1]
    # influence measures
    diff <- object$coef - cf[i,]
    cooks[i] <- sum(diff * c(solve(robust, diff)))
    venezuelas[i] <- sum(diff * c(solve(naive, diff)))
    # 
    aug[,] <- 0
  }

  # output object
  obj <- list(COEF = cf, gammas = gammas, cooks = cooks, venezuelas = venezuelas, lengths = len)
  obj
}

MSOM.GEEnorm <- function(object, x)
{ ## Gradient- and score-type test statistic for outlier detection

  # extracting info
  x <- as.matrix(x)
  nobs <- nrow(x)
  p <- ncol(x)
  ID <- object$id
  len <- table(ID)
  R <- object$working
  res <- object$resid
  scale <- object$scale
  
  #
  p1 <- p + 1
  W <- matrix(0, nobs, nobs)
  Lambda <- matrix(0, nobs, nobs)

  i <- 1
  l <- 1
  while (l < nobs) {
    idx <- l:(l + len[i] - 1)
    ss <- 1:len[i]
    W[idx,idx] <- solve(R[ss,ss])
    Lambda[idx,idx] <- outer(res[idx], res[idx])
    l <- l + len[i]
    i <- i + 1
  }
  xx <- crossprod(x, W %*% x)
  mid <- crossprod(x, (W %*% Lambda %*% W) %*% x)

  # result containers
  gone <- rep(0, nobs)
  scores <- rep(0, nobs)
  GT <- rep(0, nobs)

  #
  i <- 1
  j <- 1
  l <- 1
  while (l < nobs) {
    idx <- l:(l + len[i] - 1)
    ss <- 1:len[i]
    Xi <- as.matrix(x[idx,])
    Wi <- W[idx,idx]
    Li <- Lambda[idx,idx]
    ri <- res[idx]
    Hi <- Xi %*% solve(xx, t(Xi) %*% Wi)
    Id <- diag(len[i])
    LW <- Li %*% Wi
    Pi <- Id - LW %*% Xi %*% solve(mid, crossprod(Xi, Wi))
    Mi <- Wi %*% Pi %*% LW
    B <- rep(0, len[i])
    for (k in ss) {
      B[k] <- 1
      Psi1 <- crossprod(B, Wi %*% ri)
      prod <- crossprod(B, Wi %*% Pi) %*% B
      div <- crossprod(B, Wi %*% (Id - Hi) %*% B)
      gone[j] <- Psi1 / div
      div <- crossprod(B, Mi %*% B)
      scores[j] <- Psi1^2 / div
      GT[j] <- (Psi1 * prod * gone[j]) / div
      #
      j <- j + 1
      B[k] <- 0
    }
    l <- l + len[i]
    i <- i + 1
  }

  obj <- list(gamma = gone, scores = scores, gradient = GT)
  obj
}

influenceGEE.norm <- function(object, x)
{ # influence measures for single observations in GEE
  ID <- object$id
  len <- table(ID)
  nID <- length(len)

  naive <- object$naive # inverse of the sensitivity matrix
  robust <- object$robust # inverse of the godambe information

  x <- as.matrix(x)
  nobs <- nrow(x)
  p <- ncol(x)
  p1 <- p + 1

  y <- object$y
  R <- object$working
  mu <- object$fitted
  eta <- object$linear
  res <- object$resid
  scale <- object$scale

  W <- matrix(0, nobs, nobs)
  Lambda <- matrix(0, nobs, nobs)

  i <- 1
  l <- 1
  while (l < nobs) {
    idx <- l:(l + len[i] - 1)
    ss <- 1:len[i]
    W[idx,idx] <- solve(R[ss,ss])
    Lambda[idx,idx] <- W[idx,idx] %*% outer(res[idx], res[idx]) %*% W[idx,idx]
    l <- l + len[i]
    i <- i + 1
  }
  xx <- crossprod(x, W %*% x)
  mid <- crossprod(x, Lambda %*% x)

  aug <- rep(0, nobs)
  aug <- as.data.frame(aug)
  cf <- matrix(0, nrow = nobs, ncol = p1)
  
  cooks <- rep(0, nobs)
  venezuelas <- rep(0, nobs)
  gammas <- rep(0, nobs)
  Wald <- rep(0, nobs)
  score <- rep(0, nobs)
  GT <- rep(0, nobs)

  S <- matrix(0, nrow = p1, ncol = p1)
  V <- matrix(0, nrow = p1, ncol = p1)

  S[1:p,1:p] <- xx
  V[1:p,1:p] <- mid

  u <- W %*% res

  for (i in 1:nobs) {
    aug[i,] <- 1
    db <- cbind(Lima, aug) # USING specific dataset from Appendix E
    faug <- gee(log(error) ~ block + aug, id = Subject, data = db, family = gaussian, corstr = "AR-M", Mv = 1)
    # saving coefficient without ij-observation and gamma estimates
    cf[i,] <- faug$coef
    gammas[i] <- cf[i,3]
    # influence measures
    diff <- object$coef - cf[i,1:2]
    cooks[i] <- sum(diff * c(solve(robust, diff)))
    venezuelas[i] <- sum(diff * c(solve(naive, diff)))
    #
    P1 <- c(aug[,1] %*% u)
    B <- as.matrix(aug)
    # 
    V[p1,p1] <- crossprod(B, Lambda %*% B)
    V[p1,1:p] <- crossprod(B, Lambda %*% x)
    V[1:p,p1] <- t(V[p1,1:p])
    U <- solve(V)
    S[p1,p1] <- crossprod(B, W %*% B)
    S[p1,1:p] <- crossprod(B, W %*% x)
    S[1:p,p1] <- t(S[p1,1:p])
    Q <- U %*% S
    # test statistics
    Wald[i] <- gammas[i]^2 / faug$naive[3,3]
    score[i] <- (P1^2) / U[3,3]
    GT[i] <- (P1 * gammas[i]) / Q[3,3]
    # 
    aug[,] <- 0
  }

  # output object
  obj <- list(COEF = cf[,1:2], gammas = gammas, cooks = cooks, venezuelas = venezuelas, Wald = Wald, score = score, GT = GT)
  obj
}
