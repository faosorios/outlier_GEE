BD.1step <- function(object, x)
{ # GEE diagnostics for binomial with logit-link
  # one-step approximation to the bilinear form distance
  ID <- object$id
  len <- table(ID)
  n <- length(len)
  N <- nrow(x)

  y <- object$y
  cf <- object$coef
  R <- object$working
  mu <- object$fitted
  res <- object$resid
  rsp <- res / sqrt(mu * (1 - mu)) # Pearson residual
  naive <- object$naive # inverse of the sensitivity matrix
  robust <- object$robust # inverse of the godambe information
  p <- length(cf)
  cnames <- names(cf)
  rm.cf <- matrix(0, nrow = n, ncol = p)
  colnames(rm.cf) <- cnames
  distances <- rep(0, n)

  i <- 1
  l <- 1
  while (l < N) {
    idx <- l:(l + len[i] - 1)
    ss <- 1:len[i]
    if (length(idx) != 1)
      Psi <- t(x[idx,]) %*% diag(sqrt(mu[idx] * (1 - mu[idx]))) %*% solve(R[ss,ss], rsp[idx])
    else
      Psi <- x[idx,] * sqrt(mu[idx] * (1 - mu[idx])) * rsp[idx] / R[ss,ss]
    deleted <- cf - naive %*% Psi
    diff <- cf - deleted
    distances[i] <- crossprod(Psi, naive) %*% solve(robust, diff)
    rm.cf[i,] <- deleted
    l <- l + len[i]
    i <- i + 1
  }

  distances
}

BD.distance <- function(object, x)
{ # GEE diagnostics for binomial with logit-link
  # bilinear form distances using complete iteration of the GEE estimation algorithm
  ID <- object$id
  len <- table(ID)
  n <- length(len)
  N <- nrow(x)
  fID <- as.factor(ID)
  lev <- as.integer(levels(fID))

  y <- object$y
  cf <- object$coef
  R <- object$working
  mu <- object$fitted
  res <- object$resid
  rsp <- res / sqrt(mu * (1 - mu)) # Pearson residual
  naive <- object$naive # inverse of the sensitivity matrix
  robust <- object$robust # inverse of the godambe information
  p <- length(cf)
  cnames <- names(cf)
  rm.cf <- matrix(0, nrow = n, ncol = p)
  colnames(rm.cf) <- cnames
  distances <- rep(0, n)

  i <- 1
  l <- 1
  while (l < N) {
    idx <- l:(l + len[i] - 1)
    ss <- 1:len[i]
    if (length(idx) != 1)
      Psi <- t(x[idx,]) %*% diag(sqrt(mu[idx] * (1 - mu[idx]))) %*% solve(R[ss,ss], rsp[idx])
    else
      Psi <- x[idx,] * sqrt(mu[idx] * (1 - mu[idx])) * rsp[idx] / R[ss,ss]
    rm.obs <- guide$practID != lev[i]
    fm <- gee(bothered ~ gender + age + dayacc + severe + toilet, id = practID, subset = rm.obs, data = guide, family = binomial("logit"), corstr = "exchangeable", scale.fix = TRUE, scale.value = 1.)
    deleted <- fm$coef
    diff <- cf - deleted
    distances[i] <- crossprod(Psi, naive) %*% solve(robust, diff)
    rm.cf[i,] <- deleted
    l <- l + len[i]
    i <- i + 1
  }

  distances
}

BD.msom <- function(object, x)
{ # Bilinear form influence measures for a single observation (i.e 'patients')
  ID <- object$id
  len <- table(ID)
  n <- length(len)
  N <- nrow(x)
  p <- ncol(x)

  y <- object$y
  cf <- object$coef
  R <- object$working
  mu <- object$fitted
  eta <- object$linear
  res <- object$resid

  #
  A <- diag(mu * (1 - mu), N) # denoted D in manuscript
  Omega <- matrix(0, N, N) # Sigma in the manuscript
  invOmega <- matrix(0, N, N)
  res <- rep(0, N) # 'Pearson' residual
  Lambda <- matrix(0, N, N)

  i <- 1
  l <- 1
  while (l < N) {
    idx <- l:(l + len[i] - 1)
    ss <- 1:len[i]
    Omega[idx,idx] <- sqrt(A[idx,idx]) %*% R[ss,ss] %*% sqrt(A[idx,idx])
    invOmega[idx,idx] <- solve(Omega[idx,idx])
    res[idx] <- solve(A[idx,idx], y[idx] - mu[idx])
    Lambda[idx,idx] <- outer(res[idx], res[idx])
    l <- l + len[i]
    i <- i + 1
  }

  #
  W <- A %*% invOmega %*% A
  xx <- crossprod(x, W %*% x)
  mid <- crossprod(x, (W %*% Lambda %*% W) %*% x)

  #
  Smat <- matrix(0, p+1, p+1)
  Vmat <- matrix(0, p+1, p+1)
  Smat[2:(p+1),2:(p+1)] <- xx
  Vmat[2:(p+1),2:(p+1)] <- mid
  gone <- rep(0, N)
  Wald <- rep(0, N)
  score <- rep(0, N)
  BF <- rep(0, N)
  cooks <- rep(0, N)
  venezuelas <- rep(0, N)
  lever <- rep(0, N)
  i <- 1
  j <- 1
  l <- 1
  while (l < N) {
    idx <- l:(l + len[i] - 1)
    ss <- 1:len[i]
    Xi <- as.matrix(x[idx,])
    Wi <- W[idx,idx]
    Li <- Lambda[idx,idx]
    ri <- res[idx]
    if (length(idx) == 1)
      Xi <- t(Xi)
    if (length(idx) != 1)
      Hi <- Xi %*% solve(xx, t(Xi) %*% Wi)
    else {
      lev <- Xi %*% solve(xx, t(Xi))
      Hi <- lev * Wi
    }
    Bi <- rep(0, len[i])
    lever[idx] <- diag(Hi)
    for (k in ss) {
      Bi[k] <- 1
      Id <- diag(len[i])
      Smat[1,1] <- crossprod(Bi, Wi %*% Bi)
      BWx <- crossprod(Bi, Wi %*% Xi)
      Smat[1,2:(p+1)] <- BWx
      Smat[2:(p+1),1] <- t(Smat[1,2:(p+1)])
      Vmat[1,1] <- crossprod(Bi, (Wi %*% Li %*% Wi) %*% Bi)
      Vmat[1,2:(p+1)] <- crossprod(Bi, (Wi %*% Li %*% Wi) %*% Xi)
      Vmat[2:(p+1),1] <- t(Smat[1,2:(p+1)])
      Gmat <- t(chol(Vmat))
      Umat <- solve(Vmat)
      Prod <- Umat %*% Smat
      Jmat <- crossprod(Smat, Prod)
      Kmat <- solve(Jmat)
      Psi1 <- crossprod(Bi, Wi %*% ri)
      if (length(idx) != 1) {
        gone[j] <- Psi1 / (crossprod(Bi, Wi %*% (Id - Hi)) %*% Bi)
      } else {
        rhs <- Wi * ri
        gone[j] <- rhs / c(Wi * (1 - Hi))
      }
      Wald[j] <- (gone[j]^2) / Kmat[1,1]
      score[j] <- (Psi1^2) * Umat[1,1]
      BF[j] <- Psi1 * Prod[1,1] * gone[j]
      cooks[j] <- (gone[j]^2) * c(BWx %*% solve(mid, t(BWx)))
      venezuelas[j] <- (gone[j]^2) * c(BWx %*% solve(xx, t(BWx)))
      j <- j + 1
      Bi[k] <- 0
    }
    l <- l + len[i]
    i <- i + 1
  }

  obj <- list(gamma = gone, Wald = Wald, score = score, BF = BF, cooks = cooks, venezuelas = venezuelas, leverages = lever)
  obj
}
