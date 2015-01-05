batchFit <- function(y, probe, sigma, tau, maxiter=100){
  nchip <- ncol(y)
  nprobe <- nrow(y)
  y.new <- as.vector(y - probe)  #kill any names
  n <- length(y)
  wt <- matrix(rep(1/as.vector(sigma), nchip), nrow=nprobe)
  
  #rlm setup, taken from the rlm function
  irls.delta <- function(old, new)
    sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
  psi.huber <- function(u, k=1.145)  pmin(1, k/abs(u))
  
  # Setup for linear regression
  xtx <- matrix(0., ncol=(nchip+nprobe), nrow=nchip+nprobe)
  i1 <- seq(1, length=nchip+nprobe, by=nchip+nprobe+1) #index of diagonal
  i2a <- which(row(xtx) <= nprobe & col(xtx) >nprobe)
  i2b <- c(outer(seq(nprobe, length=nprobe, by=nrow(xtx)), 1:nchip,  '+'))
  i3 <- as.vector(row(y))          #gamma coefs
  i4 <- as.vector(col(y)) + nprobe #theta coefs
  sqwt <- sqrt(wt) 
  w <- rep(1., length(y))
  
  # Initial linear regression
  xtx[i1] <- c(rowSums(wt) + 1/tau, colSums(wt))
  xtx[i2a] <- wt
  xtx[i2b] <- wt         
  xty <- c(rowSums(wt*y.new), colSums(wt*y.new))
  
  beta <- solve(xtx, xty)
  resid <- sqwt * (beta[i3] + beta[i4] - y.new)  # weighted residuals
  
  for (iter in 1:maxiter) {
    oldresid <- resid
    scale <- median(abs(resid))/.6745
    w <- matrix(psi.huber(resid/scale, k=1.345), nrow=nprobe)
    w2 <- w*wt
    xtx[i1] <- c(rowSums(w2) + 1/tau, colSums(w2))
    xtx[i2a] <- w2
    xtx[i2b] <- w2
    temp <- w2* y.new
    xty <- c(rowSums(temp), colSums(temp))
    
    beta <- solve(xtx, xty)
    resid <- sqwt * (beta[i3] + beta[i4] - y.new)  # weighted residuals
    convi <- irls.delta(oldresid, resid)
    if (convi < 1e-4) break
  }
  return(list(Estimates=beta[-(1:nprobe)], Weights=as.vector(t(w)), Gamma=beta[1:nprobe]))
}