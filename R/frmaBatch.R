frmaBatch <- function(object, background, normalize, input.vecs, output.param, verbose){

  cdfname <- cleancdfname(cdfName(object))
  tmp <- gsub("cdf","",cdfname)
  platform <- gsub("..entrezg", "", tmp)
  vecdataname <- paste(tmp, "frmavecs", sep="")
  
  if(background == "rma"){
    if(verbose) message("Background Correcting ...\n")
    object <- bg.correct.rma(object)
  }

  if(is.null(input.vecs$normVec) | is.null(input.vecs$probeVec) | is.null(input.vecs$probeVarWithin) | is.null(input.vecs$probeVarBetween)){
    pkg <- paste(platform, "frmavecs", sep="")
    require(pkg, character.only=TRUE, quietly=TRUE) || stop(paste(pkg, "package must be installed first"))
    data(list=eval(vecdataname))

    if(is.null(input.vecs$normVec)) input.vecs$normVec <- get(vecdataname)$normVec
    if(is.null(input.vecs$probeVec)) input.vecs$probeVec <- get(vecdataname)$probeVec
    if(is.null(input.vecs$probeVarWithin)) input.vecs$probeVarWithin <- get(vecdataname)$probeVarWithin
    if(is.null(input.vecs$probeVarBetween)) input.vecs$probeVarBetween <- get(vecdataname)$probeVarBetween
  }

  if(normalize == "quantile"){
    if(verbose) message("Normalizing ...\n")
    pm(object) <- normalize.quantiles.use.target(pm(object), input.vecs$normVec)
  }

  if(verbose) message("Summarizing ...\n")
  pns <- probeNames(object)
  pms <- log2(pm(object))

  N <- 1:dim(pms)[1]
  S <- split(N, pns)

  fit <- lapply(1:length(S), function(i) {
    s <- S[[i]]
    batchFit(pms[s,], input.vecs$probeVec[s], input.vecs$probeVarWithin[s], input.vecs$probeVarBetween[s])
  })
  
  exprs <- matrix(unlist(lapply(fit, function(x) x$exprs)), ncol=ncol(pms), byrow=TRUE)
  rownames(exprs) <- names(fit)
  colnames(exprs) <- colnames(pms)

  if("weights" %in% output.param){
    weights <- matrix(unlist(lapply(fit, function(x) x$w)), ncol=ncol(pms), byrow=TRUE)
    rownames(weights) <- pns
    colnames(weights) <- sampleNames(object)
  } else weights <- NULL

  if("residuals" %in% output.param){
    residuals <- apply(exprs,2,function(x) rep(x, table(pns)))
    residuals <- (pms-input.vecs$probeVec) - residuals
    rownames(residuals) <- pns
    colnames(residuals) <- sampleNames(object)
  } else residuals <- NULL

  if("random_effects" %in% output.param){
    gammas <- matrix(unlist(lapply(fit, function(x) x$gamma)), ncol=1)
    rownames(gammas) <- pns
    colnames(gammas) <- NULL
  } else gammas <- NULL
  
  colnames(exprs) <- sampleNames(object)
  rownames(exprs) <- unique(pns)
  
  return(list(exprs=exprs, stderr=NULL, weights=weights, residuals=residuals, gammas=gammas))
}

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
  return(list(exprs= beta[-(1:nprobe)], w=as.vector(t(w)), gamma=beta[1:nprobe]))
}
