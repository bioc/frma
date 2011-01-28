frmaBatch <- function(object, background, normalize, input.vecs, output.param, verbose){

  cdfname <- cleancdfname(cdfName(object))
  platform <- gsub("cdf","",cdfname)
  
  if(background == "rma"){
    if(verbose) message("Background Correcting ...\n")
    object <- bg.correct.rma(object)
  }

  if(is.null(input.vecs$normVec) | is.null(input.vecs$probeVec) | is.null(input.vecs$probeVarWithin) | is.null(input.vecs$probeVarBetween)){
    pkg <- paste(platform, "frmavecs", sep="")
    require(pkg, character.only=TRUE, quietly=TRUE) || stop(paste(pkg, "package must be installed first"))
    data(list=eval(pkg))

    if(is.null(input.vecs$normVec)) input.vecs$normVec <- get(pkg)$normVec
    if(is.null(input.vecs$probeVec)) input.vecs$probeVec <- get(pkg)$probeVec
    if(is.null(input.vecs$probeVarWithin)) input.vecs$probeVarWithin <- get(pkg)$probeVarWithin
    if(is.null(input.vecs$probeVarBetween)) input.vecs$probeVarBetween <- get(pkg)$probeVarBetween
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
  
  colnames(exprs) <- sampleNames(object)
  rownames(exprs) <- unique(pns)
  
  return(list(exprs=exprs, stderr=NULL, weights=weights, residuals=residuals))
}

batchFit <- function(x1, x2, x3, x4){
  y.new <- x1 - x2
  yy.new <- as.vector(t(y.new))
  X <- matrix(rep(diag(ncol(y.new)),nrow(y.new)), nrow=length(yy.new), byrow=TRUE)
  Z <- diag(nrow(y.new))
  Z <- Z[rep(seq(nrow(Z)), each=ncol(y.new)),]
  G <- diag(x4)
  R <- diag(rep(x3, each=ncol(y.new)))
  V <- Z%*%G%*%t(Z) + R
  Vinv <- solve(V)
  e <- eigen(Vinv, symmetric=TRUE)
  C <- (e$vectors)%*%diag(sqrt(e$values))%*%t(e$vectors)
  yy.trans <- C%*%yy.new
  x.trans <- C%*%X

  fit.batch <- rlm(yy.trans ~ -1 + x.trans, maxit=100)
  b <- fit.batch$coef
  names(b) <- NULL
  u <- G%*%t(Z)%*%Vinv%*%(yy.new-X%*%b)
  offset.u <- mean(u)
  exprs <- b + offset.u
  w <- fit.batch$w
  return(list(exprs=exprs, w=w))
}
