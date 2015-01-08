barcode <- function(object, platform=NULL, mu=NULL, tau=NULL, cutoff=6.5, output="binary"){

  if(!class(object) %in% c("matrix", "ExpressionSet", "frmaExpressionSet")
     & !is.vector(object)){
      stop("Object must be one of: vector, matrix, ExpressionSet, or frmaExpressionSet.")
  }

  if(!output %in% c("binary", "z-score", "p-value", "lod")){
      stop("Output must be one of: binary, z-score, p-value, lod.")
  }

  if(is.vector(object)) object <- as.matrix(object)
  
  if(is.matrix(object)){
    if(is.null(platform) & (is.null(mu) | is.null(tau))){
        stop("If object is of class matrix or vector, platform cannot be NULL if mu or tau is NULL.")
    }
  }

  if(class(object) %in% c("ExpressionSet", "frmaExpressionSet")){
    if(is.null(platform)){
      if(annotation(object) %in% c("hgu133a", "hgu133a2", "hgu133plus2", "mouse4302", "hugene.1.0.st.v1", "mogene.1.0.st.v1")){
        if(cleancdfname(annotation(object)) == "hgu133acdf") platform <- "GPL96"
        if(cleancdfname(annotation(object)) == "hgu133a2cdf") platform <- "GPL571"
        if(cleancdfname(annotation(object)) == "hgu133plus2cdf") platform <- "GPL570"
        if(cleancdfname(annotation(object)) == "mouse4302cdf") platform <- "GPL1261"
        if(annotation(object) == "hugene.1.0.st.v1") platform <- "GPL6244"
        if(annotation(object) == "mogene.1.0.st.v1") platform <- "GPL6246"
      } else stop("Microarray platform given by object annotation not recognized -- please supply platform type.")
    } else {
      if(!platform %in% c("GPL96", "GPL570", "GPL571", "GPL1261", "GPL6244", "GPL6246")){
          stop("Platform must be one of: GPL96, GPL570, GPL571, GPL1261", "GPL6244", "GPL6246")
      }
    }
    object <- as.matrix(exprs(object))
  }

  if(is.null(mu) | is.null(tau)){
    if(!platform %in% c("GPL96", "GPL570", "GPL571", "GPL1261", "GPL6244", "GPL6246")){
        stop("mu and tau must be non-NULL unless platform is one of: GPL96, GPL570, GPL571, GPL1261, GPL6244, GPL6246")
    }
    if(platform=="GPL96") pkg <- "hgu133afrmavecs"
    if(platform=="GPL570") pkg <- "hgu133plus2frmavecs"
    if(platform=="GPL571") pkg <- "hgu133a2frmavecs"
    if(platform=="GPL1261") pkg <- "mouse4302frmavecs"
    if(platform=="GPL6244") pkg <- "hugene.1.0.st.v1frmavecs"
    if(platform=="GPL6246") pkg <- "mogene.1.0.st.v1frmavecs"
    require(pkg, character.only=TRUE, quietly=TRUE) || stop(paste(pkg, "package must be installed first"))
    data(list=gsub("frma", "barcode", pkg))
    bcparams <- get(gsub("frma", "barcode", pkg))
    mu <- bcparams$mu
    tau <- bcparams$tau
    names(mu) <- names(tau) <- names(bcparams$entropy)
  } else{
    if(length(mu)!=length(tau)) stop("Lengths of mu and tau must be equal.")
  }

  if(is.null(rownames(object))){
    warning("Object does not contain rownames; therefore, matching between object and barcode parameters is not guarenteed.")
    if(nrow(object)!=length(mu)) stop("Number of rows of object does not equal length of barcode parameters (mu and tau).")
  }

  if(is.null(names(mu)) | is.null(names(tau))){
    warning("Either mu or tau does not contain names; therefore, matching between object and barcode parameters is not guarenteed.")
  } else{
    if(!identical(names(mu),names(tau))) stop("Names of mu and tau must be identical.")
  }

  if(nrow(object)!=length(mu)){
      if(platform %in% c("GPL6244", "GPL6246")) message("The following error can likely be fixed by setting target=\"core\" in the frma function.")
      stop("Number of rows / features of object does not equal length of barcode parameters (mu and tau).")
  }

  if((!is.null(rownames(object))) & (!is.null(names(mu))) & (!is.null(names(tau)))){
      map <- match(rownames(object),names(mu))
      mu <- mu[map]
      tau <- tau[map]
      if(!identical(rownames(object), names(mu))){
          stop("Mismatch between rownames of data object and names of mu and unable to create unique mapping.")
      }
  }
  
  if(output %in% c("p-value", "lod", "binary")){
    pval <- pnorm(object, mean=mu, sd=tau, lower.tail=FALSE)
    if(output == "p-value"){
      colnames(pval) <- colnames(object)
      rownames(pval) <- rownames(object)
      return(pval)
    } else{
      lod <- -log10(pval)
      if(output == "lod"){
        colnames(lod) <- colnames(object)
        rownames(lod) <- rownames(object)
        return(lod)
      } else{
        bc <- matrix(as.integer(lod > cutoff), ncol=ncol(object))
        colnames(bc) <- colnames(object)
        rownames(bc) <- rownames(object)
        return(bc)
      }
    }
  }

  if(output == "z-score"){
    ##z <- sweep(sweep(e, 1, mu), 1, tau, FUN="/")
    z <- (object-mu)/tau
    colnames(z) <- colnames(object)
    rownames(z) <- rownames(object)
    return(z)
  }
}
