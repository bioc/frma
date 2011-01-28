GNUSE <- function(object,medianSE=NULL,type=c("plot","values","stats","density"),...){
  type <- match.arg(type)
  if(length(se.exprs(object))==0) stop("Object does not contain standard errors.")

  if(is.null(medianSE)){
    platform <- annotation(object)
    pkg <- paste(platform, "frmavecs", sep="")
    require(pkg, character.only=TRUE, quiet=TRUE) || stop(paste(pkg, "package must be installed first"))
    data(list=eval(pkg))
    medianSE <- get(pkg)$medianSE
  }
  
  gnuses <- sweep(se.exprs(object),1,medianSE,FUN='/')
  if(any(is.na(gnuses))) message("Some GNUSE values are NA. This is often due to probesets with only 1 probe.")
      
  if (type == "values"){
    return(gnuses)
  } else if (type == "density"){
    plotDensity(gnuses, ...)
  } else if (type=="stats"){
    Medians <- apply(gnuses,2,median, na.rm=TRUE)
    Quantiles <- apply(gnuses,2,quantile,prob=c(0.25,0.75,0.95,0.99), na.rm=TRUE)
    nuse.stats <- rbind(Medians,Quantiles[2,] - Quantiles[1,],Quantiles[3,],Quantiles[4,])
    rownames(nuse.stats) <- c("median","IQR","95%","99%")
    return(nuse.stats)
  }
  if (type == "plot"){
    boxplot(data.frame(gnuses), ...)
  }
}





