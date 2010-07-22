GNUSE <- function(object,type=c("plot","values","stats","density"),...){
  type <- match.arg(type)
  if(length(se(object))==0) stop("Object does not contain standard errors.")

  cdfname <- cleancdfname(cdfName(object))
  platform <- gsub("cdf","",cdfname)
  pkg <- paste(platform, "frmavecs", sep="")
  require(pkg, character.only=TRUE, quiet=TRUE) || stop(paste(pkg, "package must be installed first"))
  data(list=eval(pkg))

  gnuses <- sweep(se(object),1,get(pkg)$medianSE,FUN='/')

  if (type == "values"){
    return(gnuses)
  } else if (type == "density"){
    plotDensity(gnuses, ...)
  } else if (type=="stats"){
    Medians <- apply(gnuses,2,median)
    Quantiles <- apply(gnuses,2,quantile,prob=c(0.25,0.75,0.95,0.99))
    nuse.stats <- rbind(Medians,Quantiles[2,] - Quantiles[1,],Quantiles[3,],Quantiles[4,])
    rownames(nuse.stats) <- c("median","IQR","95%","99%")
    return(nuse.stats)
  }
  if (type == "plot"){
    boxplot(data.frame(gnuses), ...)
  }
}


