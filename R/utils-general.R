## Selectors for probes - useful for exon/gene arrays

getFidProbeset <- function(object){
  conn <- db(object)
  sql <- "SELECT fid, fsetid FROM pmfeature"
  featureInfo <- dbGetQuery(conn, sql)
  featureInfo <- featureInfo[order(featureInfo[["fsetid"]]),]
  rownames(featureInfo) <- NULL
  return(featureInfo)
}

getFidMetaProbesetCore <- function(object){
  conn <- db(object)
  sql <- "SELECT fid, meta_fsetid as fsetid FROM pmfeature INNER JOIN core_mps USING(fsetid)"
  featureInfo <- dbGetQuery(conn, sql)
  featureInfo <- featureInfo[order(featureInfo[["fsetid"]]),]
  rownames(featureInfo) <- NULL
  return(featureInfo)
}

getFidMetaProbesetFull <- function(object){
  conn <- db(object)
  sql <- "SELECT fid, meta_fsetid as fsetid FROM pmfeature INNER JOIN full_mps USING(fsetid)"
  featureInfo <- dbGetQuery(conn, sql)
  featureInfo <- featureInfo[order(featureInfo[["fsetid"]]),]
  rownames(featureInfo) <- NULL
  return(featureInfo)
}

getFidMetaProbesetExtended <- function(object){
  conn <- db(object)
  sql <- "SELECT fid, meta_fsetid as fsetid FROM pmfeature INNER JOIN extended_mps USING(fsetid)"
  featureInfo <- dbGetQuery(conn, sql)
  featureInfo <- featureInfo[order(featureInfo[["fsetid"]]),]
  rownames(featureInfo) <- NULL
  return(featureInfo)
}


