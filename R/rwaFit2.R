rwaFit2 <- function(x1, x2, x3, x4=NULL){
  ncols <- ncol(x1)
  w.tmp <- x2/max(x2)
  w.tmp <- matrix(rep(w.tmp, ncols), ncol=ncols)
  pe.tmp <- x3
  pe.tmp[1] <- pe.tmp[1]-sum(pe.tmp)
  rcModelWPLM(y=x1, w=w.tmp, row.effects=pe.tmp, input.scale=x4)
}