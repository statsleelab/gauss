#' FDR Inverse Quantile Transformation (FIQT)
#' 
#' @param z association Z-score vector
#' @param min.p minimum p-value admitted (to avoid zero p-values/adjusted p-values which give troubles with inverse cdf)
#' @return R vector containing the Winner's Curse adjusted Z-scores

fiqt <- function(z, min.p=10^-300){
  pvals<-2*pnorm(abs(z),low=F)
  pvals[pvals<min.p]<- min.p
  adj.pvals<-p.adjust(pvals,method="fdr")
  mu.z<-sign(z)*qnorm(adj.pvals/2,low=F)
  mu.z[abs(z)>qnorm(min.p/2,low=F)]<-z[abs(z)>qnorm(min.p/2,low=F)]
  mu.z
}