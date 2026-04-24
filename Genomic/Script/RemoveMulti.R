
# prune multi-allelic sites to one allele
# with different strategies (keep most frequent allele, keep most significant allele)
prune_multi <- function(dt.af,dt.dp=NULL,dt.pval=NULL,
                        drop.af.cols=5L,
                        drop.dp.cols=2L,
                        how="most_frequent") {
  how <- match.arg(how,c("most_frequent","most_significant"))
  stopifnot(all(key(dt.af)[1:2]==c("CHROM","POS")))
  dt.af.key <- key(dt.af)
  dt.af[,N:=.N,by=.(CHROM,POS)]
  fixpos <- dt.af[,N>1]
  dt.af[,N:=NULL]
  setkeyv(switch(how,
                 `most_frequent`={
                   stopifnot(!is.null(dt.dp))
                   stopifnot(all(key(dt.dp)==c("CHROM","POS")))
                   af.mat <- as.matrix(dt.af[fixpos,-(1:drop.af.cols)])
                   dp.mat <- as.matrix(dt.dp[setkey(dt.af[fixpos,.(CHROM,POS)],
                                                    CHROM,POS)][
                                                      ,-(1:drop.dp.cols)])
                   rbind(dt.af[!fixpos],
                         dt.af[fixpos][
                           ,COUNT:=rowSums2(af.mat * dp.mat)][
                             ,.SD[which.max(COUNT)],keyby=.(CHROM,POS)][
                               ,COUNT:=NULL])},
                 `most_significant`={
                   stopifnot(!is.null(dt.pval))
                   stopifnot("PVAL"%in%names(dt.pval))
                   stopifnot(nrow(dt.pval)==nrow(dt.af))
                   rbind(cbind(dt.af[!fixpos],dt.pval[!fixpos,.(PVAL)]),
                         cbind(dt.af[fixpos],dt.pval[fixpos,.(PVAL)])[
                           ,.SD[which.min(PVAL)],keyby=.(CHROM,POS)])}),dt.af.key)
}
