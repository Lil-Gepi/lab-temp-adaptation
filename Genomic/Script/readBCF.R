# read allele frequencies (what=AF) or coverage depths (what=DP)
# from a given variant file using bcftools-query
read_bcf <- function(vcf,what) {
  what <- match.arg(what,c("AF","XF","DP"))
  switch(what,
         AF={
           cols <- as.character(fread(cmd=paste("bcftools view -h --threads 8",
                                                vcf,
                                                "| tail -1 | cut -f1,2,4,5,8,10- | sed 's/^#//'"),
                                      header=FALSE))
           cols[5] <- "N_ALT"
           dt.af <- fread(cmd=paste("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%N_ALT[\t%AF]\n'",vcf),
                          header=FALSE,col.names=cols,key=c("CHROM","POS","ALT"), nThread = 8,stringsAsFactors = F)
           setnames(dt.af,names(dt.af)[-(1:5)],sub("^.*_F","",names(dt.af)[-(1:5)]))
           },
         XF={
           cols <- as.character(fread(cmd=paste("bcftools view -h --threads 8",
                                                vcf,
                                                "| tail -1 | cut -f1,2,4,5,8,10- | sed 's/^#//'"),
                                      header=FALSE))
           cols[5] <- "N_ALT"
           dt.af <- fread(cmd=paste("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%N_ALT[\t%XF]\n'",vcf),
                          header=FALSE,col.names=cols,key=c("CHROM","POS","ALT"), nThread = 8,stringsAsFactors = F)
           setnames(dt.af,names(dt.af)[-(1:5)],sub("^.*_F","",names(dt.af)[-(1:5)]))
         },
         DP={
           cols <- as.character(fread(cmd=paste("bcftools view -h --threads 8",
                                                vcf,
                                                "| tail -1 | cut -f1,2,10- | sed 's/^#//'"),
                                      header=FALSE))
           dt.dp <- fread(cmd=paste("bcftools query -f '%CHROM\t%POS[\t%SAD]\n'",vcf,"| uniq"),
                          header=FALSE,col.names=cols,key=c("CHROM","POS"), nThread = 8,stringsAsFactors = F)
           setnames(dt.dp,
                    names(dt.dp)[-(1:2)],
                    sub("^.*_F","",names(dt.dp)[-(1:2)]))
         })
}