library(vcfR)
###read sample file
vcf <- read.vcfR(file.choose(), verbose = FALSE )
fix <- getFIX(vcf)
fix <- as.data.frame(fix)
head(fix[1:2,])
###read parent file
vcf_parent <- read.vcfR(file.choose(), verbose = FALSE )
fix[,1] <- as.factor(fix[,1])
fix2 <- getFIX(vcf_parent)
fix2 <-as.data.frame(fix2)
rm(f)
### remove "chrom" from coloumn
fix$CHROM <- gsub("chr","",as.character(fix$CHROM))
fix$CHROM <- as.numeric(fix$CHROM)

fix2$CHROM <- gsub("chr","",as.character(fix2$CHROM))
fix2$CHROM <- as.numeric(fix2$CHROM)

### snp number for each chromosome
len <- c()
for (i in 1:12) {
  len[i] <- nrow(fix[fix$CHROM == i,])
}

####intersection
ind_2 <- list()

for (i in 1:12) {
  a <- fix[fix$CHROM == i,]
  b <- fix2[fix2$CHROM == i,]
  ind_2[[i]] <- c(intersect(a$POS,b$POS))
  
}
###count the length of total matched snps
sum(sapply(ind_2, length)) 

### get the actual position of each snp

pos <- list()
for (i in 3:12) {
  p <- unlist(ind_2[[i]])
  pos[[i]] <- NaN*seq(length(p))
  for (j in 1:length(p)) {
    pos[[i]][j] <- which(fix$CHROM == i & fix$POS == p[j])
  }
}


fix$POS[629731]
pos_test <- unlist(pos)

final <- vcf[pos_test,1:25]
write.vcf(final,file = "final.vcf.gz")

