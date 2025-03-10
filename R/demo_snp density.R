library(vcfR)
###########
largest_density_size <- c()
largest_density_pos <- c()
b <- c()
for (i in 1:12) {
  data <- read.vcfR(paste0("C:/Users/lingo1st/Onedrive/桌面/NOISYmputer/NOISYmputer_data/gbs.chr",i,".vcf"))
  data_fix <- data@fix
  data_fix <- as.data.frame(data_fix)
  data_fix$POS <- as.numeric(data_fix$POS)
  
  for (j in seq(1,max(data_fix$POS),100000)) {
    sub_data <- data_fix[data_fix$POS<= (j+99999) & data_fix$POS >= j,]
    b <- append(b,as.numeric(nrow(sub_data)))
  }
  largest_density_pos[i] <- c(100000*which.max(b))
  largest_density_size[i] <- max(b)
}

####每條染色體密度最高的區間(0.1Mb)，以及裡面的marker數量
snp_den <- data.frame(chr = as.character(c(1:12)), size = largest_density_size, pos = largest_density_pos)

####計算平均snp密度
chr_len <- c(43260640,35954074,36189985,35489479,29733216,30731386,29643843,28434680,22692709,22683701,28357783,27561960)
snp_amount <- c(6160,4695,3630,3948,3091,1532,3479,3201,2543,3027,3749,2282)
round(sum(chr_len)/sum(snp_amount),1)

###個別密度
p <- which.max(chr_len/snp_amount)
q <- which.min(chr_len/snp_amount)
(chr_len/snp_amount)[q]
###snp之間距離
for (j in 1:12 ) {
  data <- read.vcfR(paste0("C:/Users/lingo1st/Onedrive/桌面/NOISYmputer/NOISYmputer_data/gbs.chr",j,".vcf"),verbose = F)
  data_fix <- data@fix
  data_fix <- as.data.frame(data_fix)
  data_fix$POS <- as.numeric(data_fix$POS)
  distance <- c()
  for (i in 1:nrow(data_fix)-1) {
    distance[i] <- data_fix$POS[i+1] - data_fix$POS[i]
  }
  cat("\n","chr",j," pos: ",data_fix$POS[which.max((distance))],"~",data_fix$POS[which.max((distance))+1],"\n",
      "chr",j,"distance: ",max(distance),"\n")
}

#########################
