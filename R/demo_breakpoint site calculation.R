#################################################################################################
#This code is made for determining the location of genome-wide recombination breakpoint
#site using ABH form csv file. 
#Input: a ABH form file
#Output:  dataframe with chromosome, breakpoint position, left/right-flanking markers
################################################################################################

library(dplyr)
library(tidyverse)
library(vcfR)
getwd()
###使用abh_genotype.csv 計算 
data <-  read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/gbs_imputation_all individual/all_gbs_abhgenotype.csv",header = T)
breakpoint <- data.frame(Chr = NA,Pos = NA,left = NA,right = NA)
col_name <- rownames(data)
#i為染色體, j為snp數量
for (i in sprintf("%02d", c(1:12))) {
  ##選出該條染色體之子集-->dat
  dat <- data %>%
    select(matches(paste0("X.",i, "_")))
  ###使用apply，針對每個row(個體)計算co發生位置，並記錄在co_site中
  result <- apply(dat,1,function(x){
    co_site <- c()
    co_left <- c()
    co_right <- c()
    start <- ifelse(x[1] == "-",x[2],x[1])
    for (j in 2:length(x)) {
      if( x[j] == start | x[j] == "-" ){
        next
      }else{
        co_site <- c(co_site,(as.numeric(strsplit(names(x[j]),"_")[[1]][2]) + as.numeric(as.numeric(strsplit(names(x[j-1]),"_")[[1]][2])) )/2  )
        co_left <- c(co_left,as.numeric(strsplit(names(x[j-1]),"_")[[1]][2]) )
        co_right <- c(co_right,as.numeric(as.numeric(strsplit(names(x[j]),"_")[[1]][2])) )
        start <- x[j]
        
      }
    }
    return(c(co_site,co_left,co_right))
  })
  site_all <- c()
  site_left<- c()
  site_right <- c()
  for (k in 1:length(result)) {
    if(length(result[[k]])==0){
      next
    }
    temp_index <- seq(1,length(result[[k]]),by = length(result[[k]])/3)
    site_all <- c(site_all,result[[k]][ temp_index[1]:(temp_index[2]-1) ]) 
    site_left <- c(site_left,result[[k]][ temp_index[2]:(temp_index[3]-1 )])
    site_right <- c(site_right,result[[k]][ temp_index[3]:length(result[[k]]) ])
  }
  
  result_df <- data.frame(Chr = rep(i,length(site_all)),Pos = site_all,left = site_left,right = site_right)
  breakpoint <- rbind(breakpoint,result_df)
}
breakpoint <- breakpoint[-1,]
write.csv(breakpoint,file = "wgs_breakpoint.csv",row.names = F)

#######################visualization ####################################
p1 <- ggplot(data = breakpoint,aes(x = Chr,y=Pos/1000000))
for (i in 1:12) {
  p1<- p1+
    geom_point(size = 0.01,position = position_dodge(width = 2))
  #geom_segment(x = i , y = 1, xend = i, yend = chr_len[i]/1000000)
}
p1 <- p1 + ylab("Position (Mb)")+xlab("Chromosome") +ggtitle("WGS")
####try to add violin plot
p1 <- p1+ geom_violin(aes(x = Chr,y = Pos/1000000,group = Chr),fill = "orange",alpha = 0.2)
####add centromere region
cent_min <- c(15.445668,12.625206,17.883934,7.882,11.15,13.201,9.104,11.98,0.993,
              7.62,11.34,11.06)
cent_max <- c(18.05,15.48,20.51,10.06,13.54,17.84,12.71,14.41,3.93,8.69,13.5,
              12.2)

cent <- data.frame(xmin = seq(1,12,by = 1)-0.2,xmax = seq(1,12,by = 1)+0.2,ymin = cent_min,ymax = cent_max)
for (i in 1:12) {
  p1<- p1+
    geom_rect(xmin = cent$xmin[i], xmax = cent$xmax[i], ymin = cent$ymin[i], ymax = cent$ymax[i],alpha = 0.5)
}
p1 
