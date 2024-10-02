##################################################################################
# This code is made to search existing genes id within certain range in RAP-DB
##################################################################################
library("rvest")
library("magrittr")
library("RSelenium")
library("dplyr")

###熱點區間
###snp matrix
abh_geno <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_abhgenotype.csv",header = T)

############叫出bkpt
breakpoint <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv",header = T)


#########建立熱點區間 (加入overlapped interval)
overlap_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_overlap_hotspot.csv",header = T)
poisson_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file//poisson hotspot/gbs_poisson_hotspot.csv",header = T)
local_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/local_recomb_hotspot/gbs_nonoverlap_local_hotspot.csv",header = T)

#names(local_hotspot)[names(local_hotspot) == "map"] <- "chr"
#local_hotspot$chr <- gsub("chr","",local_hotspot$chr)
###綜合成hotspot_interval
hotspot_interval <- rbind(local_hotspot[,c(1,2,3)],poisson_hotspot[,c(1,4,5)],overlap_hotspot[,c(1,2,3)])
hotspot_interval$type <- rep(c("local","poisson","overlap"),times = c(25,93,21))
hotspot_interval$chr <- as.numeric(hotspot_interval$chr)


# 打開瀏覽器
##在命令提示字元 將路徑改成d槽(d:) 輸入java -Dwebdriver.chrome.driver=D:chromedriver.exe -jar selenium-server-standalone-3.141.59.jar
##如果遇到出錯指令:
#打開"C:\Users\lingo1st\Dropbox\其他\動態網頁爬蟲.docx"，裡面有超連結可以重新下載最新版本的chromedriver.exe
remDr <- remoteDriver(
  remoteServerAddr = "localhost",
  port = 4444,
  browserName = "chrome")
###rapdb
remDr$open()
remDr$navigate("https://rapdb.dna.affrc.go.jp/index.html")
path = "//*[@id='header']/ul/li[2]/a/img"
btn <- remDr$findElement(using = "xpath", value = path)
btn$clickElement()
###jbrowse
path2 <-"//*[@id='header']/ul/li[2]/ul/li[1]/a" 
btn2 <- remDr$findElement(using = "xpath", value = path2)
btn2$clickElement()
Sys.sleep(5)

for (i in 1:nrow(hotspot_interval)) {
  
  
  ###input for chr postion
  print("step3")
  path3 <- "//*[@id='location']"
  btn3 <- remDr$findElement(using = "xpath", value = path3)
  btn3$clickElement()
  btn3$clearElement()
  text <- list(paste0("chr",hotspot_interval$chr[i],":",hotspot_interval$start[i],"..",hotspot_interval$end[i]) , key = "enter")
  btn3$sendKeysToElement(text)
  Sys.sleep(4)
  #把gene locus顯示功能打開
  if(i == 1){
    path4 <- "//*[@id='dijit_TitlePane_1_pane']/label[2]/span"
    btn4 <- remDr$findElement(using = "xpath", value = path4)
    btn4$clickElement()
  }
  

  ###點選gene locus下拉選單
    print("step5")
    path5 <- "//*[@id='label_Gene locus']/span"
    btn5 <- remDr$findElement(using = "xpath", value = path5)
    #將滑鼠移動至該按鈕
    remDr$mouseMoveToLocation(webElement = btn5)
    ##點擊右鍵
    remDr$executeScript("arguments[0].dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, cancelable: true, view: window }));", list(btn5))
    #btn5$clickElement()
    
  
    
  ### save track data(需根據element中的Save字串找出按鈕)

  print("step6")
  ##找出aria-label具有save track data文字的tr(表格row)
  path6 <- "//tr[contains(@aria-label, 'Save track data ')]"
  button <- remDr$findElement(using = "xpath", value = path6)
  button$getElementLocation()
  remDr$mouseMoveToLocation(webElement =button)
  
  button$clickElement()
  button$click()
  Sys.sleep(2)
  ###更改格式成bed檔
  print("step7")
  path7 <- "//*[@id='exportDialog']/div[2]/form/fieldset[2]/div[2]"
  btn7 <- remDr$findElement(using = "xpath", value = path7)
  btn7$clickElement()
  Sys.sleep(1)
  ### 下載檔案(需根據element中的Save字串找出按鈕)
  print("step8")
  button_xpath_2 <- "//span[@class='dijitReset dijitInline dijitButtonText' and contains(text(), 'Save')]"
  btn8 <- remDr$findElement(using = "xpath", value = button_xpath_2)
  btn8$clickElement()
  Sys.sleep(2)
  ###重新整理 (非常之重要，沒有重新整理的話，會定位到前面已消失的元素，使程式沒法進行)
  remDr$refresh()
  Sys.sleep(3)
}


################# 將下載的檔案讀取到r並擷取基因代號 ###################################
path_wd <- "C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/gene/poisson"
setwd("C:/Users/lingo1st/Dropbox/林冠瑜/gene id/longer local(10kb)_plus poisson/poisson_gene")
files <- list.files(path_wd)
matching_files <- grep("^Gene locus", files, value = TRUE)
####將基因取出
all_genes <- c()
chr <- c()
gene_amount <- c()
for (i in 1:length(matching_files)) {
  #skip = 1 代表略過第一行
  p <- read.csv(paste0(path_wd,"/",matching_files[i]),header = F, sep="\t") 
  if(ncol(p)>1){
    all_genes <- c(all_genes,p$V4[-1])
    chr <- c(chr,substr(p$V1[2],4,5) )
    gene_amount <- c(gene_amount,(nrow(p)-1) )
  }else{next}
  
}

poisson_gene_df <- data.frame(chr = rep(chr,times = gene_amount),gene = all_genes)

write.csv(all_genes,file = paste0(path_wd,"/","gene id(poisson).txt"),row.names = F,quote = F)


