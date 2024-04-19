
remove(list = ls())
setwd("~/project/scCNSignature/ACE_NEW/")


library(data.table)
library(dplyr)
library(stringr)
library(ACE)


filename <- list.files("~/project/scCNSignature/ACE_NEW/pcawg_low/pcawg_QDNAseq/")
path <- "~/project/scCNSignature/ACE_NEW/pcawg_low/pcawg_QDNAseq/"
ACE_path <- "~/project/scCNSignature/ACE_NEW/pcawg_low/ACE_result/"


sc_ACE_para <- data.frame(matrix(data = NA,nrow =160,ncol = 3))
colnames(sc_ACE_para) <- c("sample","ploidy","cellularity")
rownames(sc_ACE_para) <- filename

for(i in filename){
  load(paste0(path,i))
  object <- copyNumbersSegmented
  template <- objectsampletotemplate(object)
  sqmodel <- squaremodel(template, prows = 100, ptop = 3, pbottom = 1, penalty = 0.6, penploidy = 0.6)

  parameter <- sqmodel$minimadf

  para1 <- parameter %>% arrange(error)
  ploidy <- para1$ploidy[1]
  cellularity <- para1$cellularity[1]

  sc_ACE_para[i,1] <- i
  sc_ACE_para[i,2] <- ploidy
  sc_ACE_para[i,3] <- cellularity
  
  m<-ACEcall(copyNumbersSegmented, QDNAseqobjectsample = TRUE, cellularity = cellularity,
             ploidy = ploidy, standard =1, plot = FALSE, title ="",
             qcutoff = -3, subcutoff = 0.2, trncname = FALSE,
             cap = 12, bottom = 0, chrsubset = 1:22, onlyautosomes = TRUE,
             sgc = c())
  
  b <- nrow(m)
  c <- copyNumbersSegmented@featureData@data[1:b,1:5]
  c[,4] <- m[,3]
  c[,5]<- i
  c <- na.omit(c)
  c[,4] <- round(c[,4])
  c <- c[,c(5,1:4)]
  colnames(c) <- c("sample","chromosome","start","end","segVal" )
  fwrite(c,file = paste0(ACE_path,i,".txt"))
  
}


fwrite(sc_ACE_para,"./pcawg_low/result/pcawg_ACE_para.txt")



filename <- list.files("~/project/scCNSignature/ACE_NEW/pcawg_low/ACE_result/")
ACEpath <- "~/project/scCNSignature/ACE_NEW/pcawg_low/ACE_result/"
P_path <- "~/project/scCNSignature/ACE_NEW/pcawg_low/Processed_result/"



for (i in filename) {
  
  data1 <- fread(paste0(ACEpath,i),data.table = F)
  data <- data1
  

  smooth_data <- smooth(data$segVal,"3RS3R")
  
  data$segVal <- smooth_data
  
  m <- 1
  finall_data <- data[1,]
  for(j in 1:c(nrow(data)-1)){
    if(data[j,5]==data[j+1,5]&data[j,2]==data[j+1,2]){
      finall_data[m,4]<-data[j+1,4]
    }
    else{
      m<-m+1
      finall_data[m,]<-data[j+1,]
    }
  }
  
  finall_data<-na.omit(finall_data)
  
  fwrite(finall_data,file = paste0(P_path,i))
}





remove(list = ls())
filename <- list.files("~/project/scCNSignature/ACE_NEW/pcawg_low/Processed_result/")
path <- "~/project/scCNSignature/ACE_NEW/pcawg_low/Processed_result/"
data <- fread("~/project/scCNSignature/ACE_NEW/pcawg_low/Processed_result/mini-001642979e361dc2c60f0128740b258b.rds.txt")
all_data <- data[1,]
all_data[1,] <- NA


for (i in filename) {
  
  data <- fread(paste0(path,i))
  all_data <- rbind(all_data,data)
  
}

scHCC_absCNA <- na.omit(all_data)
a <- as.data.frame(table(scHCC_absCNA$sample))
fwrite(scHCC_absCNA,file = "~/project/scCNSignature/ACE_NEW/pcawg_low/result/pcawg_ACE_smooth_segment.txt")


