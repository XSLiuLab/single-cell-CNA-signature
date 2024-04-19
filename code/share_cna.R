

remove(list = ls())
setwd("~/project/scCNSignature/")


cna_data <- fread("./ACE_NEW/result/sc_tumor_ACE_smooth_segment.txt",data.table = F)
cna_infor <- fread("./ACE_NEW/result/scHCC_infor.txt",data.table = F)


patient_name <- names(table(cna_infor$patient1))
all_shared_cna <- list()

for (p1 in patient_name) {
  
  use_infor <- cna_infor %>% dplyr::filter(patient1 == p1)
  use_cna <- cna_data %>% dplyr::filter(sample %in% use_infor$sample_ID)
  all_sample <- names(table(use_cna$sample))
  
  shared_cna  <- as.data.frame(matrix(0,ncol = 2,nrow = length(all_sample)))
  colnames(shared_cna) <- c("sample", "shared_cna")
  shared_cna$sample <- all_sample
  rownames(shared_cna) <- all_sample
  
  for (s1 in all_sample) {
    for (s2 in all_sample) {
     
      if (s1 != s2) {
        
        data1 <- use_cna %>% dplyr::filter(sample == s1)
        data2 <- use_cna %>% dplyr::filter(sample == s2)
        shared_point <- sum(paste0("chr_",data1$chromosome,"_",data1$start) %in% paste0("chr_",data2$chromosome,"_",data2$start))
        shared_cna[s1,2] <- shared_cna[s1,2] + shared_point 
        
      }
      
    }
  }
  
  all_shared_cna[[p1]] <- shared_cna
  
}

saveRDS(all_shared_cna,"./result/sc_all_shared_cna.rds")



library(dplyr)
library(data.table)


all_shared_cna <- readRDS("./result/sc_all_shared_cna.rds")
P09 <- all_shared_cna$P09
len <- nrow(P09)
P09$average <-round( P09$shared_cna/len)
all_sc_sample <- fread("./ACE_NEW/scData/sc_sig_sample_new.txt",data.table = F)
use_data <- all_sc_sample %>% dplyr::filter(patient1 == "P09") %>% dplyr::select(1,8:14)
colnames(use_data)[1] <- "sample"
use_data2 <- dplyr::inner_join(P09,use_data) %>% dplyr::arrange(shared_cna)
colnames(use_data2)[4:10] <- paste0("scSig",1:7)


library(tidyr)

draw_data <- tidyr::pivot_longer(
  use_data2,
  cols = -c(average, sample, shared_cna),  
  names_to = "Sig", 
  values_to = "value" 
)


draw_data$shared_cna <- round(draw_data$shared_cna/len)
draw_data$sample <- factor(draw_data$sample,levels = use_data2$sample)
draw_data$value <- draw_data$value
draw_data$Sig <- factor(draw_data$Sig,levels = paste0("scSig",1:7))


library(ggplot2)
library(ggsci)
library(scales)
library(RColorBrewer)



ggplot(draw_data, aes(x = sample, y = value, group = Sig, fill = Sig)) + 
  geom_bar(stat="identity",position="stack", width=0.7,size=0.25)+
  scale_fill_manual(values =  pal_npg("nrc", alpha = 1)(10)[c(1:6,9)])+
  ggprism::theme_prism(border = T)+
  scale_x_discrete(breaks = c("HRR227487", "HRR227472"),labels = c("23", "43"))+
  labs(x = "Shared CNA", y = "The proportion of scSig")





remove(list = ls())
setwd("~/project/scCNSignature/")
cna_data <- fread("./ACE_NEW/result/sc_tumor_ACE_smooth_segment.txt",data.table = F)
cna_infor <- fread("./ACE_NEW/result/scHCC_infor.txt",data.table = F)

patient_name <- names(table(cna_infor$patient1))
shared_break <- list()

for (p1 in patient_name) {
  
  use_infor <- cna_infor %>% dplyr::filter(patient1 == p1)
  use_cna <- cna_data %>% dplyr::filter(sample %in% use_infor$sample_ID)
  all_sample <- names(table(use_cna$sample))
  
  use_cna$break_point <- paste0("chr_",use_cna$chromosome,"_",use_cna$start)
  sample_point <- as.data.frame(table(use_cna$break_point))
  
  
  shared_cna  <- as.data.frame(matrix(0,ncol = 2,nrow = length(all_sample)))
  colnames(shared_cna) <- c("sample", "shared_cna")
  shared_cna$sample <- all_sample
  rownames(shared_cna) <- all_sample
  
  for (s1 in all_sample) {

    data1 <- use_cna %>% dplyr::filter(sample == s1)
    shared_point <- sample_point %>% dplyr::filter(Var1 %in% data1$break_point)
    
    shared_cna[s1,2] <- round(sum(shared_point$Freq)/nrow(shared_point)/nrow(shared_cna),3)
  }
  
  shared_break[[p1]] <- shared_cna
  
}


saveRDS(shared_break,"./result/shared_break.rds")



shared_break <- readRDS("./result/shared_break.rds")
all_shared_cna <- readRDS("./result/sc_all_shared_cna.rds")

all_shared_data <- list()

for (p1 in names(shared_break)) {

  P09 <- all_shared_cna[[p1]]
  len <- nrow(P09)
  P09$average <-round( P09$shared_cna/len)
  all_sc_sample <- fread("./ACE_NEW/scData/sc_sig_sample_new.txt",data.table = F)
  use_data <- all_sc_sample %>% dplyr::filter(patient1 == p1) %>% dplyr::select(1,8:14)
  colnames(use_data)[1] <- "sample"
  use_data2 <- dplyr::inner_join(P09,use_data) %>% dplyr::arrange(shared_cna)
  colnames(use_data2)[4:10] <- paste0("scSig",1:7)

  shared_point <- shared_break[[p1]]
  colnames(shared_point)[2] <- "shared_break"

  shared_data <- dplyr::left_join(shared_point, use_data2) %>% dplyr::arrange(-shared_break)

  all_shared_data[[p1]] <- shared_data
}


saveRDS(all_shared_data,"./result/all_shared_data.rds")


library(tidyr)

all_shared_data <- readRDS("./result/all_shared_data.rds")

file_name <- names(all_shared_data)

for (j in file_name) {
  
use_data1 <- all_shared_data[[j]]
use_data1$sample <- factor(use_data1$sample,levels = use_data1$sample)

draw_data <- tidyr::pivot_longer(
  use_data1,
  cols = -c(average, sample, shared_cna, shared_break),  
  names_to = "Sig",  
  values_to = "value" 
)



draw_data$sample <- factor(draw_data$sample,levels = use_data1$sample)
draw_data$value <- draw_data$value
draw_data$Sig <- factor(draw_data$Sig,levels = paste0("scSig",1:7))

library(ggplot2)
library(ggsci)
library(scales)
library(RColorBrewer)

p <- ggplot(draw_data, aes(x = sample, y = value, group = Sig, fill = Sig)) + 
  geom_bar(stat="identity",position="stack", width=0.7,size=0.25, show.legend = FALSE)+
  scale_fill_manual(values =  pal_npg("nrc", alpha = 1)(10)[c(1:6,9)])+
  ggprism::theme_prism(border = T)+
  scale_x_discrete(breaks = c(""),labels = c(""))+
  labs(x = "Ranking based on the sum of shared CNA", y = "Proportion of scSig")

p

p2 <- ggplot(data = use_data1, aes(x = sample, y = shared_break, group = 1)) +
  geom_line(color = "red") +
  labs(y = "Average proportion of shared CNA") +
  ylim(0,1)+
  ggprism::theme_prism(border = T)+
  scale_x_discrete(breaks = c(""),labels = c(""))

p2

library(cowplot)
combined_plot <- plot_grid(p + theme(axis.title.x = element_blank()),
                           p2 + theme(axis.title.x = element_text()),
                           ncol = 1, align = "v",rel_heights = c(2, 1))
combined_plot

pdf(file.path(paste0("./ACE_NEW/figure/shared_point/",j,".pdf")),
    width = 8, height = 6)
print(combined_plot, newpage = FALSE)
dev.off()
}
