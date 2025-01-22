
setwd("~/project/scCNSignature/")
remove(list = ls())

library(sigminer)
library(data.table)
library(dplyr)


sc_sigs <- readRDS("./data/Signature/signature/sc_ploidy_sigs_signature.rds")



sc_tally <- readRDS("data/Signature/NMF Matrix/scWGS_tally_ploidy.rds")
sc_feture <- sc_tally$nmf_matrix
sc_feture1 <- t(sc_feture)
sc_feture2 <- sc_feture1

data2_2 <- sc_sigs$Signature

result_matrix <- matrix(0, nrow = ncol(sc_feture2), ncol = 7)
colnames(result_matrix) <- colnames(sc_sigs$Signature)
rownames(result_matrix) <- colnames(sc_feture1)



for (j in 1:ncol(sc_feture2)) {
  vector1 <- sc_feture2[,j]
  
  for (i in 1:7) {
    vector2 <- data2_2[,i]
    
    result_matrix[j,i] <- cosine(vector1, vector2)
    
  }
  
}

result_matrix1 <- na.omit(result_matrix)
result_matrix1 <- result_matrix1


result_matrix2 <- result_matrix1


clin <- fread("./ACE/result/scHCC_infor.txt")
clin1 <- clin %>% filter(sample_ID %in% rownames(result_matrix2))


sc_sample <- as.data.frame(result_matrix2)
sc_sample$sample_ID <- rownames(sc_sample)
sc_sig_sample <- left_join(clin1,sc_sample)
fwrite(sc_sig_sample,"./ACE/scData/sc_sig_sample.txt")
sc_sig_sample <- fread("./ACE/scData/sc_sig_sample.txt",data.table = F)

P06 <- sc_sig_sample %>% filter(patient1 == "P06")

cna <- fread("./ACE/result/sc_ACE_smooth_segment.txt")

use_cn <- cna %>% filter(sample %in% c("HRR226897","HRR226897"))
cn <- read_copynumber(use_cn, seg_cols = c("chromosome", "start", "end", "segVal"), genome_build = "hg38", complement = F, verbose = TRUE) 
show_cn_profile(cn, nrow = 2, ncol = 1,show_title = T,ylim = c(0,max(use_cn$segVal)))


use_cn <- cna %>% filter(sample %in% c("HRR226874","HRR226874"))
cn <- read_copynumber(use_cn, seg_cols = c("chromosome", "start", "end", "segVal"), genome_build = "hg38", complement = F, verbose = TRUE) 
show_cn_profile(cn, nrow = 2, ncol = 1,show_title = T,ylim = c(0,5))

use_cn <- cna %>% filter(sample %in% c("HRR226828","HRR226828"))
cn <- read_copynumber(use_cn, seg_cols = c("chromosome", "start", "end", "segVal"), genome_build = "hg38", complement = F, verbose = TRUE) 
show_cn_profile(cn, nrow = 2, ncol = 1,show_title = T,ylim = c(0,max(use_cn$segVal)))





### envolution

all_sc_sample <- fread("./ACE/scData/sc_sig_sample.txt",data.table = F)
all_sc_sample <- fread("./ACE/scData/sc_sig_sample_new.txt",data.table = F)
sc_sample_infor <- all_sc_sample %>% dplyr::select("sample_ID","patient1")


envol_data <- all_sc_sample %>% dplyr::select(paste0("scSignature",1:7)) %>% as.data.frame()
rownames(envol_data) <- all_sc_sample$sample_ID

main_sig <- data.frame(apply(envol_data,1,FUN = function(x){which(x == max(x))}))
main_sig <- as.data.frame(t(main_sig[1,]))

colnames(main_sig) <- "main_sig"
main_sig$sample_ID <- rownames(main_sig)

envol_data1 <- envol_data
envol_data1$sample_ID <- rownames(envol_data1)
envol_data2 <- left_join(envol_data1,main_sig) 
envol_data3 <- left_join(envol_data2,sc_sample_infor)

envol_data3$main_sig <- as.factor(envol_data3$main_sig)

test <- envol_data3 %>% filter(patient1 == "P01")
table(test$main_sig)




all_result <- data.frame(matrix(NA,nrow = 10,ncol = 7))
colnames(all_result) <- c(paste0("scSig",1:7))
rownames(all_result) <- c(paste0("P0",1:9),"P10")


for (sample in c(paste0("P0",1:9),"P10")) {
  
  test <- envol_data3 %>% filter(patient1 == sample)
  all_result[sample,] <- as.numeric(round(table(test$main_sig)/nrow(test),3))
    

}




draw_figure <- all_result
draw_figure$sample <- rownames(draw_figure)

draw_figure1 <- reshape2::melt(draw_figure)


library(ggplot2)
library(ggsci)
library(scales)
library(RColorBrewer)
fwrite(all_result,"~/BioXCG/geom_bar_data.txt")

mypal <- pal_npg("nrc", alpha = 0.7)(10)

mypal
show_col(mypal,ncol = 5)

ggplot(data=draw_figure1,aes(sample,value,fill=variable))+ 
  geom_bar(stat="identity",position="stack", width=0.7,size=0.25)+
  scale_fill_manual(values = pal_npg("nrc", alpha = 1)(10)[c(1:6,9)])+
  ggprism::theme_prism(border = T)+
  labs(y="The proportion of single cells")
  

a <- draw_figure1 %>% dplyr::filter(variable == "scSig3")
a$sur <- c(21,52,45,30,27,NA,28,28,16,28)
b <- na.omit(a)
cor.test(b$value,b$sur)


test <- envol_data3 %>% filter(patient1 == "P10")

tree_test <- test %>% select(paste0("scSignature",1:7)) %>% as.data.frame()


library(ape)
library(phangorn)
library(ggtree)
library(ggstance)

dist_matrix <- dist(tree_test)
tree <- nj(dist_matrix)

table(test$main_sig)
node_root <- which(tree_test$scSignature3 == max(tree_test$scSignature3))
tree <- root(tree, outgroup = node_root[1])


colors_df <- as.data.frame(test[,8:10])
colors_df <- colors_df %>% dplyr::mutate(color = case_when(main_sig == 1 ~ pal_npg("nrc", alpha = 1)(7)[1],
                                                           main_sig == 2 ~ pal_npg("nrc", alpha = 1)(7)[2],
                                                           main_sig == 3 ~ pal_npg("nrc", alpha = 1)(7)[3],
                                                           main_sig == 4 ~ pal_npg("nrc", alpha = 1)(7)[4],
                                                           main_sig == 5 ~ pal_npg("nrc", alpha = 1)(7)[5],
                                                           main_sig == 6 ~ pal_npg("nrc", alpha = 1)(7)[6],
                                                           main_sig == 7 ~ pal_npg("nrc", alpha = 1)(10)[9]))


colors_df <- colors_df %>% dplyr::mutate(shape = case_when(main_sig == 1 ~ 3,
                                                           main_sig == 2 ~ 4,
                                                           main_sig == 3 ~ 8,
                                                           main_sig == 4 ~ 15,
                                                           main_sig == 5 ~ 16,
                                                           main_sig == 6 ~ 17,
                                                           main_sig == 7 ~ 18))


p <- ggtree(tree) + geom_tippoint(color=colors_df$color, shape= colors_df$shape, size=3)

p

a <- na.omit(p$data) %>% dplyr::arrange(y)
draw_1 <- tree_test[a$label,]
draw_1$cell <- rownames(draw_1)
draw_2 <- reshape2::melt(draw_1)
draw_2$cell <- factor(draw_2$cell,levels = a$label)


library(ggstance)
facet_plot(p, panel="scSig of single cell", data=draw_2, geom=geom_barh, mapping=aes(value,fill=variable), 
           stat='identity',width=0.7,size=0.25) + 
  scale_fill_manual(values =pal_npg("nrc", alpha = 1)(10)[c(1:6,9)])+
  theme(legend.position = "right")+
  theme_classic()



