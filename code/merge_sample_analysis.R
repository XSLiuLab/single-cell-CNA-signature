

setwd("~/project/scCNSignature/")

remove(list = ls())

library(sigminer)
library(data.table)
library(dplyr)


sc_sigs <- readRDS("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")

sc_tally <- readRDS("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/sc_tumor_tally_ploidy.rds")
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

result_matrix2 <- result_matrix1[result_matrix1[,3] < 0.95 ,]


clin <- fread("./CNS_result/scHCC_infor.csv")
clin1 <- clin %>% filter(sample_ID %in% rownames(result_matrix2))

row_annotation <- data.frame(clin1$patient1)
rownames(row_annotation) <- clin1$sample_ID
colnames(row_annotation) <- "patient"
row_annotation$Sample_id <- rownames(row_annotation)


library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(RColorBrewer)

colors <- colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100) 
values <- seq(0, 1.05, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)


Heatmap(result_matrix2, name = "Cosine",
        col = col_fun,cluster_rows = T,split = c(row_annotation$patient),
        cluster_columns = F,show_row_names = F)




test_data <- data.frame(result_matrix2)
test_data$sample_ID <- rownames(test_data)
test_data2 <- left_join(test_data,clin1)
test_data3 <- test_data2 %>% dplyr::select(paste0("scSignature",1:7),patient1)



test_result <- summarise(group_by(test_data3, patient1), mean(scSignature1),mean(scSignature2),
                         mean(scSignature3),mean(scSignature4),mean(scSignature5),mean(scSignature6),mean(scSignature7)) 

test_result$patient1 <- paste0("sc",test_result$patient1)
colnames(test_result) <- c("sample",paste0("scSig",1:7))

test_result1 <- data.frame(test_result[,-1])

rownames(test_result1) <- test_result$sample
test_result2 <- apply(test_result1, 2, function(x){round(x,1)})


colors <- colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)
values <- seq(0, 1.05, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)



Heatmap(test_result2, name = "Cosine",col = col_fun, rect_gp = gpar(type = "none"), border = TRUE,
        column_title = "scSignature Cosine of scSample",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.circle(x = x, y = y, r = abs(test_result2[i, j])/2.5 * min(unit.c(width, height)), 
                      gp = gpar(fill = col_fun(test_result2[i, j]), col = NA))
          grid.rect(x = x, y = y, width = width, height = height,
                    gp = gpar(col = "grey", fill = NA, lwd = 0.5))
        },
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE, row_names_side = "left")





sc_sigs <- readRDS("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")


sc_tally <- readRDS("./ACE_NEW/merge_data/sc_merge_tally_ploidy.rds")
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

rownames(result_matrix1) <- substr(rownames(result_matrix1),1,3)


colors <- colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)
values <- seq(0, 1.05, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)

Heatmap(result_matrix1, name = "Cosine",col = col_fun, rect_gp = gpar(type = "none"), border = TRUE,
        column_title = "scSignature Consie of merge sample",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.circle(x = x, y = y, r = abs(result_matrix1[i, j])/3 * min(unit.c(width, height)), 
                      gp = gpar(fill = col_fun(result_matrix1[i, j]), col = NA))
          grid.rect(x = x, y = y, width = width, height = height,
                    gp = gpar(col = "grey", fill = NA, lwd = 0.5))
        },
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE, row_names_side = "left")



merge_feture <- readRDS("./ACE_NEW/merge_data/sc_merge_tally_ploidy.rds")
data1 <- merge_feture$nmf_matrix
data1 <- t(data1)

sc_sig <- readRDS("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")
sc_sig_data <- sc_sig$Signature

result1 <- sig_fit(data1,sc_sig)
show_sig_fit(result1)




library(ggplot2)

compute_ratio <- function(column) {
  ratio <- column / sum(column)
  return(ratio)
}


result_matrix1 <- apply(result1, 2, compute_ratio)

show_sig_fit(result_matrix1) + ylim(0,0.8)
show_sig_fit(result_matrix2) + ylim(0,0.8)

draw_data <- t(result_matrix1)

pheatmap::pheatmap(draw_data, cluster_cols = F,cluster_rows = F, show_rownames = T,display_numbers = T,main = "scSignature Proportion of merge sample")



colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(50) 

colors <- colorRampPalette(c("#2166AC","#67A9CF","#FDDBC7","#EF8A62","#B2182B"))(50) 

values <- seq(0.15, 0.65, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)


rownames(draw_data) <- substr(rownames(draw_data),1,3)

Heatmap(draw_data, name = "Proportion",col = col_fun, rect_gp = gpar(type = "none"), border = TRUE,
        column_title = "scSignature Proportion of merge sample",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.circle(x = x, y = y, r = abs(draw_data[i, j])/2.5 * min(unit.c(width, height)), 
                      gp = gpar(fill = col_fun(draw_data[i, j]), col = NA))
          grid.rect(x = x, y = y, width = width, height = height,
                    gp = gpar(col = "grey", fill = NA, lwd = 0.5))
        },
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE, row_names_side = "left")




library(ggsci)

barplot_data <- data.frame(draw_data["P04",])
colnames(barplot_data)[1] <- "P04"
barplot_data$scSigure <- paste0("scSig",1:7)

ggplot(barplot_data, aes(x = scSigure, y = P04,fill = scSigure))+
  geom_bar(stat = "identity") +
  ggprism::theme_prism(border = T)+
  scale_fill_manual(values = pal_npg("nrc", alpha = 1)(7))+
  geom_text(aes(label = round(P04, 3)), vjust = -0.5) +
  xlab("Signature") +
  ylab("scSig proportion of bulk tissue sample P04")+
  ylim(0,1)





for (i in 1:10) {

ifelse (i < 10,sample <- paste0("P0",i),sample <- paste0("P",i)) 

  
barplot_data <- data.frame(draw_data[sample,])
colnames(barplot_data)[1] <- sample
barplot_data$scSigure <- paste0("scSig",1:7)

p <- ggplot(barplot_data, aes(x = scSigure, y = get(sample),fill = scSigure))+
  geom_bar(stat = "identity") +
  ggprism::theme_prism(border = T)+
  scale_fill_manual(values = pal_npg("nrc", alpha = 1)(7))+
  geom_text(aes(label = round(get(sample), 3)), vjust = -0.5) +
  xlab("Signature") +
  ylab(paste0("scSig proportion of bulk tissue sample ",sample))+
  ylim(0,1)
  

pdf(file.path("./ACE_NEW/figure/P04_supp/other_supp_bulk_sc/bulk/", paste0(sample, ".pdf")),
    width = 8, height = 6)
print(p, newpage = FALSE)
dev.off()

}



sc_sig <- readRDS("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")
sc_sig_data <- sc_sig$Exposure.norm
sc_sig_data1 <- data.frame(t(sc_sig_data))

sc_sig_data1$sample_ID <- rownames(sc_sig_data1)

scHCC_infor <- fread("./ACE_NEW/result/scHCC_infor.txt")
scHCC_infor1 <- scHCC_infor %>% dplyr::select(sample_ID,patient1)

draw_data <- inner_join(sc_sig_data1,scHCC_infor1)
draw_data <- draw_data[draw_data[,3] < 0.95,]

row_annotation <- data.frame(draw_data$patient1)
rownames(row_annotation) <- draw_data$sample_ID
colnames(row_annotation) <- "patient"
row_annotation$Sample_id <- rownames(row_annotation)

draw_data2 <- draw_data[,-c(8:9)]

library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(RColorBrewer)

colors <- colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100) 

colors <- colorRampPalette(c("#67A9CF","#FDDBC7","#EF8A62","#B2182B"))(100) 
values <- seq(0, 0.8, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)


Heatmap(as.matrix(draw_data2), name = "Proportion",
        col = col_fun,cluster_rows = T,split = c(row_annotation$patient),row_gap = unit(0.2, "cm"),
        cluster_columns = F,show_row_names = F)

scatte_data <- draw_data2
scatte_data$sample<- row_annotation$patient
scatte_data$id <- row_annotation$Sample_id
scatte_data1 <- scatte_data %>% dplyr::filter(sample == "P04") %>% 
  dplyr::select(-sample) %>% data.frame()

scatte_data2 <- scatte_data1[,-8]
colnames(scatte_data2)[1:7] <- paste0("scSig",1:7)
rownames(scatte_data2) <- scatte_data1$id
scatte_data2$sample <- rownames(scatte_data2)
scatte_data3 <- melt(scatte_data2)




ggplot(scatte_data3, aes(x = variable, y = value,color = variable))+
  geom_line(aes(group = sample), color = 'gray', lwd = 0.01)+
  geom_point(size = 2) +
  ggprism::theme_prism(border = T)+
  scale_color_manual(values = pal_npg("nrc", alpha = 0.5)(7))+
  xlab("Signature") +
  ylab("scSig proportion of single cell sample 4")+
  ylim(0,1)



for (i in 1:10) {

  ifelse (i < 10,infor <- paste0("P0",i),infor <- paste0("P",i)) 
  
  
  
  scatte_data <- draw_data2
  scatte_data$sample<- row_annotation$patient
  scatte_data$id <- row_annotation$Sample_id
  scatte_data1 <- scatte_data %>% dplyr::filter(sample == infor) %>% 
    dplyr::select(-sample) %>% data.frame()
  
  scatte_data2 <- scatte_data1[,-8]
  colnames(scatte_data2)[1:7] <- paste0("scSig",1:7)
  rownames(scatte_data2) <- scatte_data1$id
  scatte_data2$sample <- rownames(scatte_data2)
  scatte_data3 <- melt(scatte_data2)
  
  
p <- ggplot(scatte_data3, aes(x = variable, y = value,color = variable))+
    geom_line(aes(group = sample), color = 'gray', lwd = 0.01)+
    geom_point(size = 2) +
    ggprism::theme_prism(border = T)+
    scale_color_manual(values = pal_npg("nrc", alpha = 1)(7))+
    xlab("Signature") +
    ylab(paste0("scSig proportion of single cell sample",sample))+
    ylim(0,1)
  
  pdf(file.path("./ACE_NEW/figure/P04_supp/other_supp_bulk_sc/sc/", paste0("sc",infor, ".pdf")),
      width = 8, height = 6)
  print(p, newpage = FALSE)
  dev.off()
  
}



test_data <- draw_data[,c(1:7,9)]

test_result <- summarise(group_by(test_data, patient1), mean(scSignature1),mean(scSignature2),
          mean(scSignature3),mean(scSignature4),mean(scSignature5),mean(scSignature6),mean(scSignature7)) 

test_result$patient1 <- paste0("sc",test_result$patient1)
colnames(test_result) <- c("sample",paste0("scSig",1:7))

test_result1 <- data.frame(test_result[,-1])

rownames(test_result1) <- test_result$sample
test_result2 <- apply(test_result1, 2, function(x){round(x,1)})
  

colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(50) 
colors <- colorRampPalette(c("#2166AC","#67A9CF","#FDDBC7","#EF8A62","#B2182B"))(50)
values <- seq(0.15, 0.65, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)



Heatmap(test_result2, name = "Proportion",col = col_fun, rect_gp = gpar(type = "none"), border = TRUE,
        column_title = "scSignature Proportion of scSample",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.circle(x = x, y = y, r = abs(test_result2[i, j])/2.5 * min(unit.c(width, height)), 
                      gp = gpar(fill = col_fun(test_result2[i, j]), col = NA))
          grid.rect(x = x, y = y, width = width, height = height,
                    gp = gpar(col = "grey", fill = NA, lwd = 0.5))
        },
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE, row_names_side = "left")







setwd("~/project/scCNSignature/")
remove(list = ls())
library(sigminer)

cna <- fread("./ACE_NEW/merge_data/sc_merge_ACE_smooth_segment.txt")

cna$sample <- substr(cna$sample,1,3)


cn <- read_copynumber(cna, seg_cols = c("chromosome", "start", "end", "segVal"), genome_build = "hg38", complement = F, verbose = TRUE) 
show_cn_profile(cn, nrow = 2, ncol = 1,show_title = T)



sel_samps <- c("P01","P02","P03","P04","P05",
               "P06","P07","P08","P09","P10")

cn_max <- c(4,7,3,4,4,7,3,5,7,7)
m <- 0

plist <- list()
for (s in sel_samps) {
  use_cn <- cna %>% filter(sample %in% c(s))
  cn <- read_copynumber(cna, seg_cols = c("chromosome", "start", "end", "segVal"), genome_build = "hg38", complement = F, verbose = TRUE) 
  

  m <- m+1
  max <- cn_max[m]
  
  p <- show_cn_profile(cn,
                       samples = s,
                       show_title = TRUE,
                       return_plotlist = TRUE,
                       ylim = c(0,max))
  
  plist[[s]] <- p[[1]] + ggplot2::labs(x = NULL)
}


p <- cowplot::plot_grid(plotlist = plist, ncol = 1)

ggplot2::ggsave(file.path("./ACE_NEW/merge_data/sample","samples.pdf"),
                plot = p, width = 12, height = 42)








setwd("~/project/scCNSignature/")
remove(list = ls())


library(data.table)
library(dplyr)
library(sigminer)
library(CNpare)


sample_data <- fread("./merge_sample/1Mb/scHCC_merge_smooth_absCNA_segment.txt",data.table = F)
cell_data <- fread("./CNS_result/scHCC_smooth_absCNA_segment.txt",data.table = F)
cell_infor <- fread("./CNS_result/scHCC_infor.csv",data.table = F)


P01_sample <- sample_data %>% filter(sample == "P01") %>% select(chromosome,start,end,segVal,sample)

P01_infor <- cell_infor %>% filter(patient1 == "P01")
P01_cell <- cell_data %>% filter(sample==P01_infor$sample_ID[2]) %>% select(chromosome,start,end,segVal,sample)



posBins <- lapply(seq_len(22),function(chr)
  getBinsStartsEnds(window=1000000, chr, lengthChr[chr]))

P01_cell$chromosome <- as.character(P01_cell$chromosome)
cells_bin <- getCNbins(posBins=posBins,
                       data=P01_cell, samples=P01_cell$sample[1])

P01_sample$chromosome <- as.character(P01_sample$chromosome)
sample_bin <- getCNbins(posBins=posBins,
                        data=P01_sample, samples="P01")


a <- CNpare::getSimilarities(sample_bin,cells_bin,method="cosine")





scdata <- cell_data %>% select(chromosome,start,end,segVal,sample)
scdata$chromosome <- as.character(scdata$chromosome)
bulkdata <- sample_data %>% select(chromosome,start,end,segVal,sample)
bulkdata$chromosome <- as.character(bulkdata$chromosome)

bin.size=1000
allchr=c(1:22) 
lengthChr=lengthChr
posBins <- lapply(allchr,function(chr) 
  getBinsStartsEnds(window=bin.size*1000, chr, lengthChr[chr]))



var_cell <- cell_infor %>% filter(CNA > 40)
var_data <- cell_data %>% filter(sample %in% var_cell$sample_ID) %>% select(chromosome,start,end,segVal,sample)
var_data$chromosome <- as.character(var_data$chromosome)


sc_var <- getCNbins(posBins=posBins, data=var_data, samples=unique(var_data$sample))
bulk_cn <- getCNbins(posBins=posBins, data=bulkdata, samples=unique(bulkdata$sample))
measures2 <- getSimilarities(dat1=bulk_cn, dat2=sc_var, method="all")


data1 <- measures2 %>% select(fileid,id,cos_sim,manhattan,euclidean)
colnames(data1) <- c("cosin_ID","sample_ID","cos_sim","manhattan","euclidean")

data2 <- left_join(data1,cell_infor)
data3 <- data2 %>% filter(patient1 == "P01")

library(ggplot2)
library(ggprism)
ggplot(data3, aes(x = cosin_ID, y = cos_sim)) +
  geom_boxplot()+theme_prism()


sel_samps <- c("P01","P02","P03","P04","P05",
               "P06","P07","P08","P09","P10")

plist <- list()
for (s in sel_samps) {
  
  data3 <- data2 %>% filter(patient1 == s)
  
  p <- ggplot(data3, aes(x = cosin_ID, y = cos_sim))+
    geom_boxplot()+theme_prism()+labs(title = paste0(s,"_sc"))
  
  plist[[s]] <- p + ggplot2::labs(x = NULL)
}



p <- cowplot::plot_grid(plotlist = plist, ncol = 1)

ggplot2::ggsave(file.path("./ACE_NEW/merge_data/sample","boxplot_cos_sim.pdf"),
                plot = p, width = 7, height = 38)




sel_samps <- c("P01","P02","P03","P04","P05",
               "P06","P07","P08","P09","P10")

data3 <- data2 %>% filter(patient1 == "P01")
a <- data3 %>% group_by(cosin_ID) %>% summarise(cosin_mean = mean(manhattan))

all_cosin_value <- a[,1]
for (sample_name in sel_samps) {
  data3 <- data2 %>% filter(patient1 == sample_name)
  a <- data3 %>% group_by(cosin_ID) %>% summarise(cosin_mean = mean(manhattan))
  
  colnames(a)[2] <- sample_name
  all_cosin_value <- left_join(all_cosin_value,a)
  
}


all_cosin_value1 <- all_cosin_value[,-1]
rownames(all_cosin_value1) <- sel_samps
colnames(all_cosin_value1) <- paste0(sel_samps,"_sc")



Norm <- apply(all_cosin_value1, 2, function(x){(x-mean(x))/sd(x)})

Norm1 <- apply(Norm, 1, function(x){(x-min(x))/(max(x)-min(x))})


library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "Gainsboro", "red"))
Heatmap(Norm1,cluster_columns = F,row_names_side = "left",show_row_names = T,cluster_rows = F,show_column_names = T,name = " ",col = col_fun)



library(pheatmap)
pheatmap::pheatmap(Norm1, cluster_cols = F, cluster_rows = F, display_numbers = TRUE,main = "Normlized manhattan distance")


