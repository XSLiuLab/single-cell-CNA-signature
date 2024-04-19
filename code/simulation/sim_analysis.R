

setwd("~/project/scCNSignature/ACE_NEW/nature_sim/")
remove(list = ls())
library(dplyr)
library(data.table)



## positive simulation analysis

data1 <- readRDS("./output/Out_2_Segmentation_tables_20each_N240.rds")

all_sample <- data1[[1]]
all_sample <- all_sample[1,]
all_sample[1,] <- NA


for (i in c(1:2400)) {
  
  sample1 <- data1[[i]]
  all_sample <- rbind(all_sample,sample1)
  
}

all_sample1 <- na.omit(all_sample)

fwrite(all_sample1,"./result/Nature_simulation_sample_raw.txt")




all_sample2 <- all_sample1 %>% dplyr::select(-signature)
a <- as.data.frame(table(all_sample2$sample))
all_sample3 <- all_sample2
all_sample3$segVal <- round(all_sample3$segVal)

result1 <- divide_feature(all_sample3)

saveRDS(result1,"./result/Nature_simulation_nmf.rds")



remove(list = ls())
library(dplyr)
library(data.table)

nature_sim <- readRDS("~/project/scCNSignature/ACE_NEW/nature_sim/result/Nature_simulation_nmf.rds")
nature_nmf <- t(nature_sim)


sc_sig <- readRDS("~/project/scCNSignature/ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")
sc_sig_data <- sc_sig$Signature
result1 <- sig_fit(nature_nmf,sc_sig)
norm_result <- apply(result1, 2, function(x){x/sum(x)})

new_fit_result <- list()
new_fit_result$fit_result <- result1
new_fit_result$norma_fit_result <- norm_result
saveRDS(new_fit_result,"./result/NatureSim_fit_result.rds")


new_fit_result <- readRDS("~/project/scCNSignature/ACE_NEW/nature_sim/result/NatureSim_fit_result.rds")
nature_fit <- data.frame(t(new_fit_result$norma_fit_result))
nature_fit$Samples <- rownames(nature_fit)

sample_infor <- fread("~/project/scCNSignature/ACE_NEW/nature_sim/output/Out_2_Overview_simulation_20each_N240.txt")
sample_infor <- sample_infor %>% dplyr::select(Samples,Signature)

table(sample_infor$Signature)
draw_data <- left_join(sample_infor,nature_fit)
draw_data1 <- draw_data %>% filter(Signature %in% c("LST","WGDearly_Chr_LST","WGDlate_LST_Chr")) %>%
  dplyr::mutate(Signature = case_when(Signature == "WGDearly_Chr_LST"~"WGDearly",
                                      Signature == "WGDlate_LST_Chr"~"WGDlate",
                                      Signature == "LST"~"LST"))

LST <- draw_data1 %>% dplyr::filter(Signature == "LST") %>% dplyr::arrange(-scSignature6) %>% dplyr::slice(1:20) 
WGDearly <- draw_data1 %>% dplyr::filter(Signature == "WGDearly") %>% dplyr::slice(1:20) 
WGDlate <- draw_data1 %>% dplyr::filter(Signature == "WGDlate") %>% dplyr::arrange(-scSignature5) %>% dplyr::slice(1:20) 



draw_data1 <- rbind(LST,WGDearly,WGDlate)



result_matrix1 <- as.matrix(draw_data1[,3:9])



colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
values <- seq(0, 0.7, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)


Heatmap(result_matrix1, name = "Cosine",column_title ="SimSample vs scSignature",
        col = col_fun,cluster_rows = T,split = c(draw_data1$Signature),
        cluster_columns = F,show_row_names = F)





nature_sim <- readRDS("~/project/scCNSignature/ACE_NEW/nature_sim/result/Nature_simulation_nmf.rds")
nature_nmf <- t(nature_sim)


sc_sig <- readRDS("~/project/scCNSignature/called_signature/tcga.result/data/tcga_cn_sigs_signature.rds")
sc_sig_data <- sc_sig$Signature
result1 <- sig_fit(nature_nmf,sc_sig)
norm_result <- apply(result1, 2, function(x){x/sum(x)})

new_fit_result <- list()
new_fit_result$fit_result <- result1
new_fit_result$norma_fit_result <- norm_result
saveRDS(new_fit_result,"./result/NatureSim_fit_TCGA_result.rds")


new_fit_result <- readRDS("./result//NatureSim_fit_TCGA_result.rds")
nature_fit <- data.frame(t(new_fit_result$norma_fit_result))
nature_fit$Samples <- rownames(nature_fit)

sample_infor <- fread("~/project/scCNSignature/ACE_NEW/nature_sim/output/Out_2_Overview_simulation_20each_N240.txt")
sample_infor <- sample_infor %>% dplyr::select(Samples,Signature)


draw_data <- left_join(sample_infor,nature_fit)
draw_data1 <- draw_data %>% filter(Signature %in% c("LST","WGDearly_Chr_LST","WGDlate_LST_Chr")) %>%
  dplyr::mutate(Signature = case_when(Signature == "WGDearly_Chr_LST"~"WGDearly",
                                      Signature == "WGDlate_LST_Chr"~"WGDlate",
                                      Signature == "LST"~"LST"))

LST <- draw_data1 %>% dplyr::filter(Signature == "LST") %>% dplyr::slice(1:20) 
WGDearly <- draw_data1 %>% dplyr::filter(Signature == "WGDearly") %>% dplyr::slice(1:20) 
WGDlate <- draw_data1 %>% dplyr::filter(Signature == "WGDlate") %>% dplyr::slice(1:20) 

draw_data1 <- rbind(LST,WGDearly,WGDlate)



result_matrix1 <- as.matrix(draw_data1[,3:10])


colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
values <- seq(0, 0.7, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)


Heatmap(result_matrix1, name = "Cosine",column_title ="SimSample vs TCGA",
        col = col_fun,cluster_rows = T,split = c(draw_data1$Signature),
        cluster_columns = F,show_row_names = F)






remove(list = ls())
nature_sim <- readRDS("~/project/scCNSignature/ACE_NEW/nature_sim/result/Nature_simulation_nmf.rds")
nature_nmf <- t(nature_sim)


sc_sig <- readRDS("~/project/scCNSignature/called_signature/pcawg.result_deep/data/pcawg_deep_cn_sigs_signature.rds")
sc_sig_data <- sc_sig$Signature
result1 <- sig_fit(nature_nmf,sc_sig)
norm_result <- apply(result1, 2, function(x){x/sum(x)})

new_fit_result <- list()
new_fit_result$fit_result <- result1
new_fit_result$norma_fit_result <- norm_result
saveRDS(new_fit_result,"./result/NatureSim_fit_WGS_result.rds")


new_fit_result <- readRDS("./result//NatureSim_fit_WGS_result.rds")
nature_fit <- data.frame(t(new_fit_result$norma_fit_result))
nature_fit$Samples <- rownames(nature_fit)

sample_infor <- fread("~/project/scCNSignature/ACE_NEW/nature_sim/output/Out_2_Overview_simulation_20each_N240.txt")
sample_infor <- sample_infor %>% dplyr::select(Samples,Signature)


draw_data <- left_join(sample_infor,nature_fit)
draw_data1 <- draw_data %>% filter(Signature %in% c("LST","WGDearly_Chr_LST","WGDlate_LST_Chr")) %>%
  dplyr::mutate(Signature = case_when(Signature == "WGDearly_Chr_LST"~"WGDearly",
                                      Signature == "WGDlate_LST_Chr"~"WGDlate",
                                      Signature == "LST"~"LST"))




LST <- draw_data1 %>% dplyr::filter(Signature == "LST") %>% dplyr::slice(1:20) 
WGDearly <- draw_data1 %>% dplyr::filter(Signature == "WGDearly") %>% dplyr::slice(1:20) 
WGDlate <- draw_data1 %>% dplyr::filter(Signature == "WGDlate") %>% dplyr::slice(1:20) 

draw_data1 <- rbind(LST,WGDearly,WGDlate)

result_matrix1 <- as.matrix(draw_data1[,3:9])


colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
values <- seq(0, 0.7, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)


Heatmap(result_matrix1, name = "Cosine",column_title ="SimSample vs WGS",
        col = col_fun,cluster_rows = T,split = c(draw_data1$Signature),
        cluster_columns = F,show_row_names = F)



#### WGS
remove(list = ls())
nature_sim <- readRDS("~/project/scCNSignature/ACE_NEW/nature_sim/result/Nature_simulation_nmf.rds")
nature_nmf <- t(nature_sim)


sc_sig <- readRDS("~/project/scCNSignature/ACE_NEW/call_signature/ploidy_pcawgSignature/result/pcawg_ploidy_sigs_signature.rds")
sc_sig_data <- sc_sig$Signature
result1 <- sig_fit(nature_nmf,sc_sig)
norm_result <- apply(result1, 2, function(x){x/sum(x)})

new_fit_result <- list()
new_fit_result$fit_result <- result1
new_fit_result$norma_fit_result <- norm_result
saveRDS(new_fit_result,"./result/NatureSim_fit_sWGS_result.rds")


new_fit_result <- readRDS("./result//NatureSim_fit_sWGS_result.rds")
nature_fit <- data.frame(t(new_fit_result$norma_fit_result))
nature_fit$Samples <- rownames(nature_fit)

sample_infor <- fread("~/project/scCNSignature/ACE_NEW/nature_sim/output/Out_2_Overview_simulation_20each_N240.txt")
sample_infor <- sample_infor %>% dplyr::select(Samples,Signature)


draw_data <- left_join(sample_infor,nature_fit)
draw_data1 <- draw_data %>% filter(Signature %in% c("LST","WGDearly_Chr_LST","WGDlate_LST_Chr")) %>%
  dplyr::mutate(Signature = case_when(Signature == "WGDearly_Chr_LST"~"WGDearly",
                                      Signature == "WGDlate_LST_Chr"~"WGDlate",
                                      Signature == "LST"~"LST"))



LST <- draw_data1 %>% dplyr::filter(Signature == "LST") %>% dplyr::slice(1:20) 
WGDearly <- draw_data1 %>% dplyr::filter(Signature == "WGDearly") %>% dplyr::slice(1:20) 
WGDlate <- draw_data1 %>% dplyr::filter(Signature == "WGDlate") %>% dplyr::slice(1:20) 

draw_data1 <- rbind(LST,WGDearly,WGDlate)

result_matrix1 <- as.matrix(draw_data1[,3:9])


colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
values <- seq(0, 0.7, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)


Heatmap(result_matrix1, name = "Cosine",column_title ="SimSample vs WGS",
        col = col_fun,cluster_rows = T,split = c(draw_data1$Signature),
        cluster_columns = F,show_row_names = F)










## negative simulation analysis

sc_cns <- fread("~/project/scCNSignature/CNS_result/scHCC_smooth_absCNA_segment.txt",data.table = F)
normal_sample <- sc_cns %>% dplyr::filter(sample == "HRR226278")
normal_nmf <- divide_feature(normal_sample)
normal_nmf <- t(normal_nmf)

sc_sig <- readRDS("~/project/scCNSignature/ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")
sc_sig <- readRDS("~/project/scCNSignature/ACE_NEW/call_signature/ploidy_pcawgSignature/result/pcawg_ploidy_sigs_signature.rds")
sc_sig <- readRDS("~/project/scCNSignature/called_signature/pcawg.result_deep/data/pcawg_deep_cn_sigs_signature.rds")
sc_sig <- readRDS("~/project/scCNSignature/called_signature/tcga.result/data/tcga_cn_sigs_signature.rds")

contrast_model <- sc_sig$Signature
max <- 7
max <- 8
result_sc <- NA

sim_nmf <- normal_nmf

for (i in c(1:max)) {
  result_cosine <- cosine(sim_nmf[,1],contrast_model[,i])
  result_sc <- append(result_sc,result_cosine)
}

result_sc <- round(result_sc[-1],3)
result_sc



draw_data <- fread("~/project/scCNSignature/ACE_NEW/sim_data/sim_result_vs_scSig_neg.txt",data.table = F)
draw_data <- fread("~/project/scCNSignature/ACE_NEW/sim_data/sim_result_vs_swgs_neg.txt",data.table = F)
draw_data <- fread("~/project/scCNSignature/ACE_NEW/sim_data/sim_result_vs_wgs_neg.txt",data.table = F)
draw_data <- fread("~/project/scCNSignature/ACE_NEW/sim_data/sim_result_vs_tcga_neg.txt",data.table = F)



draw_data <- t(draw_data)
draw_data <- cbind(result_sc,draw_data)
colnames(draw_data) <- c("0.0%","1.5%","3.0%","4.5%","6.0%","7.5%","9.0%","10.5%","12.0%","13.5%","15.0%")


library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(RColorBrewer)

colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
values <- seq(0.3, 0.8, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)


Heatmap(draw_data, name = "cosine", col = col_fun,cluster_rows = F,cluster_columns = F,column_title = "simulation CNV vs tcga signature",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", draw_data[i, j]), x, y, gp = gpar(fontsize = 10))})



data1 <- as.data.frame(matrix(NA,nrow = 11,ncol = 2))
data1$V1 <- factor(1:11)
data1$V2 <- 0:10+0.05

ggplot(data1, aes(x = V1, y = V2)) +
  geom_bar(stat = "identity", width = 0.8, fill = "gray") +
  scale_y_continuous(breaks = 0:10)+
  theme_classic()


