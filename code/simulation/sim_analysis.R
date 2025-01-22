

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

nature_sim <- readRDS("/data/sim_data/Nature_simulation_nmf.rds")
nature_nmf <- t(nature_sim)


sc_sig <- readRDS("~/project/scCNSignature/ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")
sc_sig_data <- sc_sig$Signature
result1 <- sig_fit(nature_nmf,sc_sig)
norm_result <- apply(result1, 2, function(x){x/sum(x)})

new_fit_result <- list()
new_fit_result$fit_result <- result1
new_fit_result$norma_fit_result <- norm_result

nature_fit <- data.frame(t(new_fit_result$norma_fit_result))
nature_fit$Samples <- rownames(nature_fit)

sample_infor <- fread("/data/sim_data/Out_2_Overview_simulation_20each_N240.txt")
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





nature_sim <- readRDS("/data/sim_data/Nature_simulation_nmf.rds")
nature_nmf <- t(nature_sim)


sc_sig <- readRDS("/tcga_cn_sigs_signature.rds")
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

sample_infor <- fread("/data/sim_data/Out_2_Overview_simulation_20each_N240.txt")
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





nature_sim <- readRDS("/data/sim_data/Nature_simulation_nmf.rds")
nature_nmf <- t(nature_sim)

sc_sig <- readRDS("./pcawg_deep_cn_sigs_signature.rds")
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

sample_infor <- fread("/output/Out_2_Overview_simulation_20each_N240.txt")
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

