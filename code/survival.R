
remove(list = ls())
setwd("~/project/scCNSignature/ACE_NEW/")
library(data.table)
library(dplyr)
library(sigminer)


# pcawgDeep

sc_sig <- readRDS("~/project/scCNSignature/ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")
pcawg_nmf <- readRDS("~/project/scCNSignature/called_signature/pcawgSignatureDeep/pcawg_deep_tally.rds")
nmf <- t(pcawg_nmf$nmf_matrix)
sc_sig_data <- sc_sig$Signature
result1 <- sig_fit(nmf,sc_sig)
norm_result <- apply(result1, 2, function(x){x/sum(x)})

new_fit_result <- list()
new_fit_result$fit_result <- result1
new_fit_result$norma_fit_result <- norm_result
saveRDS(new_fit_result,"~/project/scCNSignature/GDSC/survival_data/pcawgDeep_fit_result.rds")

new_fit_result <- readRDS("~/project/scCNSignature/GDSC/survival_data/pcawgDeep_fit_result.rds")
norm_result <- new_fit_result$fit_result
data1 <- data.frame(t(norm_result))
data1$sample <- rownames(data1)

pcawg_sample_tidy_info <- readRDS("~/project/scCNSignature/called_signature/pcawg.result_deep/sample_infor/pcawg_sample_tidy_info.rds")
clin_pcawgLow1 <- pcawg_sample_tidy_info %>% select(sample,donor_survival_time,donor_vital_status,wgd_status,n_amplicon_ecDNA)

clin_pcawgLow1 <- left_join(data1,clin_pcawgLow1)

library(ezcox)
data3 <- clin_pcawgLow1 %>% select(sample,donor_survival_time,donor_vital_status,starts_with("sc"))  %>% na.omit() %>% 
  filter(donor_survival_time < 3650 & donor_survival_time != 0) %>% 
  dplyr::rename(time = donor_survival_time, status = donor_vital_status) %>%
  dplyr::mutate(status = ifelse(status == "deceased", 1, 0)) %>%
  dplyr::mutate_at(vars(starts_with("sc")), ~ . / 10)



show_forest(data3, covariates = paste0("scSignature", 1:7),
            merge_models = TRUE, add_caption = FALSE, point_size = 2)




remove(list = ls())
tcgaSig <- readRDS("~/project/scCNSignature/called_signature/tcga.result/data/tcga_cn_sigs_signature.rds")
pcawg_nmf <- readRDS("~/project/scCNSignature/called_signature/pcawgSignatureDeep/pcawg_deep_tally.rds")
nmf <- t(pcawg_nmf$nmf_matrix)
sc_sig_data <- tcgaSig$Signature
result1 <- sig_fit(nmf,tcgaSig)
norm_result <- apply(result1, 2, function(x){x/sum(x)})

new_fit_result <- list()
new_fit_result$fit_result <- result1
new_fit_result$norma_fit_result <- norm_result
saveRDS(new_fit_result,"~/project/scCNSignature/GDSC/survival_data/pcawgDeep_fit_tcga_result.rds")

new_fit_result <- readRDS("~/project/scCNSignature/GDSC/survival_data/pcawgDeep_fit_tcga_result.rds")
norm_result <- new_fit_result$fit_result
data1 <- data.frame(t(norm_result))
data1$sample <- rownames(data1)

pcawg_sample_tidy_info <- readRDS("~/project/scCNSignature/called_signature/pcawg.result_deep/sample_infor/pcawg_sample_tidy_info.rds")
clin_pcawgLow1 <- pcawg_sample_tidy_info %>% select(sample,donor_survival_time,donor_vital_status,wgd_status,n_amplicon_ecDNA)

clin_pcawgLow1 <- left_join(data1,clin_pcawgLow1)

library(ezcox)
data3 <- clin_pcawgLow1 %>% select(sample,donor_survival_time,donor_vital_status,starts_with("tcga"))  %>% na.omit() %>% 
  filter(donor_survival_time < 3650 & donor_survival_time != 0) %>% 
  dplyr::rename(time = donor_survival_time, status = donor_vital_status) %>%
  dplyr::mutate(status = ifelse(status == "deceased", 1, 0)) %>%
  dplyr::mutate_at(vars(starts_with("sc")), ~ . / 10)


show_forest(data3, covariates = paste0("tcgaSignature", 1:8),
            merge_models = TRUE, add_caption = FALSE, point_size = 2)




remove(list = ls())

pacwgDeepSig <- readRDS("~/project/scCNSignature/called_signature/pcawg.result_deep/data/pcawg_deep_cn_sigs_signature.rds")
norm_result <- pacwgDeepSig$Exposure
data1 <- data.frame(t(norm_result))
data1$sample <- rownames(data1)

pcawg_sample_tidy_info <- readRDS("~/project/scCNSignature/called_signature/pcawg.result_deep/sample_infor/pcawg_sample_tidy_info.rds")
clin_pcawgLow1 <- pcawg_sample_tidy_info %>% select(sample,donor_survival_time,donor_vital_status,wgd_status,n_amplicon_ecDNA)

clin_pcawgLow1 <- left_join(data1,clin_pcawgLow1)

library(ezcox)
data3 <- clin_pcawgLow1 %>% select(sample,donor_survival_time,donor_vital_status,starts_with("pcawg")) %>% 
  filter(donor_survival_time < 3650) %>% na.omit() %>% 
  dplyr::rename(time = donor_survival_time, status = donor_vital_status) %>%
  dplyr::mutate(status = ifelse(status == "deceased", 1, 0)) %>%
  dplyr::mutate_at(vars(starts_with("pcawg")), ~ . / 10)


show_forest(data3, covariates = paste0("pcawg_deep_Signature", 1:7),
            merge_models = TRUE, add_caption = FALSE, point_size = 2)






remove(list = ls())
tcga_sigs <- readRDS("~/project/scCNSignature/called_signature/tcgaSignature3/LIHC_tcga_tally3.rds")
nmf <- t(tcga_sigs$nmf_matrix)

sc_sig <- readRDS("~/project/scCNSignature/ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")
sc_sig_data <- sc_sig$Signature
result1 <- sig_fit(nmf,sc_sig)
norm_result <- apply(result1, 2, function(x){x/sum(x)})

new_fit_result <- list()
new_fit_result$fit_result <- result1
new_fit_result$norma_fit_result <- norm_result
saveRDS(new_fit_result,"~/project/scCNSignature/GDSC/survival_data/tcga_fit_result.rds")


new_fit_result <- readRDS("~/project/scCNSignature/GDSC/survival_data/tcga_fit_result.rds")
norm_result <- new_fit_result$fit_result

data1 <- data.frame(t(norm_result))
data1$sample <- rownames(data1)

clin_tcga <- readRDS("~/project/scCNSignature/TCGA/tcga_cli.rds")
clin_tcga1 <- clin_tcga %>% select(sample,OS,OS.time)
clin_tcga1 <- left_join(data1,clin_tcga1)

drug_tcga <- fread("~/project/scCNSignature/TCGA/GDCdata/liver_drug_tcga_CN.txt")

data3 <- clin_tcga1 %>% select(sample,OS,OS.time,starts_with("sc")) %>% na.omit() %>% filter(OS.time < 3650 & OS.time != 0) %>% 
  dplyr::rename(time = OS.time, status = OS) %>%
  dplyr::mutate_at(vars(starts_with("sc")), ~ . / 10)  # %>% dplyr::filter(!(sample %in% names(table(drug_tcga$barcode))[15:35]))


show_forest(data3, covariates = paste0("scSignature", 1:7),
            merge_models = TRUE, add_caption = FALSE, point_size = 2)




remove(list = ls())
tcga_sigs <- readRDS("~/project/scCNSignature/called_signature/tcgaSignature3/LIHC_tcga_tally3.rds")
nmf <- t(tcga_sigs$nmf_matrix)

pacwgDeepSig <- readRDS("~/project/scCNSignature/called_signature/pcawg.result_deep/data/pcawg_deep_cn_sigs_signature.rds")
sc_sig_data <- pacwgDeepSig$Signature
result1 <- sig_fit(nmf,pacwgDeepSig)
norm_result <- apply(result1, 2, function(x){x/sum(x)})

new_fit_result <- list()
new_fit_result$fit_result <- result1
new_fit_result$norma_fit_result <- norm_result
saveRDS(new_fit_result,"~/project/scCNSignature/GDSC/survival_data/tcga_fit_WGS_result.rds")


new_fit_result <- readRDS("~/project/scCNSignature/GDSC/survival_data/tcga_fit_WGS_result.rds")
norm_result <- new_fit_result$fit_result

data1 <- data.frame(t(norm_result))
data1$sample <- rownames(data1)

clin_tcga <- readRDS("~/project/scCNSignature/TCGA/tcga_cli.rds")
clin_tcga1 <- clin_tcga %>% select(sample,OS,OS.time)
clin_tcga1 <- left_join(data1,clin_tcga1)



data3 <- clin_tcga1 %>% select(sample,OS,OS.time,starts_with("pcawg")) %>% na.omit() %>% filter(OS.time < 3650 & OS.time != 0) %>% 
  dplyr::rename(time = OS.time, status = OS) %>%
  dplyr::mutate_at(vars(starts_with("pcawg")), ~ . / 10)  



show_forest(data3, covariates = paste0("pcawg_deep_Signature", 1:7),
            merge_models = TRUE, add_caption = FALSE, point_size = 2)





tcgaSig <- readRDS("~/project/scCNSignature/called_signature/tcga.result/data/tcga_cn_sigs_signature.rds")
norm_result <- tcgaSig$Exposure
data1 <- data.frame(t(norm_result))
data1$sample <- rownames(data1)

clin_tcga <- readRDS("~/project/scCNSignature/TCGA/tcga_cli.rds")
clin_tcga1 <- clin_tcga %>% select(sample,OS,OS.time)
clin_tcga1 <- left_join(data1,clin_tcga1)

data3 <- clin_tcga1 %>% select(sample,OS,OS.time,starts_with("tcga")) %>% filter(OS.time < 3650 & OS.time != 0) %>% na.omit() %>% 
  dplyr::rename(time = OS.time, status = OS) %>%
  dplyr::mutate_at(vars(starts_with("tcga")), ~ . / 10)  

show_forest(data3, covariates = paste0("tcgaSignature", 1:8),
            merge_models = TRUE, add_caption = FALSE, point_size = 2)






remove(list = ls())

data <- fread("./figure/sur_data/sur_data.csv",data.table = F)

data1 <- data[1:22,2:3]
rownames(data1) <- data$sig[1:22]

data2 <- data1[1:7,]
data2 <- data1[8:14,]
data2 <- data1[15:22,]



library(dendextend)
library(circlize)
library(RColorBrewer)

col_fun = colorRamp2(c(0.5, 1, 1.5), c("blue", "white", "red"))
top_anno <- HeatmapAnnotation(sample = anno_barplot(c(139,325)))

ht_opt(RESET = TRUE)

ht_opt(heatmap_column_names_gp = gpar(fontface = "italic"), 
       heatmap_column_title_gp = gpar(fontsize = 10),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE)


ComplexHeatmap::Heatmap(as.matrix(data2),cluster_columns = F,cluster_rows = F,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if(data2[i, j] > 0)grid.text(sprintf("%.2f", data2[i, j]), x, y, gp = gpar(fontsize = 10))},
                        top_annotation = top_anno,col = col_fun,name = "Hazard ratio",row_names_side = "left")

ht_opt(RESET = TRUE)







new_fit_result <- readRDS("~/project/scCNSignature/GDSC/survival_data/pcawgDeep_fit_result.rds")
norm_result <- new_fit_result$fit_result
data1 <- data.frame(t(norm_result))
data1$sample <- rownames(data1)

pcawg_sample_tidy_info <- readRDS("~/project/scCNSignature/called_signature/pcawg.result_deep/sample_infor/pcawg_sample_tidy_info.rds")
clin_pcawgLow1 <- pcawg_sample_tidy_info %>% select(sample,donor_survival_time,donor_vital_status,wgd_status,n_amplicon_ecDNA)

clin_pcawgLow1 <- left_join(data1,clin_pcawgLow1)

data3 <- clin_pcawgLow1 %>% select(sample,donor_survival_time,donor_vital_status,starts_with("sc"))  %>% na.omit() %>% 
  filter(donor_survival_time < 3650 & donor_survival_time != 0) %>% 
  dplyr::rename(time = donor_survival_time, status = donor_vital_status) %>%
  dplyr::mutate(status = ifelse(status == "deceased", 1, 0)) %>%
  dplyr::mutate_at(vars(starts_with("sc")), ~ . / 10)



use_data <- data3 %>% dplyr::select(time,status,scSignature1,scSignature3,scSignature4)

library(survival)
library(survminer)
value <- surv_cutpoint(use_data, time = "time", event = "status", variables = "scSignature3") 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])

single_survival <- use_data %>% 
  dplyr::select(status,time,scSignature3) %>%
  dplyr::mutate(group = if_else(use_data[,"scSignature3"] > cut_off,"High","Low")) %>%
  mutate(time=round(time/30,2)) %>%
  na.omit() %>% arrange(group)

single_survival$group <- factor(single_survival$group,levels = c("High","Low"))

sfit <- survfit(Surv(time, status) ~ group, data = single_survival)
  ggsurvplot(sfit,
                pval = TRUE,
                conf.int = F,
                fun = "pct",
                xlab = "Time (Months)",
                palette = c("red", "black"),
                legend.title = ggplot2::element_blank(),
                legend.labs = c("High","Low"),
                break.time.by = 20,
                risk.table = T,
                tables.height = 0.2,
                ggtheme = theme_bw())


  

  
  new_fit_result <- readRDS("~/project/scCNSignature/GDSC/survival_data/tcga_fit_result.rds")
  norm_result <- new_fit_result$fit_result
  
  data1 <- data.frame(t(norm_result))
  data1$sample <- rownames(data1)
  
  clin_tcga <- readRDS("~/project/scCNSignature/TCGA/tcga_cli.rds")
  clin_tcga1 <- clin_tcga %>% select(sample,OS,OS.time)
  clin_tcga1 <- left_join(data1,clin_tcga1)
  
  drug_tcga <- fread("~/project/scCNSignature/TCGA/GDCdata/liver_drug_tcga_CN.txt")
  
  data3 <- clin_tcga1 %>% select(sample,OS,OS.time,starts_with("sc")) %>% na.omit() %>% filter(OS.time < 3650 & OS.time != 0) %>% 
    dplyr::rename(time = OS.time, status = OS) %>%
    dplyr::mutate_at(vars(starts_with("sc")), ~ . / 10)  # %>% dplyr::filter(!(sample %in% names(table(drug_tcga$barcode))[15:35]))
  
  
  use_data <- data3 %>% dplyr::select(time,status,scSignature1,scSignature3,scSignature4)
  
  library(survival)
  library(survminer)
  value <- surv_cutpoint(use_data, time = "time", event = "status", variables = "scSignature3") 
  cut_off <- as.numeric(value[["cutpoint"]][1, 1])
  
  single_survival <- use_data %>% 
    dplyr::select(status,time,scSignature3) %>%
    dplyr::mutate(group = if_else(use_data[,"scSignature3"] > cut_off,"High","Low")) %>%
    mutate(time=round(time/30,2)) %>%
    na.omit() %>% arrange(group)
  
  single_survival$group <- factor(single_survival$group,levels = c("High","Low"))
  
  sfit <- survfit(Surv(time, status) ~ group, data = single_survival)
  ggsurvplot(sfit,
             pval = TRUE,
             conf.int = F,
             fun = "pct",
             xlab = "Time (Months)",
             palette = c("red", "black"),
             legend.title = ggplot2::element_blank(),
             legend.labs = c("High","Low"),
             break.time.by = 20,
             risk.table = T,
             tables.height = 0.2,
             ggtheme = theme_bw())
  
  