

remove(list = ls())
setwd("~/project/scCNSignature/")

library(sigminer)
library(data.table)
library(dplyr)


SP_PCAWG <- sigprofiler_import("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/ploidy_scSExtract/", order_by_expo = TRUE, type = "all")
saveRDS(SP_PCAWG, file = "./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_solutions_sp.rds")

show_sig_number_survey(
  SP_PCAWG$all_stats %>%
    dplyr::rename(
      s = `Stability (Avg Silhouette)`,
      e = `Mean Cosine Distance`
    ) %>%
    dplyr::mutate(
      SignatureNumber = as.integer(gsub("[^0-9]", "", SignatureNumber))
    ),
  x = "SignatureNumber",
  left_y = "s", right_y = "e",
  left_name = "Stability",
  right_name = "Cosine Distance",
  highlight = 7
)



pcawg_sigs <- SP_PCAWG$solution_list$S7
apply(pcawg_sigs$Exposure, 1, mean)
sig_names(pcawg_sigs)
colnames(pcawg_sigs$Signature) <- colnames(pcawg_sigs$Signature.norm) <- rownames(pcawg_sigs$Exposure) <- rownames(pcawg_sigs$Exposure.norm) <- pcawg_sigs$Stats$signatures$Signatures <- paste0("scSignature", 1:7)
sig_names(pcawg_sigs)


saveRDS(pcawg_sigs, file = "./ACE_NEW/call_signature/new_ploidy_scSignature_1Mb/result/new_sc_ploidy_sigs_signature.rds")
pcawg_sigs <- readRDS("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")





sim <- get_sig_similarity(pcawg_sigs,pcawg_sigs)


library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
values <- seq(0.3, 0.8, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)

Heatmap(sim$similarity, name = "cosine", col = col_fun,cluster_rows = F,cluster_columns = F,column_title = "scSignature similarity",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", sim$similarity[i, j]), x, y, gp = gpar(fontsize = 10))})




sc_sig <- readRDS("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")
pcawg_sig <- readRDS("./ACE_NEW/call_signature/ploidy_pcawgSignature/result/pcawg_ploidy_sigs_signature.rds")
tcga_sigs1 <- readRDS("./called_signature/tcga.result/data/tcga_cn_sigs_signature.rds")
pcawg_deep_sigs <- readRDS("./called_signature/pcawg.result_deep/data/pcawg_deep_cn_sigs_signature.rds")

sim <- get_sig_similarity(pcawg_deep_sigs,sc_sig)
sim <- get_sig_similarity(pcawg_sig,sc_sig)
sim <- get_sig_similarity(sc_sig,tcga_sigs1)

sim <- get_sig_similarity(tcga_sigs1,sc_sig)



sim <- get_sig_similarity(pcawg_deep_sigs,tcga_sigs1)
sim <- get_sig_similarity(pcawg_sig,tcga_sigs1)
sim <- get_sig_similarity(pcawg_sig,pcawg_deep_sigs)




library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
values <- seq(0.3, 0.8, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)

Heatmap(sim$similarity, name = "cosine", col = col_fun,cluster_rows = F,cluster_columns = F,column_title = "shallow_tcga",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", sim$similarity[i, j]), x, y, gp = gpar(fontsize = 10))})





