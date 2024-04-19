

setwd("~/project/scCNSignature/")
remove(list = ls())
library(sigminer)

sc_sig <- readRDS("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")
pcawg_sig <- readRDS("./ACE_NEW/call_signature/ploidy_pcawgSignature/result/pcawg_ploidy_sigs_signature.rds")


tcga_sigs1 <- readRDS("./called_signature/tcga.result/data/tcga_cn_sigs_signature.rds")
pcawg_deep_sigs <- readRDS("./called_signature/pcawg.result_deep/data/pcawg_deep_cn_sigs_signature.rds")

sim <- get_sig_similarity(pcawg_sig,sc_sig)

sim <- get_sig_similarity(sc_sig,tcga_sigs1)
sim <- get_sig_similarity(pcawg_deep_sigs,sc_sig)

sim <- get_sig_similarity(pcawg_deep_sigs,pcawg_sig)
sim <- get_sig_similarity(pcawg_sig,tcga_sigs1)


library(pheatmap)
pheatmap::pheatmap(sim$similarity, cluster_cols = F, cluster_rows = F, display_numbers = TRUE)



sc_sig <- readRDS("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_sigs_signature.rds")
pcawg_sig <- readRDS("./ACE_NEW/call_signature/ploidy_pcawgSignature/result/pcawg_ploidy_sigs_signature.rds")
tcga_sigs1 <- readRDS("./called_signature/tcga.result/data/tcga_cn_sigs_signature.rds")
pcawg_deep_sigs <- readRDS("./called_signature/pcawg.result_deep/data/pcawg_deep_cn_sigs_signature.rds")
sim <- get_sig_similarity(pcawg_deep_sigs,pcawg_sig)
sim <- get_sig_similarity(pcawg_sig,tcga_sigs1)

library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
values <- seq(0.3, 0.8, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)

Heatmap(sim$similarity, name = "cosine", col = col_fun,cluster_rows = F,cluster_columns = F,column_title = "TCGA signature VS pcawgLow signature",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", sim$similarity[i, j]), x, y, gp = gpar(fontsize = 10))})





