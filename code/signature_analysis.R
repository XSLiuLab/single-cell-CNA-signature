
remove(list = ls())

library(sigminer)
library(data.table)
library(dplyr)


SP_PCAWG <- sigprofiler_import("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/ploidy_scSExtract/", order_by_expo = TRUE, type = "all")

saveRDS(SP_PCAWG, file = "./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/sc_ploidy_solutions_sp.rds")
SP_PCAWG <- readRDS("./called_signature/test.result/data/pcawg_cn_solutions_sp.rds")

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
  left_name = "Stability (Avg Silhouette)",
  right_name = "Mean Cosine Distance",
  highlight = 7
)



pcawg_sigs <- SP_PCAWG$solution_list$S7
apply(pcawg_sigs$Exposure, 1, mean)
sig_names(pcawg_sigs)
colnames(pcawg_sigs$Signature) <- colnames(pcawg_sigs$Signature.norm) <- rownames(pcawg_sigs$Exposure) <- rownames(pcawg_sigs$Exposure.norm) <- pcawg_sigs$Stats$signatures$Signatures <- paste0("scSignature", 1:7)
sig_names(pcawg_sigs)


saveRDS(pcawg_sigs, file = "./ACE_NEW/call_signature/new_ploidy_scSignature_1Mb/result/new_sc_ploidy_sigs_signature.rds")
pcawg_sigs <- readRDS("./ACE_NEW/call_signature/new_ploidy_scSignature_1Mb/result/new_sc_ploidy_sigs_signature.rds")


pcawg_act <- list(
  absolute = get_sig_exposure(SP_PCAWG$solution_list$S7, type = "absolute"),
  relative = get_sig_exposure(SP_PCAWG$solution_list$S7, type = "relative", rel_threshold = 0),
  similarity = pcawg_sigs$Stats$samples[, .(Samples, `Cosine Similarity`)]
)

colnames(pcawg_act$absolute)[-1] <- colnames(pcawg_act$relative)[-1] <- paste0("scSignature", 1:5)
colnames(pcawg_act$similarity) <- c("sample", "similarity")

saveRDS(pcawg_act, file = "../data/pcawg_cn_sigs_CN176_activity.rds")

a <- pcawg_act$similarity
summary(pcawg_act$similarity$similarity)
hist(pcawg_act$similarity$similarity, breaks = 10, xlab = "Reconstructed similarity", main = NA)






sim <- get_sig_similarity(pcawg_sigs,pcawg_sigs)
pheatmap::pheatmap(sim$similarity, cluster_cols = F, cluster_rows = F, display_numbers = TRUE)




df_abs <- get_sig_exposure(pcawg_sigs) %>%
  tidyr::pivot_longer(cols = starts_with("scSignature"), names_to = "sig", values_to = "activity") %>%
  dplyr::mutate(
    activity = log10(activity + 1),
    sig = factor(sig, levels = paste0("scSignature", 1:7))
  )

df_rel <- get_sig_exposure(pcawg_sigs, type = "relative", rel_threshold = 0) %>%
  tidyr::pivot_longer(cols = starts_with("scSignature"), names_to = "sig", values_to = "activity") %>%
  dplyr::mutate(
    sig = factor(sig, levels = paste0("scSignature", 1:7))
  )

show_group_distribution(
  df_abs,
  gvar = "sig",
  dvar = "activity",
  order_by_fun = FALSE,
  g_angle = 90,
  ylab = "log10(activity+1)"
)




get_sample_from_sig <- function(dt, sig) {
  res <- head(dt[order(dt[[sig]], dt[[paste0("ABS_", sig)]], decreasing = TRUE)], 2L)
  res
}



pcawg_cn_obj <- readRDS("./ACE_NEW/result/scHCC_cn_1Mb_ACE.rds")
samp_summary <- pcawg_cn_obj@summary.per.sample

rel_activity <- get_sig_exposure(pcawg_sigs, type = "relative", rel_threshold = 0)
abs_activity <- get_sig_exposure(pcawg_sigs, type = "absolute")

rel_activity <- rel_activity[, lapply(.SD, function(x) {
  if (is.numeric(x)) round(x, 2) else x
})]

colnames(abs_activity) <- c("sample", paste0("ABS_scSignature", 1:7))
act <- merge(
  rel_activity, abs_activity,
  by = "sample"
)


 


dir.create("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/enrich_samples/", showWarnings = FALSE)
for (i in paste0("scSignature", 1:7)) {
  cat(paste0("Most enriched in ", i, "\n"))
  s <- get_sample_from_sig(act, i)
  print(s)
  max_cn <- max(pcawg_cn_obj@data %>% filter(sample %in% s$sample) %>% select(segVal))
  
  if (max_cn==2) {
    max_cn <- 3
  }
  
  plist <- show_cn_profile(pcawg_cn_obj,
                           samples = s$sample,
                           show_title = TRUE,
                           return_plotlist = TRUE,ylim = c(0,max_cn)
  )
  plist <- purrr::map2(plist, s$sample, function(p, x) {
    s <- samp_summary[sample == x]
    text <- paste0(
      "n_of_seg:", s$n_of_seg, "\n",
      "n_of_amp:", s$n_of_amp, "\n",
      "n_of_del:", s$n_of_del, "\n",
      "rel:", act[sample == x][[i]], "\n",
      "abs:", act[sample == x][[paste0("ABS_", i)]], "\n"
    )
    p <- p + annotate("text",
                      x = Inf, y = Inf, hjust = 1, vjust = 1,
                      label = text, color = "gray50"
    )
    p
  })
  p <- cowplot::plot_grid(plotlist = plist, nrow = 2)
  print(p)
  ggplot2::ggsave(file.path("./ACE_NEW/call_signature/ploidy_scSignature_1Mb/result/enrich_samples/", paste0(i, ".pdf")),
                  plot = p, width = 12, height = 6
  )
}




