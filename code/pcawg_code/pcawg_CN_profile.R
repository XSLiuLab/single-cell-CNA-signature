
library(sigminer)

cna <- fread("~/project/scCNSignature/ACE_NEW/pcawg_low/result/pcawg_ACE_smooth_segment.txt")


use_cn <- cna %>% filter(sample %in% c("mini-71c2374307f4acfddef67c9a185762ab.rds","mini-ca3b3edaef8df7eb4e955d9a91ddfb4d.rds"))

cn <- read_copynumber(use_cn, seg_cols = c("chromosome", "start", "end", "segVal"), genome_build = "hg38", complement = F, verbose = TRUE) 

show_cn_profile(cn, nrow = 2, ncol = 1,show_title = T,ylim = c(0,max(use_cn$segVal)))


library(sigminer)

use_cn <- fread("~/project/scCNSignature/ACE_NEW/pcawg_low/result/pcawg_ACE_smooth_segment.txt",data.table = F)

cn_obj <- read_copynumber(
  use_cn,
  seg_cols = c("chromosome", "start", "end", "segVal"),
  add_loh = F,
  max_copynumber = 1000L,
  genome_build = "hg38",
  complement = FALSE,
  genome_measure = "called"
)


saveRDS(cn_obj, file = "~/project/scCNSignature/ACE_NEW/pcawg_low/result/pcawg_cn_1Mb_ACE.rds")





