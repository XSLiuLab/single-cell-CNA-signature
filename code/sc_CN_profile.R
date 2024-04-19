
library(sigminer)


HCC_infor <- fread("~/project/scCNSignature/CNS_result/scHCC_infor.csv",data.table = F)
normal <- HCC_infor %>% filter(tissue_type == "N") 
cna <- fread("~/project/scCNSignature/ACE_NEW/result/sc_ACE_smooth_segment.txt")


#Tumor
use_cn <- cna %>% filter(sample %in% c("HRR226714","HRR226275"))
#Normal
use_cn <- cna %>% filter(sample %in% c("HRR227288","HRR227310"))
use_cn <- cna %>% filter(sample %in% c("HRR227288","HRR227310"))
use_cn <- cna %>% filter(sample %in% c("HRR226720","HRR227316"))
use_cn <- cna %>% filter(sample %in% c("HRR227325","HRR227337"))


cn <- read_copynumber(use_cn, seg_cols = c("chromosome", "start", "end", "segVal"), genome_build = "hg38", complement = F, verbose = TRUE) 

show_cn_profile(cn, nrow = 2, ncol = 1,show_title = T,ylim = c(0,max(use_cn$segVal)))


library(sigminer)

use_cn <- fread("~/project/scCNSignature/ACE_NEW/result/sc_ACE_smooth_segment.txt",data.table = F)

cn_obj <- read_copynumber(
  use_cn,
  seg_cols = c("chromosome", "start", "end", "segVal"),
  add_loh = F,
  max_copynumber = 1000L,
  genome_build = "hg38",
  complement = FALSE,
  genome_measure = "called"
)


saveRDS(cn_obj, file = "~/project/scCNSignature/ACE_NEW/result/scHCC_cn_1Mb_ACE.rds")





