
library(sigminer)

tally_X <- readRDS("/public/slst/home/wuchx/wuchx/project/scCNSignature/ploidy_call_signature/new_ploidy_pcawgSignature/new_pcawg_tumor_tally_ploidy.rds")


sigprofiler_extract(tally_X$nmf_matrix,
  output = "/public/slst/home/wuchx/wuchx/project/scCNSignature/ploidy_call_signature/new_ploidy_pcawgSignature/ploidy_pcawgSExtract",
  range = 2:15,
  nrun = 100,
  init_method = "random",
  is_exome = FALSE,
  use_conda = TRUE,
  cores = 25
)

