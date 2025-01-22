
library(sigminer)

tally_X <- readRDS("data/NMF Matrix/scWES_tally_ploidy")

sigprofiler_extract(tally_X$nmf_matrix,
  output = "./",
  range = 2:15,
  nrun = 100,
  init_method = "random",
  is_exome = FALSE,
  use_conda = TRUE,
  cores = 25
)

