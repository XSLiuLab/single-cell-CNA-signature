library(data.table)
library(QDNAseq)
library(QDNAseq.hg38)
library(ACE)

filename <- list.files("/public/slst/home/wuchx/wuchx/project/scCNSignature/use_data")

for(i in filename[1:100]){
  bins <- getBinAnnotations(binSize=100,genome="hg38")
  readCounts <- binReadCounts(bins,bamfile=paste0("/public/slst/home/wuchx/wuchx/project/scCNSignature/use_data/",i,"/",i,".rdup.bam"))
  readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="log2")
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  save(copyNumbersSegmented,file=paste0("/public/slst/home/wuchx/wuchx/project/scCNSignature/CNS_result100/QDNAseq_result/",i,".rds"))
  m<-ACEcall(copyNumbersSegmented, QDNAseqobjectsample = TRUE, cellularity = 2.5,
  ploidy = 2, standard =1, plot = FALSE, title ="",
  qcutoff = -3, subcutoff = 0.2, trncname = FALSE,
  cap = 12, bottom = 0, chrsubset = 1:22, onlyautosomes = TRUE,
  sgc = c())
  save(m,file=paste0("/public/slst/home/wuchx/wuchx/project/scCNSignature/CNS_result100/ACE_result/",i,".rds"))
}


