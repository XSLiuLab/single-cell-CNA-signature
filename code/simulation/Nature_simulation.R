
setwd("~/project/scCNSignature/ACE_NEW/nature_sim/")
remove(list = ls())
## This example script creates 240 copy number profiles from 5 mutational processes (HRD/LST, ecDNA, chromosome missegregation, early WGD and late WGD)

## Libraries
# Main
library(this.path)
library(data.table)
library(reshape2)
library(stringi)
library(GenomicRanges)
# For parallelising the simulation
library(foreach)
library(doMC)
# For plotting
library(ggplot2)
library(ggthemes)
library(lemon)
library(patchwork)
library(svglite)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5),
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))


#### Paths and files
thisPath = dirname(this.path())
OUT=file.path(dirname(thisPath), "output/")

# For chromosome and gap sizes
prePath=file.path(dirname(thisPath), "input/")
fileChrSizes = file.path(prePath, "hg19.chrom.sizes.txt")
fileCentromeres = file.path(prePath, "gap_hg19.txt")


## Samples per signature
N=200
# Plot signature definitions a priori
STEP1=TRUE
# Multiple signatures per sample, some of them on the same chromosomes
STEP2=TRUE
# Plot profiles
STEP3=TRUE

# File suffix
NAME="20each_N240"

# Functions
source(file.path(thisPath, "CINGenomeSimulation_Functions.R"))

## Load and prepare files
dfChrSizes = loadChromSizes(fileChrSizes, rmChr = TRUE)
dfCentromeres = loadCentromeres(fileCentromeres, rmChr = TRUE)




###### STEP 1: Signature definitions plotted apriori
if(STEP1) {
  
  #### Plot the features of the simulated signatures
  vFeatComps = paste0(rep(c("segsize", "changepoint", "bp10MB", "bpchrarm", "osCN"), c(22, 10, 3, 5, 3)),
                      c(1:22, 1:10, 1:3, 1:5, 1:3))
  dfSigs = data.frame("component" = factor(vFeatComps, levels = vFeatComps), 
                      # "feature" = factor(rep(c("segsize", "changepoint", "bp10MB", "bpchrarm", "osCN"), c(22, 10, 3, 5, 3))),
                      "HRDLST" = c(rep(0, 8), rep(1,7), rep(0,7), # Segsize
                                   c(0,0,1,1, rep(0, 6)), # Changepoint
                                   c(0,1,0), # bp10MB
                                   c(0,1,1,0,0), # bpchrarm
                                   c(1,0,1)), # osCN
                      "ECDNA" = c(rep(0,4), 1, 1, 1, rep(0, 15), # Segsize
                                  c(rep(0, 8), c(1,1)), # Changepoint
                                  c(0, 1, 1), # bp10MB
                                  c(0, 1, 1, 0, 0), # bpchrarm
                                  c(1,1,0)), # osCN
                      "CHR" = c(rep(0, 15), rep(1,7), # Segsize
                                c(0,0, 1, 1, rep(0, 6)), # Changepoint
                                c(1, 0, 0), # bp10MB
                                c(1, rep(0,4)), # bpchrarm
                                c(1,rep(0,2))), # osCN
                      "WGDearly" = c(rep(0, 15), rep(1,7), # Segsize
                                     c(rep(0,4),1,1,rep(0, 4)), # Changepoint
                                     c(1, 0, 0), # bp10MB
                                     c(1, rep(0,4)), # bpchrarm
                                     c(1,rep(0,2))),
                      "WGDlate" = c(rep(0, 4), rep(1,18), # Segsize
                                    c(rep(0,3),1,1,1,rep(0, 4)), # Changepoint
                                    c(0, 1, 1), # bp10MB
                                    c(0,1,1,0,0), # bpchrarm
                                    c(1, 1, 1))) # osCN
  dfMelt = melt(dfSigs)
  
  ## Rename sigs
  dfMelt$variable = factor(as.character(dfMelt$variable), levels = c("HRDLST", "ECDNA", "CHR", "WGDearly", "WGDlate"),
                           labels = c("HRD/LST", "ecDNA", "CHR", "WGD (early)", "WGD (late)"))
  
  
  dfMelt$feature = gsub('[0-9]+', '', dfMelt$component)
  dfMelt$feature[ dfMelt$feature == "bpMB" ] = "bp10MB"
  dfMelt$feature = factor(dfMelt$feature, levels = c("segsize", "changepoint",
                                                     "bp10MB", "bpchrarm",
                                                     "osCN"))
  
  
  p0 = ggplot(dfMelt, aes(x = component, y = value, fill = feature))+
    geom_bar(stat = "identity") + facet_wrap(. ~ variable, ncol = 1) + 
    coord_capped_cart(left = "both", bottom = "both") +
    labs(x = "Feature components", y = "Simulated weight") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none") 
  
  ggsave(file.path(OUT, "Out_1_Simulated_signature_defs.png"), 
         p0, bg = "white", width = 12, height = 10, units = "cm")
  ggsave(file.path(OUT, "Out_1_Simulated_signature_defs.svg"), 
         p0, bg = "white", width = 12, height = 10, units = "cm")
  
}


###### STEP 2: Multiple mutational processes per genome and per chromosome
if(STEP2) {
  
  #### Produce samples - Step 1: Sample mutational processes
  set.seed(10)
  
  ## Three signatures with all of them active in all signatures
  # One dominant sig
  lSigs3LST = sampleMutActivity(N = N, SIGS = c("LST", "ECDNA", "CHR"), VALS = c(8, 2, 4), REST = "WGD", dfChrSizes)
  lSigs3ECDNA = sampleMutActivity(N = N, SIGS = c("LST", "ECDNA", "CHR"), VALS = c(2, 4, 4), REST = "WGD", dfChrSizes)
  lSigs3CHR = sampleMutActivity(N = N, SIGS = c("LST", "ECDNA", "CHR"), VALS = c(2, 2, 8), REST = "WGD", dfChrSizes)
  # Two dominant sigs
  lSigs3LSTCHR = sampleMutActivity(N = N, SIGS = c("LST", "ECDNA", "CHR"), VALS = c(8, 2, 8), REST = "WGD", dfChrSizes)
  lSigs3LSTECDNA = sampleMutActivity(N = N, SIGS = c("LST", "ECDNA", "CHR"), VALS = c(8, 4, 4), REST = "WGD", dfChrSizes)
  lSigs3CHRECDNA = sampleMutActivity(N = N, SIGS = c("LST", "ECDNA", "CHR"), VALS = c(2, 4, 8), REST = "WGD", dfChrSizes)
  
  # With early WGD and 2 signatures (changing in dominance)
  X=0.2
  lSigs2WGDLSTECDNA = sampleMutActivity(N = N*(1-X), SIGS = c("LST", "ECDNA"), VALS = c(8, 2), REST = "WGD", dfChrSizes)
  lSigs2WGDECDNALST = sampleMutActivity(N = N*(1-X), SIGS = c("LST", "ECDNA"), VALS = c(2, 4), REST = "WGD", dfChrSizes)
  lSigs2WGDECDNACHR = sampleMutActivity(N = N*(1-X), SIGS = c("ECDNA", "CHR"), VALS = c(4, 4), REST = "WGD", dfChrSizes)
  lSigs2WGDCHRECDNA = sampleMutActivity(N = N*(1-X), SIGS = c("ECDNA", "CHR"), VALS = c(2, 8), REST = "WGD", dfChrSizes)
  lSigs2WGDLSTCHR = sampleMutActivity(N = N*(1-X), SIGS = c("LST", "CHR"), VALS = c(8, 4), REST = "WGD", dfChrSizes)
  lSigs2WGDCHRLST = sampleMutActivity(N = N*(1-X), SIGS = c("LST", "CHR"), VALS = c(2, 8), REST = "WGD", dfChrSizes)
  
  # With late WGD and 2 signatures (changing in dominance)
  lSigs2LateLSTECDNA = sampleMutActivity(N = N*X, SIGS = c("LST", "ECDNA"), VALS = c(8, 2), REST = "WGD", dfChrSizes)
  lSigs2LateECDNALST = sampleMutActivity(N = N*X, SIGS = c("LST", "ECDNA"), VALS = c(2, 4), REST = "WGD", dfChrSizes)
  lSigs2LateECDNACHR = sampleMutActivity(N = N*X, SIGS = c("ECDNA", "CHR"), VALS = c(4, 4), REST = "WGD", dfChrSizes)
  lSigs2LateCHRECDNA = sampleMutActivity(N = N*X, SIGS = c("ECDNA", "CHR"), VALS = c(2, 8), REST = "WGD", dfChrSizes)
  lSigs2LateLSTCHR = sampleMutActivity(N = N*X, SIGS = c("LST", "CHR"), VALS = c(8, 4), REST = "WGD", dfChrSizes)
  lSigs2LateCHRLST = sampleMutActivity(N = N*X, SIGS = c("LST", "CHR"), VALS = c(2, 8), REST = "WGD", dfChrSizes)
  
  ## For names
  dfOverview = data.frame("Signature" = c(rep("LST", N), rep("ecDNA", N), rep("Chr", N),
                                          rep("LST_Chr", N), rep("LST_ecDNA", N), rep("Chr_ecDNA", N),
                                          rep("WGDearly_LST_ecDNA", N*(1-X)), rep("WGDearly_ecDNA_LST", N*(1-X)), rep("WGDearly_ecDNA_Chr", N*(1-X)),
                                          rep("WGDearly_Chr_ecDNA", N*(1-X)), rep("WGDearly_LST_Chr", N*(1-X)), rep("WGDearly_Chr_LST", N*(1-X)),
                                          rep("WGDlate_LST_ecDNA", N*X), rep("WGDlate_ecDNA_LST", N*X), rep("WGDlate_ecDNA_Chr", N*X),
                                          rep("WGDlate_Chr_ecDNA", N*X), rep("WGDlate_LST_Chr", N*X), rep("WGDlate_Chr_LST", N*X)))
  
  
  #### Produce samples - Step 2: Create genome profiles
  lAllMutsDiploid1 = c(lSigs3LST, lSigs3ECDNA)
  lAllMutsDiploid2 = c(lSigs3CHR, lSigs3LSTCHR)
  lAllMutsDiploid3 = c(lSigs3LSTECDNA, lSigs3CHRECDNA)
  lAllMutsTetraploid1 = c(lSigs2WGDLSTECDNA, lSigs2WGDECDNALST)
  lAllMutsTetraploid2 = c(lSigs2WGDECDNACHR, lSigs2WGDCHRECDNA)
  lAllMutsTetraploid3 = c(lSigs2WGDLSTCHR, lSigs2WGDCHRLST)
  lAllMutsLateWGD = c(lSigs2LateLSTECDNA, lSigs2LateECDNALST, lSigs2LateECDNACHR, lSigs2LateCHRECDNA, lSigs2LateLSTCHR, lSigs2LateCHRLST)
  
  ## Just so chromosome missegregation always happens when it is called (usually 20% probability only)
  LISTVARS=list("PROB" = 1, "OVERLAP" = TRUE)
  
  ## Generate profiles - parallelised to make it faster
  CORES = ifelse(X == 0, 6, 7)
  registerDoMC(CORES)
  lStep3 = foreach(i=1:7) %dopar% {
    if(i == 1){
      simulateMultipleGenomes(N=1, WHICHSIG = lAllMutsDiploid1, CHRSIZES = dfChrSizes, CENTSIZES = dfCentromeres,
                              SEED = 12, NORMALCP = 2, LISTVARS = LISTVARS, PROGRESS = TRUE)
    } else if (i == 2) {
      simulateMultipleGenomes(N=1, WHICHSIG = lAllMutsDiploid2, CHRSIZES = dfChrSizes, CENTSIZES = dfCentromeres,
                              SEED = 13, NORMALCP = 2, LISTVARS = LISTVARS, PROGRESS = TRUE)
    } else if (i == 3) {
      simulateMultipleGenomes(N=1, WHICHSIG = lAllMutsDiploid3, CHRSIZES = dfChrSizes, CENTSIZES = dfCentromeres,
                              SEED = 14, NORMALCP = 2, LISTVARS = LISTVARS, PROGRESS = TRUE)
    } else if (i == 4) {
      simulateMultipleGenomes(N=1, WHICHSIG = lAllMutsTetraploid1, CHRSIZES = dfChrSizes, CENTSIZES = dfCentromeres,
                              SEED = 15, NORMALCP = 4, LISTVARS = LISTVARS, PROGRESS = TRUE)
    } else if (i == 5) {
      simulateMultipleGenomes(N=1, WHICHSIG = lAllMutsTetraploid2, CHRSIZES = dfChrSizes, CENTSIZES = dfCentromeres,
                              SEED = 16, NORMALCP = 4, LISTVARS = LISTVARS, PROGRESS = TRUE)
    } else if (i == 6) {
      simulateMultipleGenomes(N=1, WHICHSIG = lAllMutsTetraploid3, CHRSIZES = dfChrSizes, CENTSIZES = dfCentromeres,
                              SEED = 17, NORMALCP = 4, LISTVARS = LISTVARS, PROGRESS = TRUE)
    } else if (i == 7) {
      simulateMultipleGenomes(N=1, WHICHSIG = lAllMutsLateWGD, CHRSIZES = dfChrSizes, CENTSIZES = dfCentromeres,
                              SEED = 18, NORMALCP = 2, LISTVARS = LISTVARS, PROGRESS = TRUE)
    }
  }
  
  # Late WGD
  if(X>0) {
    lLate = lStep3[[7]]
    multiplier = rep(2, length(lLate))
    # multiplier[ sample(1:length(lLate), round(length(lLate)*0.2)) ] = 3
    lLateWGD = lapply(1:length(lLate), function(i) {
      lLate[[i]]$segVal = as.numeric(lLate[[i]]$segVal) * multiplier[i]
      return(lLate[[i]])
    })
    lStep3[[7]] = lLateWGD
  }
  
  # Unify runs into one long list
  lStep3Un = unlist(lStep3)
  
  ## Somehow names disappear - quickfix for now
  listNames = sapply(lStep3Un, function(x) unique(x$sample))
  names(lStep3Un) = listNames
  
  ## Add to names data frame for documentation
  dfOverview$Samples = listNames
  
  lDFStep3 = convertGRListToDF(lStep3Un)
  saveRDS(lDFStep3, file.path(OUT, paste0("Out_2_Segmentation_tables_", NAME, ".rds")))
  saveRDS(dfOverview, file.path(OUT, paste0("Out_2_Overview_simulation_", NAME, ".rds")))
  write.table(dfOverview, file.path(OUT, paste0("Out_2_Overview_simulation_", NAME, ".txt")), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE )
  
}

# STEP 3: Plot example profies
if(STEP3) {
  ## Plot sample profiles for all combinations
  lStep3Plots = plotAllSamples(lDFStep3, dfChrSizes)
  examplePlot3 = wrap_plots(lStep3Plots[c(1, 101, 122, 231)], ncol = 2)
  ggsave(file.path(OUT, paste0("Out_2_Example_plots_", NAME, ".png")), examplePlot3, scale = 1.5, units = "cm", width = 16, height = 12,
         bg = "white")
  ggsave(file.path(OUT, paste0("Out_2_Example_plots_", NAME, ".svg")), examplePlot3, scale = 1.5, units = "cm", width = 16, height = 12,
         bg = "white")
  
}






