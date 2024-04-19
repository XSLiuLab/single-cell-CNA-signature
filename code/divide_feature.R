
library(dplyr)
library(data.table)

HCC_infor <- fread("~/project/scCNSignature/CNS_result/scHCC_infor.csv",data.table = F)
normal <- HCC_infor %>% filter(tissue_type == "N") 
sc_segment <- fread("~/project/scCNSignature/ACE_NEW/result/sc_ACE_smooth_segment.txt",data.table = F)
sc_tumor_segment <- sc_segment %>% filter(!(sample %in% normal$sample_ID))
fwrite(sc_tumor_segment,"./result/sc_tumor_ACE_smooth_segment.txt")



sc_segment <- fread("~/project/scCNSignature/ACE_NEW/result/sc_tumor_ACE_smooth_segment.txt",data.table = F)
min(sc_segment$segVal)

divide_result <- divide_feature(sc_segment)
nmf <- divide_result$nmf_matrix

divide_figure <- draw_feature(divide_result)
divide_figure

a <- divide_result$nmf_matrix

saveRDS(divide_result,"~/project/scCNSignature/ACE_NEW/result/sc_tumor_tally_ploidy.rds")

library(ggplot2)


## divide_feature
divide_feature <- function(sc_segment){
  
  sc_segment <- as.data.frame(sc_segment)
  
  ##CopyNumber
  CN <- sc_segment %>% dplyr::select(sample,segVal)
  CN$index <- rownames(CN)
  CN <- CN %>% mutate(CN = case_when(segVal < 2 ~ "0~1",
                                     segVal == 2 ~ "2",
                                     segVal == 3 ~ "3",
                                     segVal == 4 ~ "4",
                                     segVal >= 5 ~ "5+"))
  
  CN$CN <- factor(CN$CN,levels = c("0~1","2","3","4","5+"))
  CN <- CN %>% select(sample,index,CN)
  
  
  ##Length
  Length <- sc_segment
  Length$index <- rownames(Length)
  Length <- Length %>% mutate(size = round((end-start)/1000000,2))
  Length <- Length %>% mutate(Length = case_when((size < 5)~"S",
                                                 (size >=5 & size<=10)~"M",
                                                 (size >10)~"L"))
  Length <- Length %>% select(sample,index,Length)
  Length$Length <- factor(Length$Length)
  
  
  ##ContextChangeSize
  CCS <- sc_segment %>% select(sample,chromosome,segVal)
  CCS$index <- rownames(CCS)
  
  CCS$CCS <- NA
  for (i in 2:(nrow(CCS)-1)) {
    
    if (abs(CCS[i,3]-CCS[i+1,3]) > 1 || abs(CCS[i-1,3]-CCS[i,3]) > 1) {
      CCS[i,5] <- "BB"
    }else{
      CCS[i,5] <- "AA"
    }
    
  }
  
  if(abs(CCS[1,3]-CCS[2,3]) <= 1){
    CCS[1,5] <- "AA"
  }else{
    CCS[1,5] <- "BB"
  }
  
  
  if(abs(CCS[nrow(CCS),3]-CCS[nrow(CCS)-1,3]) <= 1){
    CCS[nrow(CCS),5] <- "AA"
  }else{
    CCS[nrow(CCS),5] <- "BB"
  }
  
  
  CCS <- CCS %>% select(sample,index,CCS)
  CCS$CCS <- factor(CCS$CCS)
  
  
  ####ContextChangeShape
  
  
  SCC <- sc_segment %>% select(sample,segVal,chromosome)
  SCC$index <- rownames(SCC)
  
  SCC$SCC <- NA
  for (i in 2:(nrow(SCC)-1)) {
    
    if (SCC[i,2] < SCC[i+1,2] & SCC[i-1,2] > SCC[i,2] ) {
      SCC[i,5] <- "LL"
    }else{
      if(SCC[i,2] > SCC[i+1,2] & SCC[i-1,2] < SCC[i,2]){
        SCC[i,5] <- "HH"
      }else{
        SCC[i,5] <- "OT"
      }
      
    }
    
  }
  
  
  if (SCC[1,2] < SCC[2,2]) {
    SCC[1,5] <- "LL"
  }else{
    if(SCC[1,2] > SCC[2,2]){
      SCC[1,5] <- "HH"
    }else{
      SCC[1,5] <- "OT"
    }
  }
  
  if (SCC[nrow(SCC),2] < SCC[nrow(SCC)-1,2]) {
    SCC[nrow(SCC),5] <- "LL"
  }else{
    if(SCC[nrow(SCC),2] > SCC[2,2]){
      SCC[nrow(SCC),5] <- "HH"
    }else{
      SCC[nrow(SCC),5] <- "OT"
    }
  }
  
  SCC <- SCC %>% select(sample,index,SCC)
  SCC$SCC <- factor(SCC$SCC)
  
  
  sc_tally <- list()
  
  sc_tally$data <- sc_segment
  sc_tally$feature$CN <- CN
  sc_tally$feature$Length <- Length
  sc_tally$feature$CCS <- CCS
  sc_tally$feature$SCC <- SCC
  
  
  
  
  ###feature_combine
  
  data1 <- left_join(sc_tally$feature$CN,sc_tally$feature$Length)
  data2 <- left_join(data1,sc_tally$feature$CCS)
  data3 <- left_join(data2,sc_tally$feature$SCC)
  
  
  string1 <- paste0(data3$CN,":",data3$Length,":",data3$CCS,":",data3$SCC) 
  
  a <- as.data.frame(table(string1))
  a$nor <- log10(a$Freq+1)
  
  
  
  ###all_feature
  CN_feature <- c("0~1","2","3","4","5+")
  Length_feature <- c("S","M","L")
  CCS_feature <- c("AA","BB")
  SCC_feature <- c("HH","LL","OT")
  
  
  
  combine_strings <- function(string1, string2,string3,string4) {
    result <- character()
    for (i in string1) {
      for (j in string2) {
        for (k in string3) {
          for (z in string4) {
            result <- c(result, paste0(i,":",j,":",k,":",z))
          }
        }
      }
    }
    return(result)
  }
  
  
  combinations <- combine_strings(Length_feature,SCC_feature,CN_feature,CCS_feature)
  
  
  
  
  
  ## nmf_matrix
  data3$conbin <- paste0(data3$Length,":",data3$SCC,":",data3$CN,":",data3$CCS)
  
  nmf_matrix <- matrix(nrow = nrow(table(sc_segment$sample)),ncol = length(combinations),data = 0)
  rownames(nmf_matrix) <- names(table(sc_segment$sample))
  colnames(nmf_matrix) <- combinations
  
  for (sample_id in rownames(nmf_matrix)) {
    all_feature <- data3 %>% filter(sample == sample_id)
    single_sample <- as.data.frame(table(all_feature$conbin))
    
    nmf_matrix[sample_id,as.character(single_sample$Var1)] <- single_sample$Freq
    
  }
  
  
  
  sc_tally$nmf_matrix <- nmf_matrix
  
  
  return(sc_tally)
  
}


## feature 
draw_feature <- function(sc_tally){

## draw figure
mat <- sc_tally$nmf_matrix %>% t()

mat <- mat %>%
  rowSums(na.rm = TRUE) %>%
  as.matrix()
colnames(mat) <- "Total"


Sig <-mat

mat <- as.data.frame(Sig)
mat$context <- rownames(mat)

mat$base <- substr(mat$context, 1, 4)
mat <- tidyr::gather(mat, class, signature, -c("context", "base"))

mat <- dplyr::mutate(mat, context = factor(.data$context,levels = unique(mat$context)),
                     base = factor(.data$base, levels = unique(mat$base)),
                     class = factor(class, levels = colnames(Sig)))





mat$signature <- log10(mat$signature+1)
col_df <- mat %>%
  dplyr::filter(.data$class == .data$class[1]) %>%
  dplyr::group_by(.data$base) %>%
  dplyr::summarise(N = dplyr::n())



len_base <- length(unique(mat$base))
palette <- c(2:4, 7:9, 12:14)




y_expand <- 1
scale <- 1
base_size <- 12
x_label_angle <- 90
x_label_vjust <- 0.5
x_label_hjust <- 1



.theme_ss <- theme_bw(
  base_size = base_size,
  base_family = "Arial"
)+
  theme(
    axis.text.x = element_text(
      angle = x_label_angle, vjust = x_label_vjust,
      hjust = x_label_hjust, size = (base_size - 4) * scale,
      color = "black",
      face = "bold",
      family = "Arial"
    ),
    axis.text.y = element_text(
      hjust = 0.5,
      size = base_size * scale,
      color = "black"
    ),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold")
  )


.theme_ss <- .theme_ss + theme(
  axis.line = element_line(linewidth = 0.3, colour = "black"),
  panel.spacing.x = unit(0, "line"),
  strip.background.x = element_rect(color = "white"),
  strip.background.y = element_blank(),
  strip.text.x = element_text(
    color = "black",
    face = "bold"
  ),
  strip.text.y = element_text(
    size = 12,
    vjust = 1,
    color = "black",
    face = "bold",
    angle = 0
  )
)




p1 <- ggplot(mat) +
  geom_bar(aes_string(x = "context", y = "signature", fill = "base"),
           stat = "identity", position = "identity",
           colour = "white", width = 0.7) +
  scale_fill_manual(values = palette)+
  facet_grid(class ~ base, scales = "free", space = "free_x")+
  scale_x_discrete(breaks = mat$context,labels = substring(mat$context, first = 6))+
  guides(fill = "none")+
  .theme_ss+
  theme(plot.margin = margin(30 * y_expand, 2, 2, 2, unit = "pt"))+
  ylab("log10(count +1)")+xlab("Components")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


return(p1)
}



