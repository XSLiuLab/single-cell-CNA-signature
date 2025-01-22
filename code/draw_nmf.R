

setwd("~/project/scCNSignature/")
remove(list = ls())

library(dplyr)
library(data.table)
library(ggplot2)
library(sigminer)

sc_tally <- readRDS("/result/sc_tally.rds")



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
  base_family = "Helvetica"
)+
  theme(
    axis.text.x = element_text(
      angle = x_label_angle, vjust = x_label_vjust,
      hjust = x_label_hjust, size = (base_size - 4) * scale,
      color = "black",
      face = "bold",
      family = "Helvetica"
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



p <- ggplot(mat) +
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


p

ggplot2::ggsave("/figure/scNMF_profile.pdf",plot = p, width = 12, height = 6)

