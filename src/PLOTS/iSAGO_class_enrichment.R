setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/ip_cls_tabs/")
pacman::p_load(ggplot2, dplyr, tidyr, ggpubr)

tbFiles = list.files(pattern = ".cls_prop_iSAGOs.txt")
tbList = lapply(tbFiles, read.delim)
names(tbList) = c("PPW-1", "SAGO-1", "SAGO-2")

lst2plot = list()
for (ip in names(tbList)) {
  tmp.ip = tbList[[ip]]
  tmp.ip$norm = (tmp.ip$Raw_prop + 1) / (tmp.ip$pro_bases + 1)
  tmp.ip$IP = ip
  lst2plot[[ip]] = tmp.ip
}

df = do.call(rbind, lst2plot)

features_keep = c("exons_As", "pseudo_exons_As", "DNA_As", "LTR_As", "LINE_As",
                   "Unknown", "miRNA_S", "rRNA_S", "tRNA_S", "lincRNA_As")
df = df[df$Feature %in% features_keep, ]
df$Feature = factor(df$Feature, levels = features_keep)

df$str = ifelse(grepl("_As$", df$Feature), "sense",
                ifelse(grepl("_S$", df$Feature), "antisense", "unknown"))
df$str = factor(df$str, levels = c("sense", "antisense", "unknown"))

df$exp = "Relative enrichment by genomic category"

df$Feature_str = interaction(df$Feature, df$str, sep = "_")

col = c(
  "exons_As" = "#0B3BFA",
  "pseudo_exons_As" = "#46CCA0",
  "DNA_As" = "#00C3DB",
  "LTR_As" = "#00AB14",
  "LINE_As" = "#90DB3F",
  "Unknown" = "#D2F102",
  "miRNA_S" = "#FFDA2C",
  "rRNA_S" = "#F02765",
  "tRNA_S" = "#F75616",
  "lincRNA_As" = "#8846A3"
)

color_map = setNames(col[as.character(df$Feature)], df$Feature_str)

shape_values = c("sense" = 16, "antisense" = 17, "unknown" = 15)
shape_map = setNames(shape_values[as.character(df$str)], df$Feature_str)

labels_map = setNames(as.character(df$Feature), df$Feature_str)

norms = ggplot(df, aes(x = IP, y = log2(norm), group = Feature, color = Feature_str, shape = Feature_str)) + 
  geom_line(aes(group = Feature), color = adjustcolor("black", alpha.f = 0.7), linetype = "dashed") +
  geom_point(size = 5) +
  scale_color_manual(values = color_map, labels = labels_map) +
  scale_shape_manual(values = shape_map, labels = labels_map) +
  theme_test() +
  facet_wrap(~exp) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20, angle = 90, hjust = 0.5),
    strip.text = element_text(size = 20, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "inches"),
    legend.position = "bottom",
    legend.text = element_text(size = 14)
  ) +
  guides(color = guide_legend(title = "Feature"), shape = guide_legend(title = "Feature")) +
  xlab("") +
  ylab("log2(sRNA Density in IP)")


ggsave("iSAGO_ip_adult_class_dens_ext.pdf", device = "pdf", plot = norms,
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/", width = 8.5, height = 7, dpi = 300)
