setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/ip_cls_tabs/")
pacman::p_load(ggplot2, ggforce, dplyr)
tabs = list.files(pattern = "exted.txt")
names(tabs) = c("A_ceylanicum", "PPW-1", "SAGO-1", "SAGO-2", "H_bakeri", "H_polygyrus",
                "N_brasiliensis", "T_circumcincta")

lst =  lapply(tabs, read.delim)
lst =  lapply(lst, function(x){
  tab = x[grepl("_S|_As", rownames(x)),]
  return(tab)
})

lst2plot = list()
for (i in names(lst)) {
  x = lst[[i]]
  x$IP = i
  xs = x[grepl("_S", rownames(x)),]
  xas = x[grepl("_As", rownames(x)),]
  colnames(xas)  = paste0(colnames(x), "_As")
  x = cbind(xas,xs)

  lst2plot[[i]] = x
}
names(lst2plot) = NULL
df = do.call(rbind, lst2plot)
df$cat = sub("_As", "", df$Feature_As)
df$cat = factor(df$cat, levels = c("exons", "introns", "DNA", "LTR", "LINE", "miRNA", "rRNA", "tRNA", "lincRNA"))
cols = c("#0B3BFA", "#007CFF","#00C3DB", "#00AB14", "#90DB3F", "#FFDA2C", "#F02765" , "#F75616", "#8846A3")
df$exp = "Strand bias in exWAGO guide production"
df$group <- interaction(df$cat, df$IP, sep = "_")

str_pt = ggplot(df, aes(x = Raw_prop_As, y= Raw_prop, color = cat)) + geom_point(size=4) +
  geom_smooth(method = "lm", se = F,level=0.90) +
  scale_color_manual(values = cols, name="") +
  theme_test() + 
  theme(axis.title.x = element_text(face = "bold", size = 18), 
        axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5)) +
  ylab("% of Sense sRNAs") + xlab("% of Antisense sRNAs") +
  facet_wrap(~exp) +
  theme( strip.text = element_text(size = 18, colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12)) +
  annotation_custom(grid::linesGrob(gp = grid::gpar(col = 'red', lty = 2, lwd = 3))) + ylim(c(0, 60))

str_pt + facet_zoom(xlim = c(0, 10), ylim=c(0, 10))






df =  df %>%
  mutate(StrandBias = Raw_prop - Raw_prop_As)

df_summary = df %>%
  group_by(cat) %>%
  summarise(mean_bias = mean(StrandBias, na.rm = TRUE))


cols = c("#0B3BFA", "#007CFF","#00C3DB", "#00AB14", "#90DB3F", "#FFDA2C", "#F02765" , "#F75616", "#8846A3"
)

df_summary$exp = "Strand bias in exWAGO guide production"
str_b = ggplot() +
  geom_bar(data = df_summary,
           aes(x = cat, y = mean_bias, fill = cat),
           stat = "identity", width = 0.6, alpha = 0.8) +
  geom_jitter(data = df,
              aes(x = cat, y = StrandBias, color = cat),
              width = 0.2, size = 3, alpha = 0.9) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_color_manual(values = cols, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20,),
    axis.title.x = element_text(size = 20, )
  ) +
  ylab("Strand Bias (% Sense - % Antisense)") +
  xlab("Category") +
  ylim(c(-75, 40)) +
  facet_wrap(~exp) +
  theme( strip.text = element_text(size = 18, colour = "black"))

ggsave(path = "~/storage/Data/ms_evo_exwago/analysis/figures/", plot = str_b, device = "pdf", filename = "nematode_strandplot.pdf", 
       dpi=300, width = 8, height = 6)
