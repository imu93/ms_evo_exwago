setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/rrna_norm/")
pacman::p_load(dplyr, ggplot2)
files = list.files(pattern = ".*deTbl.txt")
files = files[!grepl("caeno", files)]
de_tbls = lapply(files, read.delim)
names(de_tbls) = c("A_ceylanicum", "H_bakeri", "H_polygyrus", "N_brasiliensis",
                   "T_circumcincta")
de_tbls = de_tbls[c("H_bakeri", "H_polygyrus", "N_brasiliensis", "T_circumcincta", "A_ceylanicum")]

n_lst = list()
for (i in names(de_tbls)) {
  tmp.tbl = de_tbls[[i]]
  tmp.tbl$exp = paste(sub("_", "\\. ", i), "exWAGO vs adult total")
  sp = sub("_", "\\. ", i)
  tmp.tbl$exp <- factor(tmp.tbl$exp)
  
  tmp.tbl$exp = sub("_", ". ", i)
  #tmp.tbl$exp <- factor(tmp.tbl$exp)
  #levels(tmp.tbl$exp) <- c(expression(paste(italic("T. circumcincta "), "exWAGO IP vs Adult total")))
  
  tmp.tbl$ord = ifelse(tmp.tbl$col == "#C0C0C0", 1, 2)
  tmp.tbl = tmp.tbl[order(tmp.tbl$ord),]
  tmp.tbl = tmp.tbl[,c("Ave_CPM", "logFC", "cl", "col", "ord", "exp")]
  n_lst[[i]] = tmp.tbl
}

lst2plot= list()
for (i in names(n_lst)) {
  tmp.tbl = n_lst[[i]]
  #tmp.tbl$exp = sub("T. circumcincta", sub("_", "\\. ", i), tmp.tbl$exp)
  #tmp.tbl$exp = factor(tmp.tbl$exp)
  tmp.tbl$sp = sub("_", ". ", i)
  lst2plot[[i]] = tmp.tbl
}
df2plot = do.call(rbind, lst2plot)
df2plot$exp %>%  table()

#levels(df2plot$exp)=  c("paste(italic(\"H. bakeri \"), \"exWAGO IP vs Adult total\")" ,  "paste(italic(\"H. polygyrus \"), \"exWAGO IP vs Adult total\")",
#                        "paste(italic(\"N. brasiliensis \"), \"exWAGO IP vs Adult total\")", "paste(italic(\"T. circumcincta \"), \"exWAGO IP vs Adult total\")",
#                        "paste(italic(\"A. ceylanicum \"), \"exWAGO IP vs Adult total\")")

df2plot$sp = factor(df2plot$sp, levels = c("H. bakeri", "H. polygyrus", "N. brasiliensis", "T. circumcincta", "A. ceylanicum"))
df2plot$norm = factor("rRNA-based")
p1 = ggplot(df2plot, aes(x=Ave_CPM, y=logFC, colour = col)) +
  geom_point() +  xlab(bquote(~log[2]~' CPM')) + theme_test() +
  ylab(bquote(~log[2]~ 'FC IP/Total')) + theme(axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size =20),axis.text.y = element_text(size = 20)) +
  geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4)  +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) +  
  scale_color_manual(values = adjustcolor(c("#155AA3", "#C0C0C0","red"), alpha.f = .6),
                     labels=factor(c('exWAGO enriched', 'non-DE', "rRNA"), levels = c("non-DE","exWAGO enriched", "rRNA")), name=NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) + theme(legend.text = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size=10)))  +
  facet_grid(norm ~ sp) +
  theme(legend.position = "bottom", legend.text = element_text(size = 20)) +
  theme( strip.text = element_text(size = 18, face = "italic", colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) +
  theme(panel.spacing = unit(1.3, "lines"))

ggsave("nematode_exwago_vs_total_rRNA_norm.png", plot = p1, device = "png",
       "~/storage/Data/ms_evo_exwago/analysis/figures/",
       dpi = 300, width = 16, height = 6)
