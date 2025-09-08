setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/tmm/cele_wagos/")
pacman::p_load(dplyr, ggplot2)
files = list.files(pattern = ".*deTbl.txt")
files = files[grepl("caeno", files)]
files = files[!grepl("alg5", files)]
de_tbls = lapply(files, read.delim)
names(de_tbls) = c("ALG-1", "ALG-2", "ALG-3", "ALG-4", 
                   "CSR-1", "ERGO-1", "HRDE-1", "NRDE-3", 
                   "PPW-1", "PPW-2", "PRG-1", "RDE-1", "SAGO-1", 
                   "SAGO-2", "VSRA-1", "WAGO-1", "WAGO-10", "WAGO-4")

n_lst = list()
for (i in names(de_tbls)) {
  tmp.tbl = de_tbls[[i]]
  head(tmp.tbl)
  tmp.tbl$exp = i
  tmp.tbl$exp <- factor(tmp.tbl$exp)
  tmp.tbl$col = ifelse(grepl("rRNA_S", rownames(tmp.tbl)), "red", tmp.tbl$col)
  tmp.tbl$ord = ifelse(tmp.tbl$col == "#C0C0C0", 1, 2)
  tmp.tbl$ord = ifelse(tmp.tbl$col == "red", 3, tmp.tbl$ord)
  tmp.tbl = tmp.tbl[order(tmp.tbl$ord),]
  tmp.tbl$sp = i
  tmp.tbl = tmp.tbl[,c("Ave_CPM", "logFC", "cl", "col", "ord", "exp","sp")]
  n_lst[[i]] = tmp.tbl
}


lst2plot= list()
for (i in names(n_lst)) {
  tmp.tbl = n_lst[[i]]
  lst2plot[[i]] = tmp.tbl
}
df2plot = do.call(rbind, lst2plot)
df2plot$exp %>%  table()



df2plot$sp = factor(df2plot$sp, levels = c("ALG-1", "ALG-2", "ALG-3", "ALG-4", 
                                           "CSR-1", "ERGO-1", "HRDE-1", "NRDE-3", 
                                           "PPW-1", "PPW-2", "PRG-1", "RDE-1", "SAGO-1", 
                                           "SAGO-2", "VSRA-1", "WAGO-1", "WAGO-10", "WAGO-4"))

df2plot$col = factor(df2plot$col, levels = c("#5E2129", "#C0C0C0", "red"))
df2plot$norm = factor("TMM")

p1 = ggplot(df2plot, aes(x=Ave_CPM, y=logFC, colour = col)) +
  geom_point() +  xlab(bquote(log[2]~' CPM')) + theme_test() +
  ylab(bquote(~log[2]~ 'FC IP/Total')) +
  theme(axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20),axis.text.y = element_text(size = 20)) +
  geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
  facet_wrap(~ exp, labeller = label_parsed, nrow=7) +
  scale_color_manual(values = adjustcolor(c("#155AA3","#C0C0C0", "red"), alpha.f = .6), 
                     labels=c('IP', 'not-DE','rRNA_S'), name=NULL) +
  theme(legend.position = "bottom", legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=10))) +
  theme( strip.text = element_text(size = 18,colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) 
p1

ggsave("wagos_vs_input_tmm.png", plot = p1, device = "png",
       "~/storage/Data/ms_evo_exwago/analysis/figures/",
       dpi = 300, width = 12, height = 18)
