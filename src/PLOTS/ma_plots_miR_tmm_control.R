setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/tmm/")
pacman::p_load(dplyr, ggplot2, rtracklayer, plyranges)
files = list.files(pattern = ".*deTbl.txt")
#files= files[!grepl("caeno", files)]
de_tbls = lapply(files, read.delim)
names(de_tbls) = c("A_ceylanicum", "PPW-1", "SAGO-1", "SAGO-2","H_bakeri","H_polygyrus", "N_brasiliensis",
                   "T_circumcincta")

n_lst = list()
up_n = list()
for (i in names(de_tbls)) {
  tmp.tbl = de_tbls[[i]]
  head(tmp.tbl)
  tmp.tbl$exp = paste(sub("_", "\\. ", i), "exWAGO vs adult total")
  sp = sub("_", "\\. ", i)
  tmp.tbl$exp = sp
  tmp.tbl$exp <- factor(tmp.tbl$exp)
  
  #tmp.tbl$exp = "T. circumcincta IP vs Adult total"
  #tmp.tbl$exp <- factor(tmp.tbl$exp)
  #levels(tmp.tbl$exp) <- c(expression(paste(italic("T. circumcincta "), "exWAGO IP vs Adult total")))
  tmp.tbl$col = ifelse(grepl("miRNA_S", rownames(tmp.tbl)), "#144675", "#C0C0C0")
  tmp.tbl$ord = ifelse(tmp.tbl$col == "#C0C0C0", 1, 2)
  tmp.tbl = tmp.tbl[order(tmp.tbl$ord),]
  tmp.tbl$sp = i
  tmp.tbl = tmp.tbl[,c("Ave_CPM", "logFC", "cl", "col", "ord", "exp","sp")]
  up_n[[i]] = tmp.tbl[tmp.tbl$cl == 1,] %>% nrow()
  n_lst[[i]] = tmp.tbl
}


lst2plot= list()
for (i in names(n_lst)) {
  tmp.tbl = n_lst[[i]]
  #tmp.tbl$exp = sub("T. circumcincta", sub("_", "\\. ", i), tmp.tbl$exp)
  tmp.tbl$exp = factor(tmp.tbl$exp)
  lst2plot[[i]] = tmp.tbl
}
df2plot = do.call(rbind, lst2plot)
df2plot$exp %>%  table()


df2plot$col = factor(df2plot$col, levels = c("#144675","#C0C0C0"))
df2plot$exp= factor(df2plot$exp, levels = c("H. bakeri", "H. polygyrus", "N. brasiliensis", "T. circumcincta", "A. ceylanicum", "SAGO-1", "SAGO-2", "PPW-1"))

p1 = ggplot(df2plot, aes(x=Ave_CPM, y=logFC, colour = col)) +
  geom_point() +  xlab(bquote(log[2]~' CPM')) + theme_test() +
  ylab(bquote(~log[2]~ 'FC IP/Total')) +
  theme(axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20),axis.text.y = element_text(size = 20)) +
  geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
  facet_wrap(~ exp, labeller = label_parsed) +
  scale_color_manual(values = adjustcolor(c("#144675","#C0C0C0"), alpha.f = .6), 
                     labels=factor(c('miRNA_S', 'Other'),
                                   levels = c('rRNA_S', "Other")), name="")+                                                                                               
  theme(legend.position = "bottom", legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=10))) +
  facet_wrap(~ exp, nrow = 2) +
  theme( strip.text = element_text(size = 18, face = "italic", colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) + 
  theme(panel.spacing = unit(1.3, "lines"))
p1

ggsave("nematode_exwago_vs_total_tmm_rRNA_all.png", plot = p1, device = "png",
       "~/storage/Data/ms_evo_exwago/analysis/figures/",
       dpi = 300, width = 16, height = 10)


