pacman::p_load(ggplot2)
setwd("~/storage/Data/ms_evo_exwago/analysis/fn_cls/raw/")
tabs = list.files(pattern = "cpm_adult.txt")
tabs = lapply(tabs, read.delim)
names(tabs) = c("A_ceylanicum", "H_bakeri", "H_polygyrus", "N_brasiliensis", "T_circumcincta")
n_tabs = list()
for (i in names(tabs)) {
  tmp.tab= tabs[[i]]
  tmp.tab$sp = sub("_", ". ", i)
  n_tabs[[i]] = tmp.tab
}
names(n_tabs) = NULL
df2plot = do.call(rbind, n_tabs)

df2plot$Class = factor(df2plot$Class,
                       levels =rev(c("exons", "introns", "DNA", "LTR",
                                     "LINE", "Unknown","miRNA", "rRNA", 
                                     "tRNA", "lincRNA", "Other_ncRNA", 
                                     "Other_repeat", "Unannotated")))



colors = rev(c("#0B3BFA","#007CFF", "#00C3DB","#00AB14", "#90DB3F", "#D2F102",
               "#FFDA2C", "#F02765", "#F75616" ,"#8846A3","#601573", "#807B77", "#DBDBDB"))

df2plot$sp = factor(df2plot$sp, levels = c("H. bakeri", "H. polygyrus", "N. brasiliensis", "T. circumcincta", "A. ceylanicum"))
df2plot$Length = factor(df2plot$Length)


ggplot(df2plot, aes(x=Length, y=Reads, fill=Class)) + 
  geom_bar(stat = "identity") + theme_test() +
  scale_fill_manual(values = colors, guide = guide_legend(reverse = TRUE), name="") + facet_wrap(~sp, nrow=1) +
  ylab("CPM of aligned reads") + xlab("Length distribution (nt)") + 
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20, angle = 90, hjust = .5)) + 
  theme( strip.text = element_text(size = 20, colour = "black", face =  "italic")) +
  theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  scale_x_discrete(breaks = seq(18, 27, by = 3)) +
  theme(panel.spacing = unit(1.3, "lines"))

ggsave("strongylis_fn_cls_simple_adult.pdf", device = "pdf", width = 16, height = 6, dpi = 300,
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/")    
