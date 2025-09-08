setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/rRNA_noIP_tabs/")
pacman::p_load(ggplot2, dplyr, tidyr, ggpubr)
list.files()
tbFiles = list.files(pattern = ".*.txt")
tbList = lapply(tbFiles, read.delim)
names(tbList) = c("A. ceylanicum", "PPW-1", "SAGO-1", "SAGO-2", "H. bakeri", "H. polygyrus", "N. brasiliensis", "T. circumcincta")
##############################################################################################################################################################
lst2plot = list()
for (ip in names(tbList)) {
  tmp.ip = tbList[[ip]]
  tmp.ip$norm = (tmp.ip$Raw_prop+1)/(tmp.ip$pro_bases+1)
  tmp.ip$IP = ip
  lst2plot[[ip]] = tmp.ip
}
##############################################################################################################################################################
df = do.call(rbind, lst2plot)
df$IP = factor(df$IP, levels = c("H. bakeri", "H. polygyrus", "N. brasiliensis", "T. circumcincta", "A. ceylanicum", "PPW-1", "SAGO-1", "SAGO-2"))
df$Feature = factor(df$Feature, levels = c("exons_As", "exons_S", "introns_As", "introns_S",
                                           "DNA_As", "DNA_S","LTR_As", "LTR_S", "LINE_As", 
                                           "LINE_S", "Unknown","miRNA_S", "rRNA_S", "tRNA_S", "lincRNA_As",
                                           "other_ncRNA", "Other_repeat", "intergenic") %>% rev)
df

df$str = ifelse(grepl("_As", df$Feature), 1, ifelse(grepl("_S", df$Feature), 2, 3))
df$str = factor(df$str)

##############################################################################################################################################################
col =  c("#0B3BFA","#0B3BFA","#007CFF","#007CFF", "#00C3DB","#00C3DB",
         "#00AB14",  "#00AB14", "#90DB3F","#90DB3F", "#D2F102", "#FFDA2C",
         "#F02765", "#F75616","#8846A3","#601573", "#807B77", "#DBDBDB") %>% rev()

df$exp = "Percentage of counts from IP-enriched regions"
prop = ggplot(df, aes(x=IP, y=Raw_prop, group = Feature, color = Feature, shape = str)) + geom_line(colour=adjustcolor("black", alpha.f = .7), linetype="dashed")  + geom_point(size=5) +
  theme_light() + scale_color_manual(values= col) +
  geom_point(aes(color=Feature, shape = str), size = 5) + theme_test() +
  scale_color_manual(values= col) +
  theme(axis.title.x = element_text(face = "bold", size = 16), axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 16, angle = 90, hjust = .5)) +
  xlab("") + ylab("% of IP CPM") +
  theme(legend.position = "bottom") +
  facet_wrap(~exp) +
  theme( strip.text = element_text(size = 18, face = "bold", colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) 

df$type = sub("_As|_S", "", df$Feature)  
lst = split(df, df$IP)

collapsed = lapply(lst, function(x){
  tmp.df = group_by(x,type) %>%
    summarise(Total_Raw_prop = sum(Raw_prop, na.rm = TRUE), 
              IP_list = toString(unique(IP)),
              exp_list = toString(unique(exp))) %>% 
    as.data.frame()
  
})
names(collapsed) = NULL
ndf = do.call(rbind, collapsed)

col2 =  c("#0B3BFA","#007CFF", "#00C3DB",
          "#00AB14",   "#90DB3F", "#D2F102", "#FFDA2C",
          "#F02765", "#F75616","#8846A3","#601573", "#807B77", "#DBDBDB") %>% rev()

ndf$IP_list = factor(ndf$IP_list, levels = c("H. bakeri", "H. polygyrus", "N. brasiliensis", "T. circumcincta", "A. ceylanicum", "PPW-1", "SAGO-1", "SAGO-2"))
ndf$type = factor(ndf$type, levels = c("exons", "introns",
                                       "DNA", "LTR", "LINE", "Unknown","miRNA", "rRNA", "tRNA", "lincRNA",
                                       "other_ncRNA", "Other_repeat", "intergenic") %>% rev)

ndf$exp_list = "Percentage of non-enriched exWAGO sRNAs"
stak_prop = ggplot(ndf, aes(x=IP_list, y=Total_Raw_prop, fill= type)) + geom_bar(stat = "identity") +
  theme_light() + scale_fill_manual(values= col2, name=NULL) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_test() +
  theme(axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20, angle = 90, hjust = .5)) +
  xlab("") + ylab("% of IP CPM") +
  theme(legend.position = "bottom") +
  facet_wrap(~exp_list) +
  theme( strip.text = element_text(size = 20,colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) 


ggsave("non_enriched_rRNA_exWAGO.pdf", device = "pdf", plot = stak_prop,
       width = 12, height = 7, dpi = 300, 
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/")


hb_df = ndf[ndf$IP_list == "H. bakeri",]
hb_stak_prop = ggplot(hb_df, aes(x=IP_list, y=Total_Raw_prop, fill= type)) + geom_bar(stat = "identity") +
  theme_light() + scale_fill_manual(values= col2, name=NULL) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_test() +
  theme(axis.title.x = element_text(face = "bold", size = 18), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5)) +
  xlab("") + ylab("% of IP CPM") +
  theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) + coord_flip()


ggsave("nxHelBake1.1_non_enriched_miR_exWAGO.pdf", device = "pdf", plot = hb_stak_prop,
       width = 12, height = 6, dpi = 300, 
       path = "~/storage/Data/ms_evo_exwago/analyses/figures/")







norm = ggplot(df, aes(x=IP, y=norm %>% log2(),  group = Feature, color = Feature, shape = str )) + geom_line(colour=adjustcolor("black", alpha.f = .7), linetype="dashed") + 
  geom_point(aes(color=Feature, shape = str), size = 5) + theme_test() +
  scale_color_manual(values= col) +
  theme(axis.title.x = element_text(face = "bold", size = 16), axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 16, angle = 90, hjust = .5)) +
  xlab("") + ylab("log2(% of IP CPM + 1 / % of genomic bases + 1)") +
  theme(legend.position = "bottom") +
  theme( strip.text = element_text(size = 12.5, face = "bold", colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) 

df = df[grepl("exons_As|DNA_As|LINE_As|LTR_As|Unknown|miRNA_S|rRNA_S|tRNA_S|lincRNA_As", df$Feature),]


df$Feature = factor(df$Feature, levels=c("exons_As", "DNA_As", "LINE_As", "LTR_As", "Unknown", 
                                         "miRNA_S", "rRNA_S", "tRNA_S","lincRNA_As"))

col = c("#0B3BFA","#00C3DB", "#00AB14", "#90DB3F", "#D2F102", "#FFDA2C", "#F02765" , "#F75616", "#8846A3")
df$exp = "Relative enrichment by genomic category"
norms = ggplot(df, aes(x=IP, y=norm %>% log2(),  group = Feature, color = Feature, shape = str )) + geom_line(colour=adjustcolor("black", alpha.f = .7), linetype="dashed") + 
  geom_point(aes(color=Feature, shape = str), size = 5) + theme_test() +
  scale_color_manual(values= col) +
  theme(axis.title.x = element_text(face = "bold", size = 16), axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5)) +
  xlab("") + ylab("log2(sRNA Density in IP)") +
  theme(legend.position = "bottom") +
  facet_wrap(~exp) +
  theme( strip.text = element_text(size = 16, face = "bold", colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) 



ggplot(df, aes(x=IP, y = Raw_prop, group = Feature, fill = Feature)) +
  geom_bar(stat = "identity") +
  theme_light() + scale_fill_manual(values= col, name=NULL) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_test() +
  theme(axis.title.x = element_text(face = "bold", size = 18), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5)) +
  xlab("") + ylab("% of IP CPM") +
  theme(legend.position = "bottom") +
  facet_wrap(~exp_list) +
  theme( strip.text = element_text(size = 18, face = "bold", colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) 


