setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/ip_cls_tabs/")
pacman::p_load(ggplot2, dplyr, tidyr, ggpubr, ggtext)
list.files()
tbFiles = list.files(pattern = ".*prop.txt")
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

df$exp = "Percentage of CPM from IP-enriched regions"
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


stak_prop = ggplot(ndf, aes(x=IP_list, y=Total_Raw_prop, fill= type)) + geom_bar(stat = "identity") +
  theme_light() + scale_fill_manual(values= col2, name=NULL) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_test() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5)) +
  xlab("") + ylab("% of IP CPM") +
  theme(legend.position = "bottom") +
  facet_wrap(~exp_list) +
  theme( strip.text = element_text(size = 16, colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) 

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

df = df[grepl("exons_As|DNA_As|LTR_As|LINE_As|Unknown|miRNA_S|rRNA_S|tRNA_S|lincRNA_As", df$Feature),]


df$Feature = factor(df$Feature, levels=c("exons_As", "DNA_As", "LTR_As", "LINE_As","Unknown", 
                                     "miRNA_S", "rRNA_S", "tRNA_S","lincRNA_As"))

col = c("#0B3BFA","#00C3DB", "#00AB14", "#90DB3F", "#D2F102", "#FFDA2C", "#F02765" , "#F75616", "#8846A3")
df$exp = "Relative enrichment by genomic category"



# Step 1: Set factor levels for ordering
df$Feature <- factor(df$Feature, levels = c(
  "exons_As", "DNA_As", "LTR_As", "LINE_As", "Unknown", 
  "miRNA_S", "rRNA_S", "tRNA_S", "lincRNA_As"
))

# Step 2: Create combined key and legend label
df$Feature_str <- interaction(df$Feature, df$str, sep = "_")
df$legend_label <- as.character(df$Feature)

# Step 3: Unique combinations, preserve Feature order
legend_df <- df %>%
  distinct(Feature_str, Feature, str, legend_label) %>%
  mutate(Feature = factor(Feature, levels = levels(df$Feature))) %>%
  arrange(Feature)

# Step 4: Build a named vector — colors named by Feature_str
colormap <- setNames(
  col[match(as.character(legend_df$Feature), levels(df$Feature))],
  legend_df$Feature_str
)

# Step 5: Shapes by strand
shapemap <- setNames(
  c(16, 17, 15)[as.numeric(legend_df$str)],
  legend_df$Feature_str
)



norms = ggplot(df, aes(x = IP, y = log2(norm), group = Feature_str)) +
  geom_line(color = adjustcolor("black", alpha.f = 0.5), linetype = "dashed") +
  geom_point(aes(color = Feature_str, shape = Feature_str), size = 4) +
  theme_test() +
  scale_color_manual(
    values = colormap,
    breaks = legend_df$Feature_str,
    labels = legend_df$legend_label,
    name = "Feature"
  ) +
  scale_shape_manual(
    values = shapemap,
    breaks = legend_df$Feature_str,
    labels = legend_df$legend_label,
    name = "Feature"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4)),
    shape = guide_legend(override.aes = list(size = 4))
  ) +
  facet_wrap(~exp) +
  labs(x = NULL, y = "log2(sRNA Density in IP)") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20, angle = 90, hjust = .5),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "inches")
  ) +
  theme(
    axis.text.x = element_text(
      size = 16,
      angle = 45,
      hjust = 1
    )
  ) +
  theme(
    axis.text.x = element_markdown(size = 16, angle = 30, hjust = 1)
  )

ggsave("nematode_ip_adult_class_dens_paper.pdf", device = "pdf", plot = norm,
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/", width = 14, height = 8, dpi = 300)

ggsave("nematode_ip_adult_class_proprs_paper.pdf", device = "pdf", plot = prop,
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/", width = 14, height = 7, dpi = 300)

ggsave("nematode_ip_adult_class_stacked_proprs_paper.pdf", device = "pdf", plot = stak_prop,
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/", width = 12, height = 7, dpi = 300)


ggsave("nematode_ip_adult_class_dens_simple.pdf", device = "pdf", plot = norms,
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/", width = 11, height = 7, dpi = 300)

##############################################################################################################################################################
setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/camera/class/")
camera_mir_files = list.files(pattern = ".*.txt")
camera_mir = lapply(camera_mir_files, read.delim)
names(camera_mir) =  c("A. ceylanicum", "PPW-1", "SAGO-1", "SAGO-2", "H. bakeri", "H. polygyrus", "N. brasiliensis", "T. circumcincta")
tmp.ids = camera_mir$`SAGO-1` %>% rownames()
c_lst = list()
for (i in names(camera_mir)) {
  tmp.tb = camera_mir[[i]]
  tmp.tb = na.omit(tmp.tb)
  if (length(setdiff(tmp.ids, rownames(tmp.tb))) >= 1) {
    missClass = setdiff(tmp.ids, rownames(tmp.tb))
    tmp.df = data.frame("NGenes"=0, "Direction"="Down", "PValue"=1,"FDR"=1)
    rownames(tmp.df) = missClass
    tmp.tb = rbind(tmp.tb,tmp.df)  
  }
  tmp.tb$FDR = ifelse(tmp.tb$FDR == 0 & tmp.tb$Direction == "Down", 1e-200, tmp.tb$FDR)
  tmp.tb$Trans_FDR = -log10(tmp.tb$FDR)
  tmp.tb$Trans_FDR = ifelse(tmp.tb$Direction == "Down", -1*(tmp.tb$Trans_FDR),  tmp.tb$Trans_FDR)
  tmp.tb$Feature = rownames(tmp.tb)
  tmp.tb$species = i
  #rownames(tmp.tb) = sub("other_ncRNA", "Other_ncRNA", rownames(tmp.tb))
  setdiff(as.character(df$Feature), rownames(tmp.tb))
  tmp.tb = tmp.tb[levels(df$Feature),]
  c_lst[[i]] = tmp.tb
}


camera_df = do.call(rbind, c_lst)

ndf = camera_df
#####################################################################################################################################################
ndf$Feature = factor(ndf$Feature, levels =  c("exons_As", "exons_S", "introns_As", "introns_S",
                                              "DNA_As", "DNA_S", "LTR_As", "LTR_S", "LINE_As", "LINE_S", "Unknown",
                                              "miRNA_S", "rRNA_S", "tRNA_S", "tRNA_As", "lincRNA_S", "lincRNA_As",
                                              "Other_repeat", "other_ncRNA" ,"intergenic")) 


ndf$species = factor(ndf$species, levels = c("SAGO-1", "SAGO-2", "PPW-1","A. ceylanicum", "T. circumcincta", "N. brasiliensis","H. bakeri", "H. polygyrus"))
ndf$sig = ifelse(ndf$FDR < .2  & ndf$Direction == "Up", ".", NA)
ndf$sig = ifelse(ndf$FDR < .1 & ndf$Direction == "Up", "*",ndf$sig)
ndf$sig = ifelse(ndf$FDR < .05 &  ndf$Direction == "Up", "**",ndf$sig)
ndf$sig = ifelse(ndf$FDR < .01 &  ndf$Direction == "Up", "***" ,ndf$sig)

ndf$log_fdr = -log10(ndf$FDR)
ndf$log_fdr = ifelse(ndf$Direction == "Down", ndf$log_fdr*-1, ndf$log_fdr)
campera_p = ggplot(ndf, aes(x = Feature , y = species, fill = log_fdr)) + theme_light() + geom_tile() +
  labs(title = "",  x = "",  y = "IP") + 
  scale_fill_gradient2(high = "#E52350", low="#00C3E5", mid="#FFFFFF", limits=c(-2,2), oob = scales::squish) +
  geom_text(aes(label = sig), vjust=.5) +
  theme(axis.title.y = element_text(face = "bold", size = 16),axis.text.y = element_text(face = "italic", size = 16)) +
  theme(axis.title.x = element_text(face = "bold", size = 16),axis.text.x = element_text( size = 16, angle = 90, hjust = 1, vjust = .5))


ggsave("nematode_ip_adult_class_camera_final.png", device = "png", plot = campera_p,
       path = "~/storage/Data/ms_exwago_evo/analysis/figures/", width = 10, height = 6, dpi = 300)



