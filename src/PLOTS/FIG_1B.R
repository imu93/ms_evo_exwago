setwd("~/backup/mac/Data/ms_exwago_evo/analysis/repeat_comparison/anns2plot/")

pacman::p_load(ape, ggplot2, pals, ggpubr, stringr, ggtree, dplyr)
cpunts = list.files(pattern = ".*.class_repeats.txt")
tbls = lapply(cpunts, function(x){read.table(x,header = T)})
names(tbls) = sub("\\..*", "", cpunts)
#tree = read.tree("~/storage/Data/manuscript_2024/analysis/phylo/nematoda_phylo_sort.txt")
#tbls  = tbls[rev(tree$tip.label)]

tb2df = list()
unk_list = list() 
for (i in names(tbls)) {
  print(i)
  
  x= tbls[[i]]
  NR = unique(x$GS)- x[x$cat=="Total","Mb"]
  tmp.gp = (NR*100)/unique(x$GS)
  x = x[!x$cat == "Total",]
  nr = data.frame("cat"="non-repeat","Mb"=NR,"Cp"=NA,"RP"=NA,
                  "GP"=tmp.gp, "GS"=NA)
  x = rbind(x,nr)
  rownames(x) = x$cat
  x$cat= NULL
  x$sp = NULL
  DNA = colSums(rbind(x["DNA",], x["RC",]))
  Other_repet = colSums(rbind(x["Other_repeat",], x["Satellite",]))
  tmp.x =  x[!grepl("DNA|RC|Other_rep|Satellite", rownames(x)),]
  x = rbind(tmp.x, "DNA"=DNA, "Other_repeat"=Other_repet)
  x = x[c("non-repeat", "DNA","SINE","LINE","LTR", "Unknown","Other_repeat"),]
  
  x$cat = rownames(x)
  x$cat = factor(x$cat, levels = c("non-repeat", "DNA","SINE","LINE",
                                   "LTR", "Unknown", 
                                   "Other_repeat"))
  
  x$sp = factor(i)
  unk_list[[i]] = x[,"GP"][6]
  tb2df[[i]] = x
}
lapply(tb2df, function(x){x[,"GP"] %>% sum()})

names(tb2df) = NULL
df = do.call(rbind, tb2df)
df$species = as.character(df$sp)
abbreviations = str_to_title(df$species) %>% str_extract_all("\\b\\w") %>% 
  unlist()

# Extract the species name after underscore
species = sub(".*_", "", df$species)

# Combine the abbreviations and species
output_string = paste(abbreviations, species, sep = ". ")
df$species = output_string
df$species = factor(df$species, levels = rev(c("T. circumcincta", "N. brasiliensis", "H. polygyrus", "H. bakeri", "A. ceylanicum", "C. elegans")))
df$cat = factor(df$cat, levels = c("non-repeat","Other_repeat","SINE","Unknown", "LINE","LTR", "DNA"))

colors = rev(c(c("#0092EFFF","#026210","#61C93A","#D2F202","#46006E"), "#575757", "#DBDBDB"))

p1 = ggplot(df, aes(x=Mb, y= species, fill = cat)) + geom_bar(stat = "identity") + theme_test() +
  theme(axis.text.y = element_text(angle = 90, hjust = .5)) +
  scale_fill_manual(values = colors)  +  xlab("Genome size (Mb)") + ylab("") +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +  theme(legend.position = "none") + theme(legend.title=element_blank())
p1

p2 = ggplot(df, aes(x= GP, y= species, fill = cat)) + geom_bar(stat = "identity") + theme_test() +
  theme(axis.text.y = element_text(angle = 90, hjust = .5)) +
  scale_fill_manual(values = colors) +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  xlab("Genomic %") + ylab("")  +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + theme(legend.title=element_blank()) +
  theme(legend.position = "none")
p2


tree = read.tree("~/backup/mac/Data/manuscript_2024/analysis/phylo/nematoda_phylo_sort.txt")
tree = drop.tip.phylo(phy = tree, tip = tree$tip.label[!grepl("helig|ceyla|circ|brasi|eleg", tree$tip.label)])
plot(tree)
nodelabels()
tr = ggtree(tree, ladderize = FALSE) +
  geom_tiplab(color="black", align = T, font.face = "italic") +
  geom_treescale(x=0, y=4) +
  coord_cartesian(clip = 'off') +
  xlim(0, 1)

ggarrange(tr, p1, p2, ncol=3, align = "h", widths = c(0.8,.8,.8))
ggsave("nematode_repeat_span_3.pdf", device = "pdf", dpi = 300,
       width = 10, height = 4, plot = last_plot(),
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/")
