setwd("~/storage/Data/ms_evo_exwago/analysis/orthofinder/ltr_ms/")

pacman::p_load(ggtree, dplyr, ggplot2, stringr, ape, treeio, ape, ggplot2, pals, ggpubr, stringr)

tree = read.tree("nematoda_phylo_sort.txt")

order = ifelse(grepl("pri|mic|mes|oscheius|hete|cae", tree$tip.label),
                "Rhabditida", "Strongylida")

spe = sub(".*_", "", tree$tip.label)
genus = str_to_title(tree$tip.label) %>% str_replace_all("[^A-Z]+", ".")
n_tip_name = paste(genus, spe)
names(order) = n_tip_name
tree$tip.label = n_tip_name

td = data.frame(label = names(order), trait = order)
internal_nodes <- data.frame(node = 28:53,
                             trait = ifelse(28:53 <= 40, "Rhabditida", "Strongylida"))

p <- ggtree(tree, ladderize = FALSE, aes(color = trait), ) %<+% td
#p$data$label[!p$data$isTip] <- as.character(p$data$node[!p$data$isTip])
p$data <- left_join(p$data, internal_nodes, by = "node", suffix = c("", ".y"))
p$data$trait <- ifelse(is.na(p$data$trait), p$data$trait.y, p$data$trait)

nodes_to_label <- with(p$data, !isTip & grepl("^[0-9]+$", label) & as.numeric(label) < 100)

trp <- p +
  geom_tiplab(color = "black", align = TRUE, fontface = "italic", ) +
  scale_color_manual(values = c("Rhabditida" = "#0298F2", "Strongylida" = "#E0021B")) +
  geom_treescale(x = 0, y = 20) +
  coord_cartesian(clip = 'off') +
  geom_label2(data = p$data[nodes_to_label, ], aes(label = label)) +
  theme(legend.position = "none") +
  xlim(0, 1.5)

setwd("~/storage/Data/ms_evo_exwago/analysis/orthofinder/repeat_class_props/")

cpunts = list.files(pattern = ".*.class_repeats.txt")
tbls = lapply(cpunts, function(x){read.table(x,header = T)})
names(tbls) = sub("\\..*", "", cpunts)
tree = read.tree("~/storage/Data/ms_evo_exwago/analysis/orthofinder/ltr_ms/nematoda_phylo_sort.txt")
tbls  = tbls[rev(tree$tip.label)]

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


unks = unlist(unk_list)
summary(unks)


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
df$species = factor(df$species, levels = rev(unique(df$species)))
df$cat = factor(df$cat, levels = c("non-repeat","Other_repeat","SINE","Unknown", "LINE","LTR", "DNA"))

colors = rev(c(c("#0092EFFF","#028610","#8FDC3F","#D2F202","#46006E"), "#575757", "#C7C7C7"))

p1 = ggplot(df, aes(x= Mb, y= species, fill = cat)) + geom_bar(stat = "identity") + theme_light() +
  theme(axis.text.y = element_text(angle = 90, hjust = .5)) +
  scale_fill_manual(values = colors)  +  ylab("") + xlab("Mb") +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(legend.position = "none") + theme(legend.title=element_blank())

p2 = ggplot(df, aes(x= GP, y= species, fill = cat)) + geom_bar(stat = "identity") + theme_light() +
  theme(axis.text.y = element_text(angle = 90, hjust = .5)) +
  scale_fill_manual(values = colors) +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  ylab("") + xlab("Genomic %")  +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + theme(legend.title=element_blank()) +
  theme(legend.position = "none") 

#saveRDS(p1, "~/storage/Data/genomes_wbp18/phylo/species_mb.RDS")
#saveRDS(p2, "~/storage/Data/genomes_wbp18/phylo/species_gp.RDS")

phylo1 = ggarrange(trp, p1,p2, nrow=1, align = "h")


#dir.create("../../figures")
ggsave("nematode_repeat_span.png", device = "png", path = "~/storage/Data/ms_evo_exwago/analysis/figures/", 
       dpi=300, width = 12, height = 8, plot= phylo1)

  
