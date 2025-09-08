setwd("~/storage/Data/ms_exwago_evo/analysis/repeat_comparison/")

pacman::p_load(dplyr, ggplot2, stringr)
cur_path = "curated_ann/"
wbp_path = "wbp_ann/"

cur_files = list.files(pattern = ".*.txt", path = cur_path)
wbp_files = list.files(pattern = ".*.txt", path = cur_path)
cur_tab =  lapply(cur_files, function(x){read.table(paste0(cur_path,x), header = T)})
names(cur_tab) = sub(".class_repeats.txt", "", cur_files)

wbp_tab =  lapply(wbp_files, function(x){read.table(paste0(wbp_path,x), header = T)})
names(wbp_tab) = sub(".class_repeats.txt", "", wbp_files)


sim_tabs = function(x){
  rownames(x) = x$cat
  x = x[,2:length(x)]
  DNA = x["DNA",] +  x["RC",]
  Other_repeat = x["Other_repeat",] + x["Satellite",]
  x = x %>% 
    filter(!(rownames(.) %in%
               c("DNA", "RC", "Satellite", "Other_repeat", "Total"))) 
  x = rbind(DNA, x, Other_repeat)
  x$GS = x["LINE", "GS"]
  x = x[c("DNA", "LINE", "SINE", "LTR", "Unknown"),]
  x$class = rownames(x)
  return(x)
}
cur_tab = lapply(cur_tab, sim_tabs)
wbp_tab = lapply(wbp_tab, sim_tabs)


reps_cur = list()
for (i in names(cur_tab)) {
  tmp.tab = cur_tab[[i]]
  tmp.tab$library = gsub(".*\\.(.*)_.*", "\\1", i)
  tmp.tab$library = sub(".summary", "", tmp.tab$library)
  tmp.tab$library =  str_to_title(tmp.tab$library) %>% 
    sub("(.).*_", paste0("\\1", ". "), .)
  tmp.tab$type = "Curated"
  reps_cur[[i]] = tmp.tab
}

reps_wbp = list()
for (i in names(wbp_tab)) {
  tmp.tab = wbp_tab[[i]]
  tmp.tab$library = gsub(".*\\.(.*)_.*", "\\1", i)
  tmp.tab$library = sub(".summary", "", tmp.tab$library)
  tmp.tab$library =  str_to_title(tmp.tab$library) %>% 
    sub("(.).*_", paste0("\\1", ". "), .)
  tmp.tab$type = "Uncurated"
  reps_wbp[[i]] = tmp.tab
}

names(reps_cur) = NULL
names(reps_wbp) = NULL
wbp = do.call(rbind,reps_wbp)
cur = do.call(rbind, reps_cur)
colnames(wbp) = paste0("wbp_", colnames(wbp))
colnames(cur) = paste0("cur_", colnames(cur))

df = cbind(wbp, cur)
df$species = factor(df$wbp_library, levels = c("H. bakeri", 
                                               "N. brasiliensis", 
                                               "T. circumcincta",
                                               "A. ceylanicum"))

df$class = factor(df$wbp_class, levels = c("DNA","LTR", "LINE", 
                                           "Unknown", "SINE"))

colors = c("#0092EFFF", "#00860F", "#90DB3F", "#D2F102", "#46026E")

df %>% ggplot(aes(wbp_GP, cur_GP, color = class)) + geom_point(size=4) + 
  geom_abline(linetype ="dashed") + facet_wrap(~species) + theme_test() +
  theme( strip.text = element_text(size = 13, face = "bold.italic", colour = "black")) +
  xlab("Genomic % using uncurated library") +
  ylab("Genomic % using curated library") +
  theme(axis.title.x = element_text(face = "bold", size = 16), 
        axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 16)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 16)) +
  scale_color_manual(values=colors, name= NULL) 



ggsave("Fig1C.pdf", device = "pdf",  plot = last_plot(), 
       width = 8, height = 6.7, dpi=300,
       path = "~/storage/Data/ms_exwago_evo/analysis/figures/")  
