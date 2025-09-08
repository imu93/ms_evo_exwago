# This script aims to produce the enrichment plots I'll use in my PhD thesis
setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/")
pacman::p_load(edgeR, rtracklayer, dplyr)
################################################################################################################################################################################
# As input I'm gonna use dge objects for each DE analysis using rRNA_S to estimate norm factors
# My modified tables called deTblnorm.txt
# And my annotation
dgeFiles = list.files(path = "dge_tmm/", pattern = ".*.Rds")
deTabFiles = list.files(path="tmm/", pattern = ".*deTbl.txt")
dgeList = lapply(dgeFiles, function(x){readRDS(paste0("dge_tmm/",x))})
ann2PropFiles = list.files(path = "ann2prop/", pattern = ".*.Rds") 

names(dgeList) = c("ancylostoma_ceylanicum", "caenorhabditis_elegans_ppw1", "caenorhabditis_elegans_sago1", "caenorhabditis_elegans_sago2",
                   "heligmosomoides_bakeri", "heligmosomoides_polygyrus", "nippostrongylus_brasiliensis", "teladorsagia_circumcincta")

deTabList = lapply(deTabFiles, function(x){read.delim(paste0("tmm/",x))})
names(deTabList) = c("ancylostoma_ceylanicum", "caenorhabditis_elegans_ppw1", "caenorhabditis_elegans_sago1", "caenorhabditis_elegans_sago2",
                     "heligmosomoides_bakeri", "heligmosomoides_polygyrus", "nippostrongylus_brasiliensis","teladorsagia_circumcincta")

ann_str_lts = lapply(ann2PropFiles, function(x){readRDS(paste0("ann2prop/", x))}) 
names(ann_str_lts) = c("ancylostoma_ceylanicum", "caenorhabditis_elegans", 
                       "heligmosomoides_bakeri", "heligmosomoides_polygyrus", "nippostrongylus_brasiliensis","teladorsagia_circumcincta")
#####################################################################################################################################################
# I'll use a for to itarate over each set of obejects
for (exp in names(dgeList)) {
  print(exp)
  #exp = "caenorhabditis_elegans_ppw1"
  # I'll define relevant variables like the table of DE elements and dge object
  dge = dgeList[[exp]]
  de_tbl = deTabList[[exp]]
  sp = strsplit(exp,"_") %>% lapply(function(x){x[c(1:2)] %>% paste(collapse = "_")}) %>% unlist()
  ann_str = ann_str_lts[[sp]]
  cond = "other"
  cond = ifelse(sp == "ancylostoma_ceylanicum", "_adult", cond)
  cond = ifelse(sp == "caenorhabditis_elegans", "input", cond)
  cond = ifelse(sp == "nippostrongylus_brasiliensis", "_adult", cond)
  cond = ifelse(sp == "heligmosomoides_bakeri", "_total", cond)
  cond = ifelse(sp == "heligmosomoides_polygyrus", "_adult", cond)
  cond = ifelse(sp == "teladorsagia_circumcincta", "_adult", cond)
  ########################################################################################################################
  # I'm going to start with the raw % of reads in the IP enriched set 
  # For this I need IP columns from the table of counts
  ip = cpmByGroup(dge, normalized.lib.size=T, prior.count=0.01)
  head(ip)
  ip = ip[,grepl(cond, colnames(ip))] %>% as.data.frame()
  # But only those DE
  des = ip[match(rownames(de_tbl[de_tbl$cl != 1,]), rownames(ip)),]
  names(des) = rownames(de_tbl[de_tbl$cl != 1,])
  # I'll split by family to then create grups by class
  # if (grepl("caeno", sp) == T) {
  #    sum_fam = split(des, gsub(":.*", "", names(des)) %>% sub("(_S|_As).*", "\\1", .) %>% sub("_rnd.*", "",.) %>% sub("(_CE|_Ce).*", "", .) ) 
  # }else{
  sum_fam = split(des, gsub(":.*", "", names(des)))
  
  #}
  
  # So here we have the % of mean expression among replicates of each family:
  fams = unlist(lapply(sum_fam, function(x){(sum(x)/sum(des))*100}))
  
  # Make sure to add $ to relevant categories to do not count twice 
  categories = c("DNA.*As|^RC.*As", "DNA.*S$|^RC.*S$", "^LINE.*As$|^PLE.*As|^Retrop.*As", "^LINE.*S$|^PLE.*S|Retrop.*S", "LTR.*As", "LTR.*S$",
                 "Unk", "Low|Simple|Sat|SINE|MITE", "lincRNA_As",  "miRNA_S", 
                 "*s_rRNA_S|LSU.*_S|^rRNA_S", "^tRNA_S", 
                 "yRNA|snRNA|snoRNA|other_ncRNA|piRNA|_rRNA_As|^rRNA_As|tRNA_As|lncRNA|miRNA_As|LSU.*_As|lincRNA_S|circRNA", "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  
  names(categories) = c("DNA_As", "DNA_S", "LINE_As", "LINE_S", "LTR_As", "LTR_S", "Unknown", 
                        "Other_repeat", "lincRNA_As",  "miRNA_S",
                        "rRNA_S", "tRNA_S", 
                        "other_ncRNA", "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  ############################################################################################################################
  # now I need to sum the percentage of CPM by category
  allFams = list()
  for (i in names(categories)) {
    #i="other_ncRNA"
    print(i)
    tmp.cat = categories[i] # extract the categoy'
    #ids = fams[grepl(tmp.cat, names(fams))] %>% names() 
    tmp.sqs = fams[grepl(tmp.cat, names(fams))] %>%  sum() # just sum
    allFams[[i]] =tmp.sqs 
  }
  #names(allFams)= NULL
  #unlist(allFams)[unlist(allFams) %>% duplicated()]
  
  
  # now create a df 
  df_cls = unlist(allFams) %>% as.data.frame()
  df_cls
  colnames(df_cls) = c("Raw_prop")
  df_cls$Feature = rownames(df_cls)
  df_cls$Fam = "Class"
  df_cls$Feature = factor(df_cls$Feature, levels = rev(df_cls$Feature))
  # Just and small control
  sum(df_cls$Raw_prop) %>% print()
  ###########################################################################################################################
  # Now I need to know the number of bases represented by elements within this categories 
  
  categories = c("DNA.*As|^RC.*As", "DNA.*S$|^RC.*S$", "^LINE.*As$|^PLE.*As|^Retrop.*As", "^LINE.*S$|^PLE.*S|Retrop.*S", "LTR.*As", "LTR.*S$",
                 "Unk", "Low|Simple|Sat|SINE|MITE", "lincRNA_As",  "miRNA_S", 
                 "*s_rRNA_S|LSU.*_S|^rRNA_S", "^tRNA_S", 
                 "yRNA|snRNA|snoRNA|other_ncRNA|piRNA|_rRNA_As|^rRNA_As|tRNA_As|lncRNA|miRNA_As|LSU.*_As|lincRNA_S|circRNA", "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  
  names(categories) = c("DNA_As", "DNA_S", "LINE_As", "LINE_S", "LTR_As", "LTR_S", "Unknown", 
                        "Other_repeat", "lincRNA_As",  "miRNA_S",
                        "rRNA_S", "tRNA_S", 
                        "other_ncRNA", "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  ###########################################################################################################################
  # I use the next for to get the number of bases represented by each class
  wdt_lst = list()
  for (i in names(categories)) {
    tmp.cat = categories[i]
    tmp.wdt = ann_str[grepl(tmp.cat, ann_str$type),] %>% width() %>% sum() # Chenge this line to specificaly ann_str to estimate the % of expressed bases
    wdt_lst[[i]] = tmp.wdt
  }
  # an add to the df
  df_cls$bases = unlist(wdt_lst)
  df_cls$pro_bases = 100*(df_cls$bases/sum(df_cls$bases))
  outFile2 = paste0("~/storage/Data/ms_evo_exwago/analysis/de_tables/tmm_noIP_tabs/", exp,".cls_prop_ub_tmm.txt")
  write.table(df_cls, outFile2, quote = F, sep = "\t")
  ##############################################################################################################################
}  
