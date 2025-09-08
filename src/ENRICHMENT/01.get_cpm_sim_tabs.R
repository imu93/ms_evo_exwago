# This script aims to produce the enrichment plots I'll use in my PhD thesis
setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/")
pacman::p_load(edgeR, rtracklayer, dplyr)
################################################################################################################################################################################
# As input I'm gonna use dge objects for each DE analysis using rRNA_S to estimate norm factors
# My modified tables called deTblnorm.txt
# And my annotation
dgeFiles = list.files(path = "dge/", pattern = ".*.Rds")
deTabFiles = list.files(path="rrna_norm/", pattern = ".*deTbl.txt")
dgeList = lapply(dgeFiles, function(x){readRDS(paste0("dge/",x))})
ann2PropFiles = list.files(path = "ann2prop/", pattern = ".*.Rds") 

names(dgeList) = c("ancylostoma_ceylanicum", "caenorhabditis_elegans_ppw1", "caenorhabditis_elegans_sago1", "caenorhabditis_elegans_sago2",
                   "heligmosomoides_bakeri", "heligmosomoides_polygyrus", "nippostrongylus_brasiliensis", "teladorsagia_circumcincta")

deTabList = lapply(deTabFiles, function(x){read.delim(paste0("rrna_norm//",x))})
names(deTabList) = c("ancylostoma_ceylanicum", "caenorhabditis_elegans_ppw1", "caenorhabditis_elegans_sago1", "caenorhabditis_elegans_sago2",
                     "heligmosomoides_bakeri", "heligmosomoides_polygyrus", "nippostrongylus_brasiliensis","teladorsagia_circumcincta")

ann_str_lts = lapply(ann2PropFiles, function(x){readRDS(paste0("ann2prop/", x))}) 
names(ann_str_lts) = c("ancylostoma_ceylanicum", "caenorhabditis_elegans", 
                       "heligmosomoides_bakeri", "heligmosomoides_polygyrus", "nippostrongylus_brasiliensis","teladorsagia_circumcincta")
#####################################################################################################################################################
# I'll use a for to itarate over each set of obejects
for (exp in names(dgeList)) {
  print(exp)
  #exp = "caenorhabditis_elegans_sago1"
  # I'll define relevant variables like the table of DE elements and dge object
  dge = dgeList[[exp]]
  de_tbl = deTabList[[exp]]
  sp = strsplit(exp,"_") %>% lapply(function(x){x[c(1:2)] %>% paste(collapse = "_")}) %>% unlist()
  ann_str = ann_str_lts[[sp]]
  cond = "other"
  cond = ifelse(sp == "ancylostoma_ceylanicum", "IP", cond)
  cond = ifelse(sp == "caenorhabditis_elegans", "Ip", cond)
  cond = ifelse(sp == "nippostrongylus_brasiliensis", "IP", cond)
  cond = ifelse(sp == "heligmosomoides_bakeri", "IP", cond)
  cond = ifelse(sp == "heligmosomoides_polygyrus", "IP", cond)
  cond = ifelse(sp == "teladorsagia_circumcincta", "IP", cond)
  ########################################################################################################################
  # I'm going to start with the raw % of reads in the IP enriched set 
  # For this I need IP columns from the table of counts
  
  ip = cpmByGroup(dge,  normalized.lib.size=T)
  head(ip)
  ip = ip[,grepl(cond, colnames(ip))] %>% as.data.frame()
  # But only those DE
  des = ip[match(rownames(de_tbl[de_tbl$cl == 1,]), rownames(ip)),]
  names(des) = rownames(de_tbl[de_tbl$cl == 1,])
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
                 "yRNA|snRNA|snoRNA|other_ncRNA|piRNA|rRNA_As|tRNA_As|lncRNA|miRNA_As|LSU.*_As|lincRNA_S|circRNA", "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  
  names(categories) = c("DNA_As", "DNA_S", "LINE_As", "LINE_S", "LTR_As", "LTR_S", "Unknown", 
                        "Other_repeat", "lincRNA_As",  "miRNA_S",
                        "rRNA_S", "tRNA_S", 
                        "other_ncRNA", "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  ############################################################################################################################
  # now I need to sum the percentage of CPM by category
  allFams = list()
  #id_lst = list()
  for (i in names(categories)) {
    #i="other_ncRNA"
    print(i)
    tmp.cat = categories[i] # extract the categoy'
    #id_lst[[i]] = fams[grepl(tmp.cat, names(fams))] %>% names() 
    tmp.sqs = fams[grepl(tmp.cat, names(fams))] %>%  sum() # just sum
    allFams[[i]] =tmp.sqs 
  }
  #names(id_lst)= NULL
  #setdiff(fams %>% names(), id_lst)
  
  
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
                 "yRNA|snRNA|snoRNA|other_ncRNA|piRNA|rRNA_As|tRNA_As|lncRNA|miRNA_As|LSU.*_As|lincRNA_S|circRNA", "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  
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
  outFile2 = paste0("~/storage/Data/ms_evo_exwago/analysis/de_tables/ip_cls_tabs/", exp,".cls_prop.txt")
  write.table(df_cls, outFile2, quote = F, sep = "\t")
  ##############################################################################################################################
  allFams = ann_str
  unstr = allFams[strand(allFams) == "*",]
  unstr = table(sub(":.*", "", unstr$ID)) 
  allFams=allFams[grepl("DNA|LINE|SINE|LTR|PLE|RC|Retro", allFams$ID),] 
  
  #if (grepl("caeno", sp) == T) {
  #  allFams_genome = split(width(allFams), sub("(_S|_As).*", "\\1", allFams$ID) %>%  sub("(:.*)", "", .) %>% sub("_CE.*|_Ce.*|_rnd.*", "", .))
  #}else{
  allFams_genome = split(width(allFams), sub("_.*(_.*)", "\\1", allFams$ID) %>%  sub("(:.*)", "", .) %>% sub("(_CE|_Ce).*", "", .))
  #}
  
  allFams_genome_wdt = lapply(allFams_genome, function(x){(sum(x)/2)/1e6}) %>% unlist() 
  allFams_genome_wdt_fx = ifelse(names(allFams_genome_wdt) %in% names(unstr), allFams_genome_wdt*2, allFams_genome_wdt)
  names(allFams_genome_wdt_fx) = names(allFams_genome_wdt)
  
  reps_wdt_mb = allFams_genome_wdt_fx[!grepl("exons|introns|miRNA|other|inter|^tRNA|s_rRNA|miRNA|piRNA|yRNA|snRNA|snoRNA|other_ncRNA|Low|Sat|Simple|inter|lincRNA|lncRNA", 
                                             names(allFams_genome_wdt_fx))]
  
  
  # I'll split by family to then create grups by class
  #if (grepl("caeno", sp) == T) {
  #  sum_fam = split(des, sub("(_S|_As).*", "\\1", names(des)) %>%  sub("(:.*)", "", .) %>% sub("_CE.*|_Ce.*|_rnd.*", "", .))
  #}else{
  sum_fam = split(des,  sub("_.*(_.*)", "\\1", names(des)) %>% sub(":.*", "",.) %>% sub("(_CE|_Ce).*", "", .) ) 
  #}
  names(sum_fam)
  
  
  # So here we have the % of mean expression among replicates of each family:
  fams = unlist(lapply(sum_fam, function(x){(sum(x)/sum(des))*100}))
  
  
  fams = fams[!grepl("8s|exons|introns|miRNA|other|inter|^tRNA|s_rRNA|miRNA|piRNA|yRNA|snRNA|snoRNA|other_ncRNA|Low|Sat|Simple|inter|lincRNA|lncRNA|Unk",
                     names(fams))]
  
  no_fams = setdiff(names(reps_wdt_mb), names(fams))
  vls = rep(0, length(no_fams))
  names(vls) = no_fams
  fams = c(fams, vls)
  fams = fams[names(reps_wdt_mb)]
  reps_wdt_mb =  reps_wdt_mb[match(names(fams), names(reps_wdt_mb))]
  reps_wdt_mb = reps_wdt_mb[names(fams)]
  fam_gprop = 100*(reps_wdt_mb/sum(reps_wdt_mb))
  fam_gprop = fam_gprop[names(fams)]
  df_fam = data.frame("raw_IP_prop"=fams, reps_wdt_mb, fam_gprop)
  
  # Let's estimate CPM by class
  cpms =  cpm(dge,normalized.lib.sizes = T, prior.count = 1)
  cpms = cpms[,grepl(cond, colnames(cpms))]
  cpms = rowMeans(cpms)
  #if (grepl("caeno", sp) == T) {
  #  fam_ave_cpm = split(cpms, sub("(_S|_As).*", "\\1", names(cpms)) %>%  sub("(:.*)", "", .) %>% sub("_CE.*|_Ce.*|_rnd.*", "", .))
  #}else{
  fam_ave_cpm = split(cpms, sub("_.*(_.*)", "\\1", names(cpms)) %>%  sub("(:.*)", "", .) %>% sub("(_CE|_Ce).*", "", .)) 
  #}
  
  
  fam_ave_cpm = fam_ave_cpm[rownames(df_fam)]
  fam_ave_cpm =  lapply(fam_ave_cpm, sum) %>% unlist()
  df_fam$Ave_CPM = fam_ave_cpm
  df_fam = df_fam[!grepl("SINE|MITE",rownames(df_fam)),]
  outFile2 = paste0("~/storage/Data/ms_evo_exwago/analysis/de_tables/ip_fam_tabs/", exp,".fam_prop.txt")
  write.table(df_fam, outFile2, quote = F, sep = "\t")
}
