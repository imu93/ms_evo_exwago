setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/rrna_norm/")
pacman::p_load(dplyr, reshape2)
files  = list.files(pattern = "deTbl.txt")
tabs = lapply(files, read.delim) 
names(tabs)=  c("A_ceylanicum", "PPW-1", "SAGO-1", "SAGO-2", "H_bakeri", "H_polygyrus",
                "N_brasiliensis", "T_circumcincta")

lapply(tabs, function(x){
  #x = tabs$`PPW-1`
  tRNA = x[grepl("tRNA_S", rownames(x)),]
  all_trna = nrow(tRNA)
  up_trna = tRNA[tRNA$cl == 1,] %>% nrow
  down_trna = tRNA[tRNA$cl != 1,] %>% nrow
  prop_up = 100*(up_trna/all_trna)
  prop_down = 100*(down_trna/all_trna)
  df = data.frame("N_features"=all_trna, "N_up"=up_trna , "N_down"=down_trna, "prop_ip"=prop_up, "prop_down"=prop_down)
}) %>% do.call(rbind, .) -> tRNA_df

tRNA_df$prop_down %>% mean()

lapply(tabs, function(x){
  #x = tabs$`PPW-1`
  tRNA = x[grepl("miRNA_S", rownames(x)),]
  all_trna = nrow(tRNA)
  up_trna = tRNA[tRNA$cl == 1,] %>% nrow
  down_trna = tRNA[tRNA$cl != 1,] %>% nrow
  prop_up = 100*(up_trna/all_trna)
  prop_down = 100*(down_trna/all_trna)
  df = data.frame("N_features"=all_trna, "N_up"=up_trna , "N_down"=down_trna, "prop_ip"=prop_up, "prop_down"=prop_down)
}) %>% do.call(rbind, .) -> miRNA_df
miRNA_df$prop_down %>% mean()

lapply(tabs, function(x){
  #x = tabs$`PPW-1`
  tRNA = x[grepl("snoRNA_S", rownames(x)),]
  all_trna = nrow(tRNA)
  up_trna = tRNA[tRNA$cl == 1,] %>% nrow
  down_trna = tRNA[tRNA$cl != 1,] %>% nrow
  prop_up = 100*(up_trna/all_trna)
  prop_down = 100*(down_trna/all_trna)
  df = data.frame("N_features"=all_trna, "N_up"=up_trna , "N_down"=down_trna, "prop_ip"=prop_up, "prop_down"=prop_down)
}) %>% do.call(rbind, .) -> snoRNA_df
snoRNA_df$prop_down %>% mean()

lapply(tabs, function(x){
  #x = tabs$`PPW-1`
  tRNA = x[grepl("snRNA_S", rownames(x)),]
  all_trna = nrow(tRNA)
  up_trna = tRNA[tRNA$cl == 1,] %>% nrow
  down_trna = tRNA[tRNA$cl != 1,] %>% nrow
  prop_up = 100*(up_trna/all_trna)
  prop_down = 100*(down_trna/all_trna)
  df = data.frame("N_features"=all_trna, "N_up"=up_trna , "N_down"=down_trna, "prop_ip"=prop_up, "prop_down"=prop_down)
}) %>% do.call(rbind, .) -> snRNA_df
snRNA_df$prop_down %>% sd()

lapply(tabs, function(x){
  #x = tabs$`PPW-1`
  tRNA = x[grepl("rRNA_S", rownames(x)),]
  all_trna = nrow(tRNA)
  up_trna = tRNA[tRNA$cl == 1,] %>% nrow
  down_trna = tRNA[tRNA$cl != 1,] %>% nrow
  prop_up = 100*(up_trna/all_trna)
  prop_down = 100*(down_trna/all_trna)
  df = data.frame("N_features"=all_trna, "N_up"=up_trna , "N_down"=down_trna, "prop_ip"=prop_up, "prop_down"=prop_down)
}) %>% do.call(rbind, .) -> snRNA_df

lst = list(tRNA_df, miRNA_df, snRNA_df, snoRNA_df)
names(lst) = c("tRNA","miRNA", "snRNA", "snoRNA")

df_lst =list()
for (i in names(lst)) {
  tmp.f = lst[[i]]
  tmp.f$Feature = i
  tmp.f$species = rownames(tmp.f)
  df_lst[[i]] = tmp.f
}

names(df_lst) = NULL

df2plot = do.call(rbind, df_lst)
write.table(df2plot,"nematode_ncrna_ub_srna_fams.txt", sep = "\t", quote = F)
