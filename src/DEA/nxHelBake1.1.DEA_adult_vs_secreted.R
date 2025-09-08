setwd("~/storage/Data/ms_evo_exwago/analysis/soloLTR/ltr_ann/heligmosomoides_bakeri/")
pacman::p_load(edgeR, rtracklayer, ggplot2, ggpubr, kableExtra, gplots, RColorBrewer, grDevices, dplyr, purrr, ggExtra, stringr)
conts = c("EV_IP-Adult_IP", "Sup_IP-Adult_IP", "Adult_IP-EV_IP-Sup_IP")
ma_colors = list(c("#3EC1E6","#BCBBBD", "#225999"), c("#3EE662","#BCBBBD","#225999"), c(c("#225999", "#BCBBBD", "#3EE69B"))) 
names(ma_colors) = conts
###################################################################################################################################################################################################################
ma_lis = list() 
de_tets_lst = list()
control_norm_regs_lst = list()
cameraRes_fam_lis = list()
cameraRes_lis = list()
for (cont in conts) {
  #cont = "Sup_IP-Adult_IP"
  comparison = cont
  libs = gsub("-", "|", comparison)
  counts = read.table("nxHelBake1.1.counts_adult_vs_secretetion_filtred_ltr_split_1827.txt")
  colnames(counts) = sub("\\.trim.*", "", colnames(counts))
  head(counts)
  counts = counts[,grepl(libs, colnames(counts))]
  lib.size = colSums(counts)
  round((lib.size)/1e6, 2)
  # I'll use a MDS plot to have an idea of the influence of the batch effect
  # But first I need to build my dge object
  group = colnames(counts) %>% strsplit("_") %>%
    map(function(x){paste(x[c(1,2)], collapse = "_")}) %>%  unlist() %>% 
    factor()
  dge = DGEList(counts=counts, group=group, lib.size = lib.size)
  keep = filterByExpr(dge, min.count= 5)
  dge = dge[keep, , keep.lib.sizes=FALSE]
  print(table(keep))
  f.exp.counts = colSums(dge$counts)
  batch = factor(paste0("batch",sub(".*(r.*)","\\1",rownames(dge$samples))))
  design = model.matrix(~0+dge$samples$group+batch)
  colnames(design) =  c(levels(dge$samples$group),levels(batch)[-1])
  ###################################################################################################################################################################################################################
  dge = calcNormFactors(dge, method = "TMM")
  # As a techincal note: use alpha in rich colors  to have transparency
  colors = rich.colors(length(levels(group)), alpha = .5)
  names(colors) = levels(group)
  par(mfrow=c(1,2))
  plotMDS(dge, col=colors[dge$samples$group], pch = 16, cex=3,main="MDS plot of H. bakeri exWAGO IP")
  plot.new()
  legend("topright", legend = names(colors), pch = 16, col = unique(colors))
  ###################################################################################################################################################################################################################
  dge = estimateDisp(dge, design=design, robust=TRUE)
  #dge$counts["introns_As:1437670",]
  plotBCV(dge)
  fit = glmFit(dge, design, dispersion = dge$common.dispersion)
  if (comparison == "Adult_IP-EV_IP-Sup_IP") {
    comparison = "Adult_IP-((EV_IP+Sup_IP)/2)"
  }
  contrast = makeContrasts(comparison,levels = dge$design)
  gt = glmTreat(fit, contrast=contrast, lfc = log2(1.2))
  topTags(gt)
  dt = decideTests.DGEExact(gt, adjust.method="BH", p.value=0.05, lfc = log2(1.2))
  de_tets_lst[[cont]] =  table(dt)
  topTable = topTags(gt, n=Inf)$table
  topTable = topTable[rownames(dt),]
  #dge$counts["introns_As:1437670",]
  ###################################################################################################################################################################################################################
  de_tbl = data.frame('Ave_CPM'=topTable$logCPM, 'log-FC'= topTable$logFC, 'cl'= dt)
  colnames(de_tbl) = c('Ave_CPM', 'logFC', 'cl')
  de_tbl$col = ifelse(de_tbl$cl == 0, "#C0C0C0", ifelse(de_tbl$cl == 1, "#5E2129", "#333333"))
  de_tbl$ord = ifelse(de_tbl$col == "#C0C0C0", 1, ifelse(de_tbl$col == "#333333", 2, 3))
  de_tbl = de_tbl[order(de_tbl$ord),]
  
  outFile = "nxHelBake1.1.ewago_other.lfc1_fdr05.deTbl_fil.txt"
  outFile = sub("other", cont, outFile)
  write.table(de_tbl, outFile, quote = F, sep = "\t")
  
  outFile2 = "nxHelBake1.1.ewago_other.lfc1_fdr05.topTbl_fil.txt"
  outFile2 = sub("other", cont, outFile2)
  write.table(topTable, outFile2, quote = F, sep = "\t")
  
  outFile4 = "nxHelBake1.1_exwago_other_filtered_dge.Rds"
  outFile4 = sub("other", cont, outFile4)
  saveRDS(dge, outFile4)
  ###################################################################################################################################################################################################################
  de_tbl$col = factor(de_tbl$col, levels = c("#5E2129", "#C0C0C0", "#333333"))
  cols = ma_colors[[cont]]
  conds = strsplit(comparison, "-") %>% unlist()
  de_tbl$cont = paste("exWAGO", comparison)
  de_tbl$cont = ifelse(grepl("Adult_IP-EV_IP-Sup_IP", cont), sub("EV_IP-Sup_IP", "Secreted_IPs", cont), de_tbl$cont)
  de_tbl$cont = ifelse(grepl("Sup_IP-Adult_IP", cont), sub("Sup_IP", "Non_Vesicular_IP", cont), de_tbl$cont)
  de_tbl$cont = ifelse(grepl("EV_IP-Adult_IP", cont), sub("EV_IP", "Vesicular_IP", cont), de_tbl$cont)
  sp_cont = strsplit(de_tbl$cont %>% unique(), "-") %>% unlist()
  
  
  p = ggplot(de_tbl, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote(~log[2]~' CPM')) + theme_test() +
    ylab(bquote(~log[2]~ ' FC')) + theme(axis.title.x = element_text( size = 20), axis.text.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),axis.text.y = element_text(size = 20)) +
    geom_hline(aes(yintercept = 0), colour = "red", linewidth = 1.4, linetype="dashed", alpha=0.4) + 
    theme( strip.text = element_text(size = 18, colour = "black"))  +
    theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) + 
    scale_color_manual(values = adjustcolor(cols, alpha.f = .6), labels=c(sp_cont[1], 'shared', sp_cont[2]), name=NULL) +
    ylim(c(-14,14)) +
    facet_wrap(~cont)
  
  fp = p + theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size=10)))
  
  ma_lis[[cont]] = fp
  ####################################################################################################################################################################
  
  ann = import("~/storage/Data/ms_evo_exwago/raw/ShortStack/hBake/heligmosomoides_bakeri.exWAGO_DE_LTRsplit.gff3")
  ann_fl = ann[grepl("_fl", ann$LTR_type)]
  rownames(de_tbl) = ifelse(rownames(de_tbl) %in% ann_fl$ID, sub("_", "_fl_", rownames(de_tbl)), rownames(de_tbl))
  cts = c("soloLTR", "LTR_fl", "int_fl")
  
  df_lst = list()
  for (i in cts) {
    df = de_tbl
    df$col = ifelse(grepl(i, rownames(df)), "#144675", "#C0C0C0")
    df$ord = ifelse(df$col == "#C0C0C0", 1,2)
    df$class = i
    df = df[order(df$ord),]
    df_lst[[i]] = df
  }
  names(df_lst) = NULL
  df = do.call(rbind, df_lst)
  
  
  df$class = factor(df$class, levels = c("soloLTR", "LTR_fl", "int_fl"))
  control_norm = ggplot(df, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point(show.legend = F) +  xlab(bquote(~log[2]~'CPM')) + theme_minimal() + facet_wrap(~class, nrow=2) + theme_light() +
    theme( strip.text = element_text(size = 11, face = "bold", colour = "black")) +
    ylab(bquote(~log[2]~ 'FC')) + theme(axis.title.x = element_text(face = "bold", size = 24), axis.text.x = element_text(face = "bold", size = 24)) +
    theme(axis.title.y = element_text(face = "bold", size = 24),axis.text.y = element_text(face = "bold", size = 24)) +
    geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
    scale_color_manual(values = adjustcolor(c("#144675", "#BCBBBD"), alpha.f = .6)) + ylim(c(-15,18))
  control_norm_regs_lst[[cont]] = control_norm
  
}


ggarrange(ma_lis$`EV_IP-Adult_IP`, 
          ma_lis$`Sup_IP-Adult_IP`,
          ncol=2, nrow = 1)

ggsave(filename = "nxHelBake1.1.exWAGO_secretion.specific.png", 
       device = "png", dpi=300,
       path = "./figures/", 
       plot = last_plot(), width = 12, height = 5)

ggarrange(ma_lis$`Adult_IP-EV_IP-Sup_IP`, ma_lis$`EV_IP-Adult_IP`, 
          ma_lis$`Sup_IP-Adult_IP`, 
          ncol=3)

ggsave(filename = "nxHelBake1.1.exWAGO_secretion.specific.png", 
       device = "png", dpi=300,
       path = "./figures/", 
       plot = last_plot(), width = 18, height = 5)


ev_spe = ma_lis$`EV_IP-Adult_IP`
sup_spe = control_norm_regs_lst$`Sup_IP-Adult_IP`
adult_spe = control_norm_regs_lst$`Adult_IP-EV_IP-Sup_IP`

ggsave(filename = "nxHelBake1.1.exWAGO_ev_spe_gcats.png", 
       device = "png", dpi=300, path = "./figures/", plot = ev_spe, width = 6.5 , height = 6)

ggsave(filename = "nxHelBake1.1.exWAGO_sup_spe_gcats.png", 
       device = "png", dpi=300, path = "./figures/", plot = sup_spe, width = 5, height = 6)

ggsave(filename = "nxHelBake1.1.exWAGO_adult_spe_gcats.png", 
       device = "png", dpi=300, path = "./figures/", plot = adult_spe, width = 5, height = 6)


names(de_tets_lst) = NULL
tbl = do.call(rbind, de_tets_lst)
