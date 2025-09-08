pacman::p_load(edgeR, ggplot2, ggpubr, ggplotify, dplyr, purrr, gplots, rtracklayer, plyranges, reshape, kableExtra, ggExtra,
               ggdist, ggbreak, Rmisc, stringr, pals)
list.files()
#####################################################################################################################################################################
inPath_1 = "~/storage/Data/ms_evo_exwago/analysis/soloLTR/ltr_ann/heligmosomoides_bakeri/"
inFile_1 = "heligmosomoides_bakeri.counts_ltr_split.1827.txt.gz" # main table


outPath_1 = "ip_ub_tmm/" %>% paste0(inPath_1, .)
outPath_2 = "ip_ub_mir/" %>% paste0(inPath_1, .)
outPath_3 = "ip_ub_dge/" %>% paste0(inPath_1, .)
outPath_4 = "class/"%>% paste0(inPath_1, .)
outPath_5 = "family/" %>% paste0(inPath_1, .)
outPath_6 = "figures/" %>% paste0(inPath_1, .)

outs = c(outPath_1, outPath_2, outPath_3, outPath_4, outPath_5, outPath_6)
lapply(outs, dir.create)

conts = list(c("EV_IP", "EV_Ub"), c("Sup_IP", "Sup_Ub"))

names(conts) = lapply(conts, function(x){paste(x, collapse = "-")})
l_lfc = c(1, 1)
l_fdr = c(.05,.05)
names(l_lfc) = names(conts)
names(l_fdr) = names(conts)
dt_lst_nonorm = list()
dt_lst_normrrna = list()
tmm_ma = list()
mir_ma = list()

for ( cont in names(conts)) {
  #cont = "EV_IP-EV_Ub"
  print(cont)
  contrast_id = cont
  libs = conts[[cont]] %>% paste(collapse = "|")
  # Read table and create groups 
  counts = read.delim(paste0(inPath_1, inFile_1))
  counts = counts[,grepl(libs, colnames(counts))]
  colnames(counts) = sub("\\.trim.*", "", colnames(counts))
  groups = colnames(counts) %>% strsplit("_") %>%
    map(function(x){paste(x[c(1,2)], collapse = "_")}) %>%  unlist() %>% 
    factor()
  ################################################################################################################################################  
  lib.size = colSums(counts)
  # How many reads do we have?
  round((lib.size)/1e6, 2)
  # Create DGE and filter by expression
  dge = DGEList(counts=counts, group=groups, lib.size = lib.size)
  keep = filterByExpr(dge, min.count= 5)
  dge = dge[keep, , keep.lib.sizes=FALSE]
  print(table(keep))
  f.exp.counts = colSums(dge$counts)
  dge = calcNormFactors(dge, method = "TMM")
  ##############################################################################################################################################
  if (grepl("EV|Sup", libs)) {
    # As Kyriaki mentioned we should use a batch effect since this are paired worms
    batch = factor(paste0("batch",sub(".*(r.*)","\\1",rownames(dge$samples))))
    design = model.matrix(~0+dge$samples$group+batch)
    colnames(design) =  c(levels(dge$samples$group),levels(batch)[-1])
  }else{
    design = model.matrix(~0+dge$samples$group)
    colnames(design) =  c(levels(dge$samples$group))
  }
  ##############################################################################################################################################
  dge = estimateDisp(dge, design=design, robust=TRUE)
  plotBCV(dge)
  # As expected, we have variance in low expressed regions and low in highly expressed 
  # Let's fit glm
  ncont = colnames(design)[c(1,2)] %>% paste(collapse = "-")
  fit = glmFit(dge, design = dge$design, dispersion = dge$common.dispersion)
  contrast = makeContrasts(ncont,levels = dge$design)
  cont_fc = l_lfc[cont]
  cont_fdr = l_fdr[cont]
  gt = glmTreat(fit,  contrast=contrast, lfc = cont_fc )
  topTags(gt)
  dt = decideTests.DGEExact(gt, adjust.method="BH",  p.value=cont_fdr, lfc = cont_fc)
  tb_nonorm = table(dt)
  dt_lst_nonorm[[ncont]] = tb_nonorm
  
  topTable = topTags(gt, n=Inf)$table
  topTable = topTable[rownames(dt),]
  de_tbl = data.frame('Ave_CPM'=topTable$logCPM, 'log-FC'= topTable$logFC, 'cl'= dt)
  colnames(de_tbl) = c('Ave_CPM', 'logFC', 'cl')
  de_tbl$col = ifelse(de_tbl$cl == 0, "#C0C0C0", ifelse(de_tbl$cl == 1, "#5E2129", "#333333"))
  de_tbl$ord = ifelse(de_tbl$col == "#C0C0C0", 1, ifelse(de_tbl$col == "#333333", 2, 3))
  de_tbl = de_tbl[order(de_tbl$ord),]
  
  outTtab = "heligmosomoides_bakeri.other_fdr05_lfc1_tmm_topTbl.txt"
  outTtab = sub("other", ncont, outTtab)
  write.table(topTable, paste0(outPath_1,outTtab), quote = F, sep = "\t")
  
  outDEtbl = "heligmosomoides_bakeri.other_fdr05_lfc1_tmm_deTbl.txt"
  outDEtbl = sub("other", ncont, outDEtbl)
  write.table(de_tbl, paste0(outPath_1,outDEtbl), quote = F, sep = "\t")
  ##############################################################################################################################################
  de_tbl$exp = cont 
  de_tbl$exp = ifelse(grepl("Adult", cont),
                      "Adult exWAGO IP vs Adult total",
                      ifelse(grepl("EV", cont), 
                             "Vesicular exWAGO IP vs Unbound", 
                             "Non-vesicular exWAGO IP vs Unbound"))
  
  
  
  de_tbl$col = factor(de_tbl$col, levels = c("#5E2129", "#C0C0C0", "#333333"))
  
  p1 = ggplot(de_tbl, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote(~log[2]~' CPM')) + theme_test() +
    ylab(bquote(~log[2]~ 'FC IP/Ub')) +
    theme(axis.title.x = element_text(face = "bold", size = 16),
          axis.text.x = element_text( size = 16)) +
    theme(axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text( size = 16)) +
    geom_hline(aes(yintercept = 0), colour = "blue", 
               linewidth = 1, linetype="dashed", alpha=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size = 24,
                                    face = "bold")) + 
    scale_color_manual(values = adjustcolor(
      c("#155AA3", "#C0C0C0", "#333333"), alpha.f = .6),
      labels=c('IP enriched', 'not-DE', 'Unbound'), name=NULL) +
    theme(legend.position = "bottom", legend.text = element_text(size = 13)) +
    guides(color = guide_legend(override.aes = list(size=10))) +
    facet_wrap(~exp) +
    theme( strip.text = element_text(
      size = 12.5, face = "bold", colour = "black")) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches")) +
    ylim(c(-18, 11))
  p1 
  tmm_ma[[cont]] = p1
  #####################################################################################################################################################################
  
  # It is important to realize that IP vs Input are not fair contrasts.
  # These libraries are inherently compositionally biased due to the nature of small RNA sequencing.
  # This bias arises because Argonaute IPs selectively enrich for sRNAs with specific properties—
  # such as defined length distributions and 5′ nucleotide preferences.
  #
  # As a result, standard normalization methods that assume similar global distributions (like TMM)
  # may fail or require adjustment. To address this, we focus on a subset of regions (e.g., rRNA fragments)
  # that are not expected to be enriched in IPs and thus serve as a more stable reference for scaling.
  
  get_trimmed_category_counts = function(dge_full, category_pattern = "_rRNA_S", groups, 
                                         min_count = 5, logratioTrim = 0.3, sumTrim = 0.05) {
    
    # 1. Get the rows corresponding to the category of interest
    category_counts = dge_full$counts[grepl(category_pattern, rownames(dge_full$counts)), ]
    
    # 2. Get the corresponding N reads for this category
    lib_sizes_cat = colSums(category_counts)
    
    # 3. Build a new dge object but just for this category
    dge_cat = DGEList(counts = category_counts, group = groups, lib.size = lib_sizes_cat)
    
    # 4. Filter out lowly expressed regions. Since I'm using the libsize only for this category
    #    I do not expect to loose to many of this regions because of this filter. 
    keep = filterByExpr(dge_cat, min.count = min_count)
    dge_cat = dge_cat[keep, , keep.lib.sizes = FALSE]
    message("Features kept after filterByExpr:"); print(table(keep))
    
    # 5. Get the reference sample based on geometric mean
    lib_sizes = dge_cat$samples$lib.size
    geo_mean_lib_size = exp(mean(log(lib_sizes)))
    ref = which.min(abs(log(lib_sizes / geo_mean_lib_size)))
    
    # 6.Now I'll compare each sample against the reference. 
    test_samples = setdiff(seq_len(ncol(dge_cat$counts)), ref)
    
    trimming_masks = lapply(test_samples, function(test) {
      # I need to compare the logFC (M-value) between my reference and each test
      # I will also include a pseudocount to avoid inf with log2(0) 
      logR = log2((dge_cat$counts[, test] + 0.5) / lib_sizes[test]) -
        log2((dge_cat$counts[, ref] + 0.5) / lib_sizes[ref])
      
      # Now the average expression value between ref and each test (A-value)
      # Again, since this is log2 y need pseudocounts
      absE = 0.5 * (log2((dge_cat$counts[, test] + 0.5) / lib_sizes[test]) +
                      log2((dge_cat$counts[, ref] + 0.5) / lib_sizes[ref]))
      
      # This is the key to behave like tmm. I will remove (30%) extreme values in both up and down based on FC 
      keep_M = rank(logR) > length(logR) * logratioTrim & rank(logR) < length(logR) * (1 - logratioTrim)
      # Now the same but for expression where I'll just shrink by 5% of the lower and higer expressed regions
      keep_A = rank(absE) > length(absE) * sumTrim & rank(absE) < length(absE) * (1 - sumTrim)
      # Get the vector per lib
      keep <- keep_M & keep_A
      return(keep)
    })
    # 7. Now only keep stable regions in at least one contrast vs ref.
    # Why?
    # If we required regions to pass the trimming in *all* comparisons (i.e., Reduce(&, ...)),
    # we would likely discard most or all features due to natural biological and technical variability.
    # This is especially true for small RNA-seq data, where Argonaute IPs often have highly skewed composition.
    #
    # By using Reduce(|, ...), we are instead selecting regions that pass trimming in *at least one* comparison.
    # This approach is less strict, but more realistic—it allows us to retain regions that are generally stable,
    # even if they appear slightly extreme in some individual comparisons due to noise or enrichment.
    keep_any <- Reduce(`|`, trimming_masks)
    filtered_counts <- dge_cat$counts[keep_any, ]
    return(filtered_counts)
  }
  filtered_counts = get_trimmed_category_counts(
    dge_full = dge,               
    category_pattern = "_rRNA_S", 
    groups = groups               
  )
  
  normCounts = filtered_counts
  sumNormRNAs = colSums(normCounts)
  ########################################################################################################################################
  # Read again the table
  counts = read.delim(paste0(inPath_1, inFile_1))
  counts = counts[,grepl(libs, colnames(counts))]
  colnames(counts) = sub("\\.trim.*", "", colnames(counts))
  groups = colnames(counts) %>% strsplit("_") %>%
    map(function(x){paste(x[c(1,2,3)], collapse = "_")}) %>%  unlist() %>% 
    factor()
  lib.size = colSums(counts)
  # How many reads do we have?
  round((lib.size)/1e6, 2)
  # Create DGE and filter by expression
  dge = DGEList(counts=counts, group=groups, lib.size = lib.size)
  keep = filterByExpr(dge, min.count= 5)
  dge = dge[keep, , keep.lib.sizes=FALSE]
  print(table(keep))
  
  # Estimate norm factors based on tmm rRNA
  totalLibsize = colSums(dge$counts)
  normFactors = sumNormRNAs / totalLibsize
  normFactors = normFactors / (prod(normFactors)^(1/length(normFactors)))
  dge$samples$norm.factors = normFactors
  dge = estimateDisp(dge, design=design, robust=TRUE)
  plotBCV(dge)
  #####################################################################################################################################################################
  # Fit models
  fit = glmFit(dge, design=design, robust=T, dispersion = dge$common.dispersion)
  contrast = makeContrasts(ncont,levels = dge$design)
  gt = glmTreat(fit,  contrast=contrast, lfc = 1)
  topTags(gt)
  dt = decideTestsDGE(gt, adjust.method="BH", p.value=0.05, lfc = 1)
  normrrna_dt = table(dt)
  normrrna_dt
  topTable = topTags(gt, n=Inf)$table
  topTable = topTable[rownames(dt),]
  #####################################################################################################################################################################
  de_tbl_norm = data.frame('Ave_CPM'=topTable$logCPM, 
                           'log-FC'= topTable$logFC, 'cl'= dt)
  colnames(de_tbl_norm) = c('Ave_CPM', 'logFC', 'cl')
  de_tbl_norm$col = ifelse(de_tbl_norm$cl == 0, "#C0C0C0", 
                           ifelse(de_tbl_norm$cl == 1, "#15C1FB", "#C0C0C0"))
  de_tbl_norm$col = factor(de_tbl_norm$col, levels = c("#15C1FB", "#C0C0C0"))
  
  de_tbl_norm$exp = ifelse(grepl("Adult", cont),
                           "Adult exWAGO IP vs Adult total", 
                           ifelse(grepl("EV", cont),
                                  "Vesicular exWAGO IP vs Unbound",
                                  "Non-vesicular exWAGO IP vs Unbound"))

  outTtab = "heligmosomoides_bakeri.other_fdr01_lfc1_mir_topTbl.txt"
  outTtab = sub("other", ncont, outTtab)
  write.table(topTable, paste0(outPath_2,outTtab), quote = F, sep = "\t")
  
  outDEtbl = "heligmosomoides_bakeri.other_fdr05_lfc1_mir_deTbl.txt"
  outDEtbl = sub("other", ncont, outDEtbl)
  write.table(de_tbl_norm, paste0(outPath_2,outDEtbl), quote = F, sep = "\t")  

  ma = ggplot(de_tbl_norm, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote(~log[2]~' CPM')) + theme_test() +
    ylab(bquote(~log[2]~ 'FC IP/Total')) + 
    theme(axis.title.x = element_text(face = "bold", size = 16), 
          axis.text.x = element_text( size = 16)) +
    theme(axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text( size = 16)) +
    geom_hline(aes(yintercept = 0), colour = "blue", 
               linewidth = 1, linetype="dashed", alpha=0.7)  +
    theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) +  
    scale_color_manual(values = adjustcolor(c("#155AA3", "#C0C0C0"), alpha.f = .6),
                       labels=factor(c('exWAGO enriched', 'Unbound'), 
                                     levels = c("Unbound","exWAGO enriched")), 
                       name=NULL) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
    theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size=10)))  +
    facet_wrap(~ exp) +
    theme( strip.text = element_text(size = 12.5, face = "bold", colour = "black")) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches")) +
    ylim(c(-15, 15))
  ma 
  
  mir_ma[[cont]] = ma
  
  outDge = "heligmosomoides_bakeri.other_miR_norm_dge.Rds"
  outDge = sub("other",ncont, outDge)
  saveRDS(dge,paste0(outPath_3, outDge))
  
  
  ################################################################################################################################################################################
  categories = c("DNA.*As|RC.*As", "DNA.*S|RC.*S", "LINE.*As|PLE.*As|Retrop*As", "LINE.*S|PLE.*S|Retrop*S", "LTR.*As", "LTR.*S",
                 "Unk", "Low|Simple|Sat|SINE|MITE", "lincRNA_As", "miRNA_S", 
                 "^tRNA_S", "_rRNA|yRNA|snRNA|snoRNA|other_ncRNA|piRNA|lncRNA|miRNA_As|lincRNA_S|^tRNA_As", 
                 "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  
  names(categories) = c("DNA_As", "DNA_S", "LINE_As", "LINE_S", "LTR_As", "LTR_S", "Unknown", 
                        "Other_repeat", "lincRNA_As", "miRNA_S", "tRNA_S", 
                        "other_ncRNA", "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  allFams = list()
  for (i in names(categories)) {
    tmp.cat = categories[i]
    tmp.sqs = dge$counts[grepl(tmp.cat, rownames(dge$counts)),] %>% rownames()
    allFams[[i]] =tmp.sqs 
  }
  design = dge$design
  allFamsIdx = lapply(allFams, function(x) rownames(dge$counts) %in% x)
  allFamsIdx = allFamsIdx[sapply(allFamsIdx, sum) > 5]
  cameraRes = camera(dge, allFamsIdx, design, contrast)
  
  outClass = "heligmosomoides_bakeri.other_camera_enr_class_simple_miR.txt"
  outClass = sub("other", ncont, outClass) 
  
  write.table(cameraRes, 
              paste0(outPath_4, outClass),
              sep = "\t", row.names = T, col.names = T, quote = F)
  
  ##########################################################################################################################################################################
  aveCPM = 2 ** aveLogCPM(dge, normalized.lib.sizes = T, prior.count = .01)
  names(aveCPM) = rownames(dge$counts)
  
  aveCPMByfam = split(aveCPM, sub(":.*", "", names(aveCPM)))
  
  categories =  names(aveCPMByfam)[grepl("DNA|LINE|LTR|Unk",names(aveCPMByfam))]
  names(categories) = categories
  
  allFams = list()
  for (i in names(categories)) {
    tmp.cat = categories[i]
    tmp.sqs = dge$counts[grepl(tmp.cat, rownames(dge$counts)),] %>% rownames()
    allFams[[i]] =tmp.sqs 
  }
  
  allFamsIdx = lapply(allFams, function(x) rownames(dge$counts) %in% x)
  allFamsIdx = allFamsIdx[sapply(allFamsIdx, sum) > 5]
  cameraRes = camera(dge, allFamsIdx, design, contrast)
  cameraRes = cameraRes[match(names(categories), rownames(cameraRes)),]
  rownames(cameraRes) = names(categories)
  cameraRes$Feature = ifelse(grepl("DNA|RC", rownames(cameraRes)), "DNA", 
                             ifelse(grepl("LINE|PLE", rownames(cameraRes)), "LINE",
                                    ifelse(grepl("LTR", rownames(cameraRes)), "LTR", "Unknown")))
  
  outClass = "heligmosomoides_bakeri.other_camera_enr_fam_simple_miR.txt"
  outClass = sub("other", ncont, outClass) 
  
  write.table(cameraRes, 
              paste0(outPath_5, outClass),
              sep = "\t", row.names = T, col.names = T, quote = F)
}
#######################################################################################################################################################################
tmm_p = ggarrange(tmm_ma$`EV_IP-EV_Ub`, tmm_ma$`Sup_IP-Sup_Ub`, 
          common.legend = T, legend = "bottom")
mir_p = ggarrange(mir_ma$`EV_IP-EV_Ub`, mir_ma$`Sup_IP-Sup_Ub`, 
          common.legend = T, legend = "bottom")

all_ma =  ggarrange(tmm_p, mir_p, nrow = 2)
ggsave(filename = "heligmosomoides_bakeri.EV_SUP_MA.png", device = "png",
       dpi = 300, width = 10, height = 10, path = outPath_6)
