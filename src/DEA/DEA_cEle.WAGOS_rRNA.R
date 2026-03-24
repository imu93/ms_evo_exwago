pacman::p_load(edgeR, ggplot2, ggpubr, ggplotify, dplyr, purrr, gplots, rtracklayer, plyranges, reshape, kableExtra, ggExtra,
               Rmisc, stringr, pals, stringi, grDevices, Biostrings, BSgenome)

inPath_1 = "~/storage/Data/ms_evo_exwago/evo_exwago_guides/count_tables/"
inFile_1 = "caenorhabditis_elegans.counts.1827.txt.gz" # main table


outPath_1 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/tmm_gbe/"
outPath_2 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/rrna_norm_gbe/"
outPath_3 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/dge_gbe/"
outPath_4 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/camera/class_gbe/" 
outPath_5 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/camera/family_gbe/"
outPath_6 = "~/storage/Data/ms_evo_exwago/analysis/figures_gbe/"
outPath_7 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/dge_tmm_gbe/"

drs = c(outPath_1, outPath_2, outPath_3, outPath_4, outPath_5, outPath_6, outPath_7)
sapply(drs, dir.create)



counts = read.delim(paste0(inPath_1, inFile_1))
counts = counts[,!grepl("total", colnames(counts))]
counts = counts[,!grepl("rnase|mock",colnames(counts))]
conds = colnames(counts) %>% str_replace("Ce_Adult_|Ce_L4_", "") %>%
  str_replace("\\..*", "") %>% 
  stri_replace_last_regex("_[^_]*$", "") %>% unique()
conts = split(conds, sub("^[^_]*_", "", conds))
conts = lapply(conts, rev)
names(conts) = lapply(conts, function(x){paste(x,collapse = "-")})
#conts = conts[!grepl("alg|^Ip_rde1", names(conts))]




l_lfc = rep(1, length(conts))
l_fdr = rep(.05, length(conts))

names(l_lfc) = names(conts)
names(l_fdr) = names(conts)
# For the final version of the paper I will use the same FC and p-value 
# For all libraris
#l_lfc["Ip_ppw1-input_ppw1"] = 2
#l_lfc["Ip_ergo1-input_ergo1"] = 2
#l_lfc["Ip_wago1-input_wago1"] = 2
#l_lfc["Ip_wago4-input_wago4"] = 2
#l_fdr["Ip_wago4-input_wago4"] = .01
#l_lfc["Ip_prg1-input_prg1"] = 2
#l_fdr["Ip_prg1-input_prg1"] = .01
#l_lfc["Ip_csr1-input_csr1"] = 2
#l_fdr["Ip_csr1-input_csr1"] = .01
#l_fdr["Ip_alg1-input_alg1"] = .01
#l_fdr["Ip_alg2-input_alg2"] = .01


dt_lst_nonorm = list()
dt_lst_normrrna = list()
tmm_ma = list()
norm_miR_simp_list = list()
mds_list = list()

for ( cont in names(conts)) {
  #cont = "Ip_ergo1-input_ergo1" 
  print(cont)
  contrast_id = cont
  libs = conts[[cont]][1] %>% str_replace("Ip_", "") %>% paste0("_")
  libs = paste0("_", libs)
  # Read table and create groups 
  counts = read.delim(paste0(inPath_1, inFile_1))
  counts = counts[,!grepl("rnase|mock",colnames(counts))]
  counts = counts[,grepl(libs, colnames(counts))]
  colnames(counts) = sub("\\.trim.*", "", colnames(counts))
  
   
  groups = colnames(counts) %>% strsplit("_") %>%
    map(function(x){paste(x[c(3,4)], collapse = "_")}) %>%  unlist() %>% 
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
  colors = rich.colors(2) %>% rep(2) %>% rev() %>% adjustcolor(alpha.f = .2)
  shapes = ifelse(grepl("input", levels(dge$samples$group)), 15, 16)
  names(colors) = levels(groups)
  names(shapes) = levels(groups)
  mds = plotMDS(dge, col=colors[dge$samples$group], pch= shapes[dge$samples$group], cex = 4)
  
  df_mds =  data.frame("dim1"=mds$x, "dim2"=mds$y, "cols"=colors[dge$samples$group])
  mds_gp = ggplot(df_mds,aes(dim1, dim2, color=cols)) + 
    geom_point(shape=16, size = 12) +
    theme_light() +
    scale_color_manual(values = adjustcolor(rep(rich.colors(2), 2), alpha.f = .7),
                       name="", labels=c("Ip", "Input")) +
    xlab(paste0("Leading logFC dim 1 ", round(mds$var.explained[1]*100), "%")) +
    ylab(paste0("Leading logFC dim 2 ", round(mds$var.explained[2]*100), "%")) +
    theme(axis.title.x = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", size = 12)) +
    theme(axis.title.y = element_text(face = "bold", size = 12),
          axis.text.y = element_text(face = "bold", size = 12),
          legend.position = "bottom",
          legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size=14))) +
    ggtitle(cont) + theme(plot.title = element_text(hjust = 0.5,
                                                    size = 15, face = "bold"))
  
  mds_list[[cont]] = mds_gp
  ###############################################################################################################################################
  design = model.matrix(~0+dge$samples$group)
  colnames(design) =  levels(dge$samples$group)
  
  ##############################################################################################################################################
  dge = estimateDisp(dge, design=design, robust=TRUE)
  plotBCV(dge)
  # As expected, we have variance in low expressed regions and low in highly expressed 
  # Let's fit glm
  ncont = colnames(design)[c(2,1)] %>% paste(collapse = "-")
  fit = glmFit(dge, design = dge$design, dispersion = dge$common.dispersion)
  contrast = makeContrasts(ncont,levels = dge$design)
  cont_fc = l_lfc[cont]
  cont_fdr = l_fdr[cont]
  gt = glmTreat(fit,  contrast=contrast, lfc = 1)
  topTags(gt)
  dt = decideTestsDGE(gt, adjust.method="BH",  p.value=0.05, lfc = 1)
  tb_nonorm = table(dt)
  dt_lst_nonorm[[cont]] = tb_nonorm
  topTable = topTags(gt, n=Inf)$table
  topTable = topTable[rownames(dt),]
  outToptable = paste0(outPath_1, "caenorhabditis_elegans.other_tmm_topTbl.txt")
  outToptable = sub("other", cont, outToptable)
  write.table(topTable,outToptable, quote = F, sep = "\t")
  
  de_tbl = data.frame('Ave_CPM'=topTable$logCPM, 'log-FC'= topTable$logFC, 'cl'= dt)
  colnames(de_tbl) = c('Ave_CPM', 'logFC', 'cl')
  de_tbl$col = ifelse(de_tbl$cl == 0, "#C0C0C0", ifelse(de_tbl$cl == 1, "#5E2129", "#C0C0C0"))
  de_tbl$ord = ifelse(de_tbl$col == "#C0C0C0", 1, ifelse(de_tbl$col == "#C0C0C0", 2, 3))
  de_tbl = de_tbl[order(de_tbl$ord),]
  de_tbl$exp = "C. elegans IP vs Input"
  de_tbl$exp = factor(de_tbl$exp)
  levels(de_tbl$exp) = c(expression(paste(italic("C. elegans "),
                                          ncont)))
  de_tbl$exp = sub("ncont", ncont, de_tbl$exp)
  
  
  
  outDEtbl = paste0(outPath_1, "caenorhabditis_elegans.other_tmm_deTbl.txt")
  outDEtbl = sub("other", cont, outDEtbl)
  write.table(de_tbl,outDEtbl, quote = F, sep = "\t")
  
  
  p1 = ggplot(de_tbl, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote(log[2]~' CPM')) + theme_test() +
    ylab(bquote(~log[2]~ 'FC IP/Total')) +
    theme(axis.title.x = element_text(size = 20),
          axis.text.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 20)) +
    geom_hline(aes(yintercept = 0), colour = "blue",
               linewidth = 1, linetype="dashed", alpha=0.4) +
    facet_wrap(~ exp, labeller = label_parsed) +
    scale_color_manual(values = adjustcolor(c("#155AA3","#C0C0C0"), alpha.f = .6), 
                       labels=c('IP', 'not-DE'), name=NULL) +
    theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    theme( strip.text = element_text(size = 18, colour = "black")) +
    guides(color = guide_legend(override.aes = list(size=10))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches"))
  p1
  outFile_tmm="caenorhabditis_elegans.other_tmm_dge.Rds"
  outFile_tmm=sub("other", cont, outFile_tmm)
  saveRDS(dge,paste0(outPath_7,outFile_tmm))
  tmm_ma[[cont]] = p1 
  
  
  
  
  
  ##############################################################################################################################################
  
  cats = c("rRNA_S", "miRNA_S", "tRNA_S")
  names(cats) = c("rRNA_S", "miRNA_S", "tRNA_S")
  df_lst = list()
  for (i in names(cats)) {
    df = de_tbl
    df$col = ifelse(grepl(i, rownames(df)), "#5E2129", "#C0C0C0")
    df$ord = ifelse(df$col == "#C0C0C0", 1,2)
    df = df[order(df$ord),]
    df$bty = i
    df_lst[[i]] = df
  }
  df = do.call(rbind, df_lst)
  cts_1 = ggplot(df, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point(show.legend = F) + facet_wrap(~ bty) +
    xlab(bquote(log[2]~' CPM')) + theme_test() +
    ylab(bquote(~log[2]~ 'FC IP/Total')) +
    theme(axis.title.x = element_text(size = 20),
          axis.text.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 20)) +
    geom_hline(aes(yintercept = 0), colour = "blue",
               linewidth = 1, linetype="dashed", alpha=0.4) +
    scale_color_manual(values = adjustcolor(c("#155AA3","#C0C0C0"), alpha.f = .6)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    theme( strip.text = element_text(size = 18, colour = "black")) +
    guides(color = guide_legend(override.aes = list(size=10))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches"))
  
  cont
  write.table(df, paste0("~/storage/Data/ms_evo_exwago/analysis/de_tables/tmm_norm_cats/",
                        paste0(cont,"_nomcats.txt")), sep = "\t", quote = F, col.names = T, row.names = T)
  
  
  ##############################################################################################################################################
  # It is important to realize that IP vs Input are not fair contrasts.
  # These libraries are inherently compositionally biased due to the nature of small RNA sequencing.
  # This bias arises because Argonaute IPs selectively enrich for sRNAs with specific properties—
  # such as defined length distributions and 5′ nucleotide preferences.
  #
  # As a result, standard normalization methods that assume similar global distributions (like TMM)
  # may fail or require adjustment. To address this, we focus on a subset of regions (e.g., rRNA fragments)
  # that are not expected to be enriched in IPs and thus serve as a more stable reference for scaling.
  
  get_trimmed_category_counts = function(dge_full, category_pattern = "^rRNA_S", groups, 
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
    category_pattern = "^rRNA_S", 
    groups = groups               
  )
  
  normCounts = filtered_counts
  sumNormRNAs = colSums(normCounts)

  ########################################################################################################################################
  
  counts = read.delim(paste0(inPath_1, inFile_1))
  counts = counts[,!grepl("rnase|mock",colnames(counts))]
  counts = counts[,grepl(libs, colnames(counts))]
  colnames(counts) = sub("\\.trim.*", "", colnames(counts))
  
  groups = colnames(counts) %>% strsplit("_") %>%
    map(function(x){paste(x[c(3,4)], collapse = "_")}) %>%  unlist() %>% 
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
  cont_fc = l_lfc[cont]
  cont_fdr = l_fdr[cont]
  gt = glmTreat(fit,  contrast=contrast, lfc = cont_fc)
  topTags(gt)
  dt = decideTestsDGE(gt, adjust.method="BH", p.value=cont_fdr, lfc = cont_fc)
  normrrna_dt = table(dt)
  dt_lst_normrrna[[cont]] = normrrna_dt
  topTable = topTags(gt, n=Inf)$table
  topTable = topTable[rownames(dt),]
  
  outToptable = paste0(outPath_2, "caenorhabditis_elegans.other_rRNA_topTbl.txt")
  outToptable = sub("other", cont, outToptable)
  write.table(topTable,outToptable, quote = F, sep = "\t")

  
  #####################################################################################################################################################################
  de_tbl_norm = data.frame('Ave_CPM'=topTable$logCPM, 
                           'log-FC'= topTable$logFC, 'cl'= dt)
  colnames(de_tbl_norm) = c('Ave_CPM', 'logFC', 'cl')
  de_tbl_norm$col = ifelse(de_tbl_norm$cl == 0, "#C0C0C0", 
                           ifelse(de_tbl_norm$cl == 1, "#15C1FB", "#C0C0C0"))
  
  de_tbl_norm$col = ifelse(rownames(de_tbl_norm) %in% rownames(filtered_counts), 
                           "#DBBF10", de_tbl_norm$col)
  de_tbl_norm$col = factor(de_tbl_norm$col, levels = c("#15C1FB", "#C0C0C0",  "#DBBF10"))
  de_tbl_norm$ord = 1
  de_tbl_norm$ord = ifelse(de_tbl_norm$col == "#15C1FB", 2, de_tbl_norm$ord)
  de_tbl_norm$ord = ifelse(de_tbl_norm$col == "#DBBF10", 3, de_tbl_norm$ord)
  de_tbl_norm$ord = factor(de_tbl_norm$ord)
  
  de_tbl_norm = de_tbl_norm[order(de_tbl_norm$ord),]
  de_tbl_norm[de_tbl_norm$col == "#DBBF10",]
  
  de_tbl_norm$exp = cont
  
  outDEtbl = paste0(outPath_2, "caenorhabditis_elegans.other_rRNA_deTbl.txt")
  outDEtbl = sub("other", cont, outDEtbl)
  write.table(de_tbl_norm,outDEtbl, quote = F, sep = "\t")
  #####################################################################################################################################################################
  ma = ggplot(de_tbl_norm, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote(~log[2]~' CPM')) + theme_test() +
    ylab(bquote(~log[2]~ 'FC IP/Total')) + 
    theme(axis.title.x = element_text(size = 20), 
          axis.text.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),
          axis.text.y = element_text( size = 20)) +
    geom_hline(aes(yintercept = 0), colour = "blue", 
               linewidth = 1, linetype="dashed", alpha=0.7)  +
    theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) +  
    scale_color_manual(values = adjustcolor(c("#155AA3", "#C0C0C0", "red"), alpha.f = .6),
                       labels=factor(c('IP', 'Unbound', "rRNA_S"), 
                                     levels = c("Unbound","IP", "rRNA_S")), 
                       name=NULL) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
    theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size=10)))  +
    facet_wrap(~ exp) +
    theme( strip.text = element_text(size = 18, colour = "black")) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches")) 
  ma 
 
  norm_miR_simp_list[[cont]]= ma
  outDGE = paste0(outPath_3, "caenorhabditis_elegans.other_rRNA_norm_dge.Rds")
  outDGE = sub("other", cont, outDGE)
  saveRDS(dge, outDGE)
  #####################################################################################################################################################################
  categories = c("DNA.*As|RC.*As", "DNA.*S|RC.*S", "LINE.*As|PLE.*As|Retrop*As", "LINE.*S|PLE.*S|Retrop*S", "LTR.*As", "LTR.*S",
                 "Unk", "Low|Simple|Sat|SINE|MITE", "lincRNA_As", "miRNA_S", "_rRNA_S",
                 "^tRNA_S","_rRNA_As|yRNA|snRNA|snoRNA|other_ncRNA|piRNA|lncRNA|miRNA_As|lincRNA_S|^tRNA_As", 
                 "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  
  names(categories) = c("DNA_As", "DNA_S", "LINE_As", "LINE_S", "LTR_As", "LTR_S", "Unknown", 
                        "Other_repeat", "lincRNA_As", "miRNA_S", "rRNA_S", "tRNA_S", 
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
  
  outFile_camClass = paste0(outPath_4, "caenorhabditis_elegans.other_camera_enr_class_simple_rRNA.txt")
  outFile_camClass = sub("other", cont, outFile_camClass)
  write.table(cameraRes, outFile_camClass , sep = "\t", row.names = T, col.names = T, quote = F)
  #####################################################################################################################################################################
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
  cameraRes$Ip = cont
  cameraRes$Feature = ifelse(grepl("DNA|RC", rownames(cameraRes)), "DNA", 
                             ifelse(grepl("LINE|PLE", rownames(cameraRes)), "LINE",
                                    ifelse(grepl("LTR", rownames(cameraRes)), "LTR", "Unknown")))
  cameraRes = cameraRes[!grepl("^NA", rownames(cameraRes)),]

  outFile_camFam = paste0(outPath_5, "caenorhabditis_elegans.other_camera_enr_family_simple_rRNA.txt")
  outFile_camFam = sub("other", cont, outFile_camFam)
  write.table(cameraRes, outFile_camFam , sep = "\t", row.names = T, col.names = T, quote = F)
}
#####################################################################################################################################################################

p_tmm = ggarrange(plotlist = tmm_ma, ncol = 3, nrow = 6,common.legend = T, legend = "bottom")

#ggsave(filename = "caenorhabditis_elegans.wagos_tmm_ma.png", plot = p_tmm, device = "png", 
#       path = outPath_6, 
#       dpi = 300, width = 18, height = 25)

norm_miR_simp_list = norm_miR_simp_list[!grepl("alg|^rde", names(norm_miR_simp_list))]

p_mir = ggarrange(plotlist = norm_miR_simp_list, ncol = 3, nrow = 7, common.legend = T, legend = "bottom")

ggsave(filename = "caenorhabditis_elegans.wagos_rRNA_ma.png", plot = p_mir, device = "png", 
       path = outPath_6, 
       dpi = 300, width = 18, height = 25)
