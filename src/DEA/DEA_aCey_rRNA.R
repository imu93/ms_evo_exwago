pacman::p_load(edgeR, ggplot2, ggpubr, ggplotify, dplyr, purrr, gplots, 
               ggExtra, ggforce, rtracklayer, Biostrings, BSgenome)
list.files()

# Define Paths and inFiles
inPath_1 = "~/storage/Data/ms_evo_exwago/raw/ShortStack/aCey/"
inFile_1 = "ancylostoma_ceylanicum.counts.1827.txt" # main table


outPath_1 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/tmm/"
outPath_2 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/rrna_norm/"
outPath_3 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/dge/"
outPath_4 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/camera/class/" 
outPath_5 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/camera/family/"
outPath_6 = "~/storage/Data/ms_evo_exwago/analysis/figures/"
outPath_7 = "~/storage/Data/ms_evo_exwago/analysis/de_tables/dge_tmm/"

# Read table and edit colnames
counts = read.delim(paste0(inPath_1, inFile_1))
colnames(counts) = sub("\\.trim.*", "", colnames(counts))
colnames(counts)

# Select IP and Adult total
counts = counts[,grepl("Acey_IP_polyP|Acey_adult_polyP", colnames(counts))]


colnames(counts)
# Define groups
groups = colnames(counts) %>% strsplit("_") %>%
  map(function(x){paste(x[c(1,2,3)], collapse = "_")}) %>%  unlist() %>% 
  factor()

# Estimate libsize per replicate and build DGE object
lib.size = colSums(counts)
dge = DGEList(counts=counts, group=groups, lib.size = lib.size)

# Filter by expression
keep = filterByExpr(dge, min.count= 5)
dge = dge[keep, , keep.lib.sizes=FALSE]
print(table(keep))
f.exp.counts = colSums(dge$counts)

# Now I'll start using TMM to estimate normFactors
dge = calcNormFactors(dge, method = "TMM")

# Now I gonna produce a MDS to check of the samples group as expected
colors = rich.colors(length(levels(groups)), alpha = .5)
names(colors) = levels(groups)
par(mfrow=c(1,2))
plotMDS(dge , col=colors[dge$samples$group], 
        pch = 16, cex=3,main="MDS plot of A. ceylanicum exWAGO vs Adult")
plot.new()
legend("left", legend = names(colors), pch = 16, col = unique(colors), cex=2)


# Create the design object
design = model.matrix(~0+dge$samples$group)
colnames(design) =  c(levels(dge$samples$group))

# Let's estimate de dispersion
dge = estimateDisp(dge, design=design, robust=TRUE, tagwise = T)
plotBCV(dge)
plot(dge$AveLogCPM,dge$trended.dispersion)

#fit = glmQLFit(dge,design = dge$design)
#plotQLDisp(fit)

# Let's estimate de dispersion
dge = estimateDisp(dge, design=design, robust=TRUE)
plotBCV(dge)
plot(dge$AveLogCPM,dge$trended.dispersion)
fit = glmFit(dge,design = dge$design,dispersion = dge$common.dispersion)

# Now I'm going to specify the contrast and perform the DEA
contrast = makeContrasts("Acey_IP_polyP-Acey_adult_polyP",levels = dge$design)
gt = glmTreat(fit, contrast=contrast, lfc = 1)
topTags(gt)
dt = decideTestsDGE(gt, adjust.method="BH", p.value=0.05, lfc = 1)
table(dt)
topTable = topTags(gt, n=Inf)$table
topTable = topTable[rownames(dt),]

write.table(topTable, paste0(outPath_1, "ancylostoma_ceylanicum.exwago_ip_vs_total_fdr05_lfc1_tmm_topTbl.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)

# Simplify table to plot 
de_tbl = data.frame('Ave_CPM'=topTable$logCPM, 'log-FC'= topTable$logFC, 'cl'= dt)
colnames(de_tbl) = c('Ave_CPM', 'logFC', 'cl')
de_tbl$col = ifelse(de_tbl$cl == 0, "#C0C0C0", 
                    ifelse(de_tbl$cl == 1, "#5E2129", "#333333"))
de_tbl$ord = ifelse(de_tbl$col == "#C0C0C0",
                    1, ifelse(de_tbl$col == "#333333", 2, 3))
de_tbl = de_tbl[order(de_tbl$ord),]
de_tbl[de_tbl$logFC <= -1,]
de_tbl$col = factor(de_tbl$col, levels = c("#5E2129", "#C0C0C0", "#333333"))
de_tbl$exp = "A. ceylanicum IP vs Adult total"
de_tbl$exp = factor(de_tbl$exp)
levels(de_tbl$exp) = c(expression(paste(italic("A. ceylanicum "),
                                        "exWAGO IP vs Adult total")))

write.table(de_tbl ,paste0(outPath_1, "ancylostoma_ceylanicum.exwago_ip_vs_total_fdr05_lfc1_tmm_deTbl.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)

p = ggplot(de_tbl, aes(x=Ave_CPM, y=logFC, colour = col)) +
  geom_point() +  xlab(bquote(log[2]~' CPM')) + theme_test() +
  ylab(bquote(~log[2]~ 'FC IP/Total')) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  geom_hline(aes(yintercept = 0), colour = "blue",
             linewidth = 1, linetype="dashed", alpha=0.4) +
  facet_wrap(~ exp, labeller = label_parsed) +
  scale_color_manual(values = adjustcolor(c("#155AA3","#C0C0C0", "#333333"), alpha.f = .6), 
                     labels=c('IP', 'not-DE','Unbound'), name=NULL) +
  theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
  theme( strip.text = element_text(size = 18, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=10))) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches"))
p
saveRDS(dge,paste0(outPath_7,"ancylostoma_ceylanicum.exwago_ip_adult_tmm_dge.Rds"))
# As expected TMM is unbalanced, thus I have lots of type II errors in DEA
#############################################################################################################################################################
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
ggsave("ancylostoma_ceylanicum.tmm_normcats.png", device = png, plot = cts_1, width = 14, height = 5, dpi = 300,
       path = "/home/isaac/storage/Data/ms_evo_exwago/analysis/de_tables/figures")

write.table(df, paste0("~/storage/Data/ms_evo_exwago/analysis/de_tables/tmm_norm_cats/",
  "aceylanicum.tmm_nomcats.txt"), sep = "\t", quote = F, col.names = T, row.names = T)

#############################################################################################################################################################

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
counts = counts[,grepl("Acey_IP_polyP|Acey_adult_polyP", colnames(counts))]

# Get the rigth columns
counts = counts[,grepl("Acey_IP_polyP|Acey_adult_polyP", colnames(counts))]
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
contrast = makeContrasts("Acey_IP_polyP-Acey_adult_polyP",levels = dge$design)
gt = glmTreat(fit,  contrast=contrast, lfc = 1)
topTags(gt)
dt = decideTestsDGE(gt, adjust.method="BH", p.value=0.01, lfc = 1)
normrrna_dt = table(dt)
normrrna_dt
topTable = topTags(gt, n=Inf)$table
topTable = topTable[rownames(dt),]

outToptable = paste0(outPath_2, "ancylostoma_ceylanicum.exwago_IP_adult_lfc1_fdr01_rRNA_topTbl.txt")
write.table(topTable,outToptable, quote = F, sep = "\t")
####################################################################################################################################################################
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
de_tbl_norm$exp = "A. ceylanicum IP vs Adult total"
de_tbl_norm$exp <- factor(de_tbl_norm$exp)
ncont = "IP vs Adult total"
levels(de_tbl_norm$exp) <- c(expression(paste(italic("A. ceylanicum "))))


outDEtbl = paste0(outPath_2, "ancylostoma_ceylanicum.exwago_IP_adult_lfc1_fdr01_rRNA_deTbl.txt")

write.table(de_tbl_norm,outDEtbl, quote = F, sep = "\t")

ma = ggplot(de_tbl_norm, aes(x=Ave_CPM, y=logFC, colour = col)) +
  geom_point() +  xlab(bquote(~log[2]~' CPM')) + theme_test() +
  ylab(bquote(~log[2]~ 'FC IP/Total')) + 
  theme(axis.title.x = element_text(face = "bold", size = 20), 
        axis.text.x = element_text(face = "bold", size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20),
        axis.text.y = element_text(face = "bold", size = 20)) +
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
  facet_wrap(~ exp, labeller = label_parsed) +
  theme( strip.text = element_text(size = 18, face = "bold", colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) 
ma 
#############################################################################################################################################################
cats = c("rRNA_S", "miRNA_S", "tRNA_S")
names(cats) = c("rRNA_S", "miRNA_S", "tRNA_S")
df_lst = list()
for (i in names(cats)) {
  df = de_tbl_norm
  df$col = ifelse(grepl(i, rownames(df)), "#5E2129", "#C0C0C0")
  df$ord = ifelse(df$col == "#C0C0C0", 1,2)
  df = df[order(df$ord),]
  df$bty = i
  df_lst[[i]] = df
}
df = do.call(rbind, df_lst)
cts_2 = ggplot(df, aes(x=Ave_CPM, y=logFC, colour = col)) +
  geom_point(show.legend = F) + facet_wrap(~ bty) +
  xlab(bquote(log[2]~' CPM')) + theme_test() +
  ylab(bquote(~log[2]~ 'FC IP/Total')) +
  theme(axis.title.x = element_text(face = "bold", size = 20),
        axis.text.x = element_text(face = "bold", size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20),
        axis.text.y = element_text(face = "bold", size = 20)) +
  geom_hline(aes(yintercept = 0), colour = "blue",
             linewidth = 1, linetype="dashed", alpha=0.4) +
  scale_color_manual(values = adjustcolor(c("#155AA3","#C0C0C0"), alpha.f = .6)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
  theme( strip.text = element_text(size = 18, face = "bold", colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=10))) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches"))


ggarrange(ma, cts_2, ncol = 1, heights =c(1,.7))

saveRDS(dge,paste0(outPath_3, "ancylostoma_ceylanicum.ip_vs_adult_rRNA_norm_dge.Rds"))


#de_tbl_norm[grepl("LTR", rownames(de_tbl_norm)) & de_tbl_norm$cl == 1,]
####################################################################################################################################################################
# Let's run CAMERA

# In this section I will run CAMERA at the class and family levels
# Let's start defining classes
# This section may not 100% necessary

categories = c("DNA.*As|RC.*As", "DNA.*S|RC.*S", "LINE.*As|PLE.*As|Retrop*As",
               "LINE.*S|PLE.*S|Retrop*S", "LTR.*As", "LTR.*S",
               "Unk", "Low|Simple|Sat|SINE|MITE", "lincRNA_As", 
               "miRNA_S", "_rRNA_S", "^tRNA_S",
               "_rRNA_As|yRNA|snRNA|snoRNA|other_ncRNA|piRNA|lncRNA|miRNA_As|lincRNA_S|^tRNA_As", 
               "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")

names(categories) = c("DNA_As", "DNA_S", "LINE_As", "LINE_S", "LTR_As",
                      "LTR_S", "Unknown", "Other_repeat", "lincRNA_As", 
                      "miRNA_S", "rRNA_S", "tRNA_S",  "other_ncRNA", "exons_As",
                      "exons_S", "introns_As", "introns_S", "intergenic")

####################################################################################################################################################################
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

write.table(cameraRes, 
            paste0(outPath_4,"ancylostoma_ceylanicum.IP_camera_enr_class_simple_rRNA.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

####################################################################################################################################################################
# Now camera but using families 
# The following lines are just for split the table by families
aveCPM = 2 ** aveLogCPM(dge, normalized.lib.sizes = T, prior.count = .01)
names(aveCPM) = rownames(dge$counts)

aveCPMByfam = split(aveCPM, sub(":.*", "", names(aveCPM)))

categories =  names(aveCPMByfam)[grepl("DNA|LINE|LTR|Unk",names(aveCPMByfam))]
names(categories) = categories

allFams = list()
for (i in names(categories)) {
  tmp.cat = categories[i]
  tmp.sqs = dge$counts[grepl(tmp.cat, rownames(dge$counts)),] %>% rownames()
  allFams[[i]] = tmp.sqs 
}

allFamsIdx = lapply(allFams, function(x) rownames(dge$counts) %in% x)
allFamsIdx = allFamsIdx[sapply(allFamsIdx, sum) > 5]
cameraRes = camera(dge, allFamsIdx, design, contrast)
cameraRes = cameraRes[match(names(categories), rownames(cameraRes)),]
rownames(cameraRes) = names(categories)
cameraRes$Feature = ifelse(grepl("DNA|RC", rownames(cameraRes)), "DNA", 
                           ifelse(grepl("LINE|PLE", rownames(cameraRes)), "LINE",
                                  ifelse(grepl("LTR", rownames(cameraRes)), "LTR", "Unknown")))
cameraRes = cameraRes[!grepl("^NA", rownames(cameraRes)),]

write.table(cameraRes, 
            paste0(outPath_5,"ancylostoma_ceylanicum.IP_camera_enr_family_simple_rRNA.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)
