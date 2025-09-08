#setwd("/home/isaac/storage/Data/ms_evo_exwago/analyses/soloLTR/ltr_ann/heligmosomoides_polygyrus/annotation")
pacman::p_load(Biostrings, rtracklayer, BSgenome, dplyr, stringr, pwalign)
args = commandArgs(trailingOnly = T)
#faFile = "heligmosomoides_polygyrus.ltr_cons.fa"
faFile = args[1]
outFile = paste0(sub("\\..*", "", faFile), ".ltr_library.fa")
selfPath = "self_blast/"
files = list.files(path = selfPath, pattern = ".*.self-blast.pairs.txt")
files_full =  paste0(selfPath, files)
names(files_full) = files
fa = readDNAStringSet(faFile)
names(fa) = sub(" .*", "", names(fa))
noLTR = DNAStringSet()
n_mods = DNAStringSet()
pid_solos = list()
for (model in names(files_full)) {
  #model = "rnd-4_family-1093_LTR.fa.self-blast.pairs.txt"
  tmp.mod = files_full[[model]]
  tmp.self = read.delim(tmp.mod, skip = 1, header = F, sep = " ")
  colnames(tmp.self) = c("ID", "qstart", "qend", "sstart", "send")
  tmp.self$ID = sub("__","#", tmp.self$ID)
  tmp.fa = fa[unique(tmp.self$ID),]
  tmp.gr_fa = GRanges(seqnames = names(tmp.fa), IRanges(start = 1, width = width(tmp.fa)), strand = "+") 
  # if I can not detect LTRs 
  if (nrow(tmp.self) < 2) {
    noLTR = c(noLTR, tmp.fa)
    next
  }
  # if this element has tirs instead of LTRs
  if (tmp.self[2,]$sstart > tmp.self[2,]$send) {
    noLTR = c(noLTR, tmp.fa)
  }
  # if I found LTRs
  else{
    # if present I need to remove directed repeats in the midle of the model 
    # First deffine model length
    tmp.self$mod_len = width(tmp.fa)
    # 
    tmp.self = tmp.self[2:nrow(tmp.self),]
    tmp.self = tmp.self[!tmp.self$qstart == tmp.self$send,]
    tmp.self$middle = tmp.self$mod_len / 2
    tmp.self = tmp.self[!(tmp.self$qstart >= (tmp.self$middle - 500) & tmp.self$qstart <= (tmp.self$middle + 500)), ]
    tmp.self = tmp.self[!(tmp.self$qend >= (tmp.self$middle - 500) & tmp.self$qend <= (tmp.self$middle + 500)), ]
    tmp.self$distance = ifelse(tmp.self$qend > tmp.self$sstart, abs(tmp.self$qstart-tmp.self$send), abs(tmp.self$sstart-tmp.self$qend))
    tmp.self = tmp.self[tmp.self$distance >= (tmp.self$mod_len / 2), ]
    if (nrow(tmp.self) <= 1) {
      noLTR = c(noLTR, tmp.fa)
      next
    }
    
    tosolo.gr = tmp.self
    gr1 = GRanges(seqnames = tosolo.gr$ID, IRanges(start = tosolo.gr$qstart, end = tosolo.gr$qend), strand = "+") 
    gr1 = reduce(gr1)
    if (length(gr1) > 2) {
      gr1 = reduce(gr1, min.gapwidth = 50)
    }
    if (length(gr1) > 2) {
      gr1 = gr1[order(seqnames(gr1), start(gr1)),]
      gr1 = gr1[c(1, length(gr1)),]
    }
    
    solo_gr_all = c(tmp.gr_fa, gr1) 
    # cleaning step to remove bases out of the LTR
    solo_gr = disjoin(solo_gr_all)
    solo_gr = solo_gr[rev(order(width(solo_gr))),]
    solo_gr = c(gr1, solo_gr[1,])
    solo_gr = solo_gr[order(start(solo_gr)),]
    solo_gr$type = c("LTR-5","CDS", "LTR-3")
    solo_fa =  getSeq(tmp.fa, solo_gr)
    names(solo_fa) = rep(names(tmp.fa), length(solo_gr)) %>%  str_replace("\\#", paste0("_",solo_gr$type, "#"))
    tst.piden = pairwiseAlignment(pattern = solo_fa[1,], subject = solo_fa[3]) %>% pid()
    if (tst.piden < 70) {
      noLTR = c(noLTR, tmp.fa)
    }
    else{
      solo_fa = solo_fa[rev(order(width(solo_fa))),][c(1:2),]
      names(solo_fa) = sub("_LTR-.*#","_LTR#", names(solo_fa))
      n_mods = c(n_mods, solo_fa)
    }
    
    pid_solos[[model]] = tst.piden
  }
}
# Save boxplot of LTR region piden
pdf("./pdf/ltr_region_piden_edges.pdf", width = 6, height = 7)
pid_solos %>% unlist() %>% boxplot(main="LTR region piden", ylab="% identity")
dev.off()

# now only the LTR region
soloLTR = n_mods[grepl("_LTR#", names(n_mods)),]
cdsLTR = n_mods[!grepl("_LTR#", names(n_mods)),]
# some controls
ltr_lib = fa[grepl("LTR", names(fa)),]
ltr_lib = ltr_lib[grepl("rnd", names(ltr_lib)),]
#Missed LTRs in the analyses 
missed_LTR = setdiff(sub("#.*","", names(ltr_lib)), c(cdsLTR, noLTR) %>% names() %>% str_replace("_CDS.*|#.*", ""))
# these repeats are fragments of LTR so does not matter
mLTR = ltr_lib[match(missed_LTR,sub("#.*","", names(ltr_lib)))]
all_LTR = c(soloLTR, cdsLTR, noLTR, mLTR)
# do I have all the LTRs
setdiff(sub("#.*","", names(ltr_lib)), sub("_CDS.*|_LTR.*|#.*", "", names(all_LTR)) %>% unique())

# As a test I wont use the whole annotation
# I'll do a first test with repeat masker
writeXStringSet(all_LTR, outFile)