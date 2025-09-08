pacman::p_load(Biostrings)
args = commandArgs(trailingOnly = T)
faFile = args[1]
outFa = paste0(sub("\\..*", "", faFile), ".ltr_cons.fa")
fa = readDNAStringSet(faFile)
ltr = fa[grepl("LTR", names(fa)),]
writeXStringSet(ltr, outFa)
