pacman::p_load(rtracklayer, Rsamtools, stringr, dplyr)
list.files()
args = commandArgs(trailingOnly = T)
# Load Files
soloLTR = import(args[1]) # soloLTR annotation
flLTR = import(args[2]) # Full length annotation
allLTR = import(args[3]) # transfered annotation
de_lts = read.delim(args[4]) # DE table
outfile = paste0(sub("\\..*", "", args[1]), ".exWAGO_DE_LTRsplit.gff3")
# First I ned to obtain LTRs that are not FL or solo 
annLTRs = c(flLTR, soloLTR)
LTR_frag = allLTR[!allLTR %over% annLTRs,]

print(paste0("overlaps solo vs other:", subsetByOverlaps(soloLTR, LTR_frag) %>% length()))
print(paste0("overlaps fl vs other:", subsetByOverlaps(flLTR, LTR_frag) %>% length()))

# Mege all LTRs
allLTR =  c(soloLTR, flLTR, LTR_frag)
table(allLTR$type)


# Now I need only the sgements of DE LTR
de_ltr = de_lts[grepl("LTR", rownames(de_lts)) & de_lts$cl == 1,]
DELTR = allLTR[allLTR$segment_ID %in% rownames(de_ltr),] 
#boxplot(DELTR$prec_div %>% as.numeric())
DELTR$ID =  DELTR$ID %>% sub("_NA:", "", .)

export(DELTR, outfile)

