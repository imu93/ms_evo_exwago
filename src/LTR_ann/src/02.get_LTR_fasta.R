# This script aims to retrieve annotated segments as LTR from the initial
# segmented annotations

pacman::p_load(rtracklayer, BSgenome)
args = commandArgs(trailingOnly = T)
inFa = args[1]
inSeg = args[2]

# Check points
if(is.na(inFa) == T){
  stop("No genome assembly was provided")
}

if(is.na(inSeg) == T){
  stop("No genome annotation was provided")
}

# Parse input
fa = readDNAStringSet(inFa)
gff = import(inSeg)
outFile = paste0(sub("\\..*", "", inFa), ".genomic.LTR_segements_final.fa")
# Get all LTRs
ltrs = gff[grepl("LTR",gff$type),]
sqs = getSeq(fa, ltrs)
names(sqs) = ltrs$ID
# Write fasta file
writeXStringSet(sqs, outFile)