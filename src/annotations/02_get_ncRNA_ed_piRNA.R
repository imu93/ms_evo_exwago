setwd("~/storage/Data/evo_exwago/analyses/annotations/annotations_28052024/nippostrongylus_brasiliensis_28052024/ncrna/")
pacman::p_load(rtracklayer, Rsamtools, R.utils)
list.files()
miFile = list.files(pattern = ".*.1e-5.90piden.bed")
mdFile = list.files(pattern = ".*.mirDeep.bed")
rRFile = list.files(pattern =".*rRNAs.gff")
tRFile = list.files(pattern =".*tRNAscan.gff3")
rfFile = list.files(pattern =".*.infernal.gff3")
piFile = list.files(pattern =".*piRNAs.gff3")
lncFile = list.files(pattern = ".*lncRNAs.gff")
#yrFile = list.files(pattern =".*yRNA.1e-10.98piden.bed")
faFile = list.files(pattern =".*.fa$")
##################################################################################################################################################################
outFile = sub(".fa", ".ncRNA_unstr.gff3", faFile)
outFile2 = sub(".fa", ".ncRNA_str.gff3", faFile)

ncFile = list(miFile, mdFile, piFile, tRFile, rfFile, lncFile)
names(ncFile) = c("mirnas", "mdeep","pirnas", "trnas", "rfam", "lncrnas")

ncList = lapply(ncFile, function(x){
  tmp.gr = import(x)
  tmp.gr = GenomicRanges::trim(tmp.gr)
  return(tmp.gr)
  })

##################################################################################################################################################################
sInfo = seqinfo(scanFaIndex(faFile))
all_ncrna = list()
for (fam in names(ncList)) {
  x = ncList[[fam]]
  seqlevels(x) = seqlevels(sInfo)
  seqinfo(x) = sInfo
  x = GenomicRanges::trim(x)
  all_ncrna[[fam]] = x 
}
##################################################################################################################################################################
# Let's start with miRNAs for mirBase
all_ncrna$mirnas$class = "miRNA"
# Now mirDeep2
all_ncrna$mdeep$class = "miRNA"
all_ncrna$mdeep = all_ncrna$mdeep[all_ncrna$mdeep$score > 1,]
##################################################################################################################################################################
# I'll keep miRNA separately
all_ncrna$rfam$class = all_ncrna$rfam$type
all_ncrna$rfam$class = sub("mir.*", "miRNA", all_ncrna$rfam$class)
##################################################################################################################################################################
# Now format rfam
all_ncrna$rfam = all_ncrna$rfam[!grepl("rRNA|tRNA|Histone", all_ncrna$rfam$type),]
#simplify spliceosomal and yRNAs
all_ncrna$rfam$class =  ifelse(grepl("SmY|ceN72-3_ceN74-2", all_ncrna$rfam$type), "yRNA", all_ncrna$rfam$class)
all_ncrna$rfam$class = ifelse(grepl("U[0-9]",all_ncrna$rfam$type), "snRNA", all_ncrna$rfam$class)
all_ncrna$rfam$class = ifelse(grepl("Afu|sn[0-9]|SNOR|ceN",all_ncrna$rfam$type), "snoRNA", all_ncrna$rfam$class)
all_ncrna$rfam$class = ifelse(!grepl("yRNA|snoRNA|snRNA|miRNA", all_ncrna$rfam$class), "other_ncRNA", all_ncrna$rfam$class)
all_ncrna$rfam = all_ncrna$rfam[!grepl("archaea|microsporidia|bacteria", all_ncrna$rfam$type),]
table(all_ncrna$rfam$class)
##################################################################################################################################################################
# Let's process rRNA
#all_ncrna$rrnas$class=  "rRNA"
rnamer_df = read.delim(rRFile, comment.char = "#", header = F)
colnames(rnamer_df) = c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "class", "name")
all_ncrna$rrnas = makeGRangesFromDataFrame(rnamer_df, keep.extra.columns = T)
##################################################################################################################################################################
# The next calss will be tRNA
all_ncrna$trnas$class = "tRNA"
##################################################################################################################################################################
# Finally piRNAs
all_ncrna$pirnas$class = "piRNA"
##################################################################################################################################################################
# Do not miss yRNAs from nar 2019
# Now I'll edit the yRNA file
#all_ncrna$yrnas$class = "yRNA"
##################################################################################################################################################################
# Finally lncRNAs
all_ncrna$lncrnas$class = all_ncrna$lncrnas$type
##################################################################################################################################################################
all_ncrna_list = list()
for (fam in names(all_ncrna)) {
  x = all_ncrna[[fam]] 
  x = x[,"class"]
  all_ncrna_list[[fam]] = x 
}
names(all_ncrna_list) = NULL
all_ncrna = do.call(c, all_ncrna_list)
all_ncrna = GenomicRanges::trim(all_ncrna)
##################################################################################################################################################################
# Merge, reduce and estimate nt
all_ncrna$type = all_ncrna$class
all_ncrna$class = NULL
all_ncList = split(all_ncrna,all_ncrna$type)
all_ncWidth = sapply(all_ncList,function(x)sum(width(x)))

all_ncRed <- lapply(all_ncList, reduce)
all_ncRedWidth = sapply(all_ncRed, function(x) sum(width(x)))
res = data.frame(before=all_ncWidth,reduced=all_ncRedWidth)
res

nclst=GRangesList()
for (i in names(all_ncRed)) {
  x = all_ncRed[[i]]
  x$type = i
  nclst[[i]]=x
}
##################################################################################################################################################################
# Finally I'm gonna create the the GRanges object
names(nclst) = NULL
ncRNAs = do.call(c, nclst)
ncRNAs = ncRNAs[order(seqnames(ncRNAs), start(ncRNAs)),]
##################################################################################################################################################################
# Add sinfo and save
sInfo = seqinfo(scanFaIndex(faFile))
seqlevels(ncRNAs) = seqlevels(sInfo)
ncRNAs = GenomicRanges::trim(ncRNAs)
seqinfo(ncRNAs) = sInfo
export(ncRNAs, outFile, "gff3")

ncRNAs$type = paste0(ncRNAs$type, "_S")
nclst = split(ncRNAs, ncRNAs$type)
ncAS <- list()
for (i in names(nclst)) {
  print(i)
  ncRNA_As <- nclst[[i]]
  strand(ncRNA_As) <- ifelse(strand(ncRNA_As) == "+","-","+")
  ncRNA_As$type = gsub("_S", "_As", ncRNA_As$type)
  ncAS[[i]] = ncRNA_As
}
names(ncAS) = NULL
ncAs = do.call(c, ncAS)
ncRNAs = c(ncRNAs, ncAs)
##################################################################################################################################################################
ncRNAs = ncRNAs[order(seqnames(ncRNAs), start(ncRNAs)),]
# Add sinfo and save
sInfo = seqinfo(scanFaIndex(faFile))
seqlevels(ncRNAs) = seqlevels(sInfo)
ncRNAs = GenomicRanges::trim(ncRNAs)
seqinfo(ncRNAs) = sInfo
export(ncRNAs, outFile2, "gff3")

