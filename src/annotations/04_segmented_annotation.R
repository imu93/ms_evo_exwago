pacman::p_load(rtracklayer, Rsamtools)
args = commandArgs(trailingOnly = T)

faFile = args[1] # fasta fiel
repFile = args[2] # repeat file stranded from script 01
codFile = args[3] # coding annotation stranded  from script 03
ncFile =  args[4] # ncrnas stranded from script 02

repList = import(repFile)
codList = import(codFile)
ncList  = import(ncFile)
outFile = sub("exon_intron_str.gff3", "segment.gff3", codFile)
outFile2 = sub("exon_intron_str.gff3","segments_mb.Rds", codFile)

sInfo = seqinfo(scanFaIndex(faFile))
seqlevels(repList) = seqlevels(sInfo)
seqinfo(repList) = sInfo

seqlevels(codList) = seqlevels(sInfo)
seqinfo(codList) = sInfo

seqlevels(ncList) = seqlevels(sInfo)
seqinfo(ncList) = sInfo

names(ncList) <- NULL

ncRNA = split(ncList, ncList$type)
names(ncRNA)
ncRNA = ncRNA[!grepl("pseudo", names(ncRNA)),]
# I will use the same order as Cei
ind.ncRNA <- c("miRNA_S","miRNA_As", "piRNA_S", "piRNA_As","yRNA_S", "yRNA_As", "tRNA_S","tRNA_As", "28s_rRNA_S",
               "28s_rRNA_As", "18s_rRNA_S", "18s_rRNA_As", "8s_rRNA_S", "8s_rRNA_As",
               "snRNA_S","snRNA_As","snoRNA_S","snoRNA_As", "lincRNA_S", "lincRNA_As", "lncRNA_S",
               "lncRNA_As","other_ncRNA_S","other_ncRNA_As")

setdiff(names(ncRNA), ind.ncRNA)
ncRNA = ncRNA[ind.ncRNA]
names(ncRNA)
# Let's process repeats
# Remove useless columns
repList$rep_name = NULL
repList$class = NULL
repList$score = NULL
repList$phase = NULL
repList$ID = NULL
repList$wdt = NULL
repList$class = repList$type
repList$class = as.character(repList$class)
repList$type = NULL
repList$source = NULL
table(repList$class)

# I will saparate unstranded elements
unstrandedList = repList[grepl("Low|Simple|Unk|SINE|MITE|Satellite", repList$class),]
unstrandedList = split(unstrandedList, unstrandedList$class)
uns_ord = rev(sort(unlist(lapply(unstrandedList, function(x){sum(width(reduce(x)))}))))
simLow = uns_ord[c("Satellite", "Simple_repeat", "Low_complexity", "Unknown")]
uns_ord = uns_ord[!grepl("Sate|Simple|Low|Unk", names(uns_ord))]
uns_ord = c(uns_ord, simLow)
unstrandedList = unstrandedList[names(uns_ord)]
names(unstrandedList)

r = repList
repList = repList[!grepl("Low|Simple|Unk|SINE|MITE|Sate", repList$class),]

# I will order repeats by their length in the genome
repList = split(repList, repList$class)
ord = rev(sort(unlist(lapply(repList, function(x){sum(width(reduce(x)))}))))
reord = sort(ord[grepl("Sat", names(ord))])
ord = ord[!grepl("Sat", names(ord))]
ord = c(ord, reord)
repList = repList[names(ord)]
names(repList)

# Before remove overlaps I will edit codList
coding = split(codList, codList$class)

# Define all in the order of priority
strandedList <- c(ncRNA, coding, repList)
names(strandedList)

reord = strandedList[c("exons_S", "exons_As", "introns_S", "introns_As")]
strandedList = strandedList[!grepl("exons_S|exons_As|introns_S|introns_As", names(strandedList))]
strandedList = c(strandedList, reord)
names(strandedList)

# Let's check that all elements remains in the list
stopifnot(all(c(names(repList), names(codList), names(ncRNA)) %in% c(names(strandedList), names(unstrandedList))))

# Calculate sum of all widths before overlap reduction
allWidth <- do.call("sum",sapply(c(strandedList, unstrandedList), width))

# Start with one, cycle through the others
print("Removing overlaps, according to sorted types")
allSegments <- strandedList[[1]]
allSegments$type <- names(strandedList)[1]
for (type in names(strandedList)[-1]) {
  print(type)
  trimmedSegs <- setdiff(strandedList[[type]], allSegments)
  if (length(trimmedSegs) > 0) {
    trimmedSegs$type <- type
    allSegments <- c(allSegments, trimmedSegs)
  }
}
for (type in names(unstrandedList)) {
  print(type)
  x =  unstrandedList[[type]]
  trimmedSegs <- setdiff(x, allSegments, ignore.strand=TRUE)
  if (length(trimmedSegs) > 0) {
    trimmedSegs$type <- type
    allSegments <- c(allSegments, trimmedSegs)
  }
}

# Calculate difference in sum of all widths after overlap reduction
remWidth <- round((allWidth-sum(width(allSegments)))/allWidth*100,2)
print(paste0("Removed ",remWidth,"% of annotated bases due to overlaps in sorted types"))

# Calculate all intergenic spaces
allGaps <- gaps(reduce(allSegments, ignore.strand=TRUE))
allGaps <- allGaps[strand(allGaps) == "*"]
allGaps$type <- "intergenic"

# Calculate the ranges that were removed due to any overlap
allSegments <- c(allSegments, allGaps)

# Check genome sizes
allTypes <- sapply(split(width(allSegments), allSegments$type), sum)

# multiply *2 those that are on the '*' strand
for (type in c(names(unstrandedList),"intergenic")) {
  allTypes[type] <- allTypes[type] * 2
}

allTypes = allTypes[!is.na(allTypes)]

print(paste("Annotated genome size=",sum(allTypes)/1e6/2,"Mb"))
print(paste("Real genome size=",sum(seqlengths(allGaps))/1e6,"Mb"))
ann_bases = allTypes/1e6/2
allSegments$class = NULL
allSegments$class = allSegments$type

export(allSegments, outFile, "gff3")
saveRDS(ann_bases, outFile2)
