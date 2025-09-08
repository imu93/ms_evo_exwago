## Get the ranges covered by exons, introns, and their anti-sense versions
pacman::p_load(rtracklayer, Rsamtools)
args = commandArgs(trailingOnly = T)
gffFile = args[1] # gene annotation 
faFile  = args[2] # genome file
outFile = sub(".fa", ".exon_intron_str.gff3", faFile)             
outFile2 = sub(".fa", ".exon_intron_unstr.gff3", faFile) 
# Using FaIndex
gff = import(gffFile)
sInfo = seqinfo(scanFaIndex(faFile))
seqlevels(gff) = seqlevels(sInfo)
seqinfo(gff) = sInfo

genomeLength = sum(seqlengths(gff))

# For Georgios/Cesare version
codingExons = gff[gff$type %in% c("five_prime_UTR", "CDS", "three_prime_UTR", "three_prime_utr", "five_prime_utr")]
# Since introns have to be defined in a reduced manner, do reduce for exons as well
codingExons = reduce(codingExons)
print(paste("All codingExons, percent =",round(sum(width(codingExons)) / genomeLength * 100,1)))

# Need to define own introns
# The following is very fast, but looses introns that overlap with exons of another transcript
allGaps = gaps(codingExons)
introns = allGaps[allGaps %within% gff[gff$type %in% c("mRNA","transcript")]]
introns = reduce(introns) # not required, but just in case
introns = introns[strand(introns) != "*"] # one gene was exactly the contig, thus also %within% '*'
print(paste("All introns, percent =",round(sum(width(introns)) / genomeLength * 100,1)))

# Define some AS categories
intronsAS = introns
strand(intronsAS) = ifelse(strand(introns) == "+","-","+")
codingExonsAS = codingExons
strand(codingExonsAS) = ifelse(strand(codingExons) == "+","-","+")

# Save the basic exon/intron GR objects for further use
intronExonList = list("exons_S"=codingExons, "exons_As"=codingExonsAS, "introns_S"=introns, "introns_As"=intronsAS)
# Merge, format and save
intexnList = list()
for (i in names(intronExonList)) {
  cat = intronExonList[[i]]
  cat$source = "rtracklayer"
  cat$type = i
  cat$score = NA
  cat$phase = NA
  cat$class = i
  intexnList[[i]] = cat
}
names(intexnList) = NULL
intexnstr = do.call(c, intexnList)
# export
export(intexnstr, outFile, "gff3")

# I need a version with out strand info in class
exIn = list("exons"=codingExons, "introns"=introns)
exInList = list()
for (i in names(exIn)) {
  cat = exIn[[i]]
  cat$source = "rtracklayer"
  cat$type = i
  cat$score = NA
  cat$phase = NA
  cat$class = i
  exInList[[i]] = cat
}
names(exInList) = NULL
exIn = do.call(c, exInList)
# export
export(exIn, outFile2, "gff3")
