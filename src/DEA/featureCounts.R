pacman::p_load(Rsubread, rtracklayer)
args = commandArgs(trailingOnly = T)
annFile = args[1]
outFile1 = sub("\\..*",".srna_1827.FC.Rds", annFile)
outFile2 = sub("\\..*", ".counts.1827.txt", annFile)
genome =  import(annFile)
table(genome$class)

# Let's build the SAF fromat
df = data.frame("GeneID"=genome$ID, "Chr"=seqnames(genome), "Start"=start(genome), "End"=end(genome), "Strand"=strand(genome))
bam_list = list.files(pattern = ".*mapped.bam$")
# And count
tb = featureCounts(files = bam_list, annot.ext = df, allowMultiOverlap= TRUE, strandSpecific = 1,
                                      fracOverlap= .7, nthreads = 10, useMetaFeatures=TRUE,
                   largestOverlap = TRUE, reportReadsPath = "./", reportReads = "CORE")

# allowMultiOverlap = TRUE means that I will allow feature-feature overlaps
# strandSpecific = 1 means that the count should be guided by strand
# fracOverlap means that at least 70% of each read must align to afeature to be assigned
# largestOverlap = TRUE specifiy that if I have reads aligning to both features, reads will be asigened to the one with the larges overlap
# Bearing in mid that all these bams lack the 'NH' tag all of them will count as unique mappers
# reportReads = info about which reads were used to count
# Finally I'll save the table

saveRDS(tb, outFile1)
write.table(tb$counts, outFile2, quote = F, sep = "\t", col.names = T, row.names = T)
