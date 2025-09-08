pacman::p_load(rtracklayer, plyranges, stringr)
list.files()
# With this peace of script I'll retriever the repetitive model of each range after having used the reduce function
args = commandArgs(trailingOnly = T)
annFile= args[1] #  segment annotation from 04_segmented_annotation
repsFile = args[2] # merged repeats from 01_assign_overlaps
outFile1 = sub("segment.gff3", "ann2prop.Rds", annFile)
outFile2 = sub("segment.gff3", "segment_final.gff3", annFile)
#################################################################################################################################################################
genome = import(annFile)
genome$ID = paste(genome$class, 1:length(genome$class), sep = ":")
rawReps = import(repsFile)
reps = genome[grepl("DNA|RC|LINE|SINE|LTR|Retro|Unk|Simple|Low", genome$class),]
table(reps$class)
#################################################################################################################################################################
# First I'll get the intersect between the rawReps and the reduced once
srtandedReps = reps
ovs = plyranges::join_overlap_intersect(rawReps, srtandedReps)
ovs = ovs[order(seqnames(ovs), start(ovs)),]
ovs$wdt = width(ovs)
ovs$rep.wdt = width(reps[match(ovs$ID.y, reps$ID),])
ovs$cov = ovs$wdt*100/ovs$rep.wdt
ovs = ovs[as.character(ovs$type.x) == as.character(ovs$type.y),]
# I need to extract all the elements with more than one element
dups = ovs[duplicated(ovs$ID.y) | duplicated(ovs$ID.y, fromLast=TRUE),]
# Also I'll save the elements with just one model to use them again in the last step of the script
ovs_ed = ovs[!ovs$ID.y %in% dups$ID,]
# Now I'll group all the elements by ID and retriever those with the max cov value
dups_red = dups[dups$cov == ave(dups$cov, dups$ID.y, FUN=max),]
# Some of them are duplicated but with the same maximum value so I'll chose randomly
dups_red = dups_red[!duplicated(dups_red$ID.y)]
# Finally I'll merge the objects to give a model to each repeat
ovs_ed= c(dups_red, ovs_ed)
srtandedReps$rep_name = ovs_ed[match(srtandedReps$ID, ovs_ed$ID.y),]$rep_name
#srtandedReps$rep_name = paste(srtandedReps$rep_name, collapse = ",")
genome = genome[!grepl("DNA|RC|LINE|SINE|LTR|Unk|Retro|Simple|Low", genome$class),]
genome = c(genome, srtandedReps)
saveRDS(genome, outFile1)
#################################################################################################################################################################
uns= genome[strand(genome) == "*",]
strand(uns) = "+"
uns_As = uns
strand(uns_As) = "-"
uns = c(uns, uns_As)
genome = genome[!strand(genome) == "*",]
genome = c(uns, genome)
genome = genome[order(seqnames(genome), start(genome)),]
genome$source = NULL
genome$score = NULL
genome$phase = NULL
#################################################################################################################################################################
genome$ID =  ifelse(grepl("DNA|LINE|SINE|LTR|PLE", genome$class), str_replace(genome$ID, "(_As:.*$|_S:.*)", paste("_", genome$rep_name,"\\1", sep = "")), genome$ID)
genome$ID =  ifelse(grepl("Unknown", genome$class), str_replace(genome$ID, "(:.*$)", paste("_", genome$rep_name,"\\1", sep = "")), genome$ID)
genome$igv =paste(seqnames(genome), ":", start(genome), "-", end(genome), sep="")
sum(width(genome))
sum(width(reduce(genome)))
export(genome, outFile2, "gff3")
