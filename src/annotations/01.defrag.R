# pacman will install the packages
if (!require("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(rtracklayer, plyranges, dplyr, stringi)
args = commandArgs(trailingOnly = T)

# Let's read the repeatCraft file
repFile = args[1] # .rmerge.gff file from RepeatCraft
tblFile = args[2] # .tbl from RepeatMasker

##################################################################################################################################################
# Add some outputs
outFile1 = paste0(sub("\\..*", "", tblFile), ".summary.class_repeats.txt")
outFile2 = paste0(sub("\\..*", "", tblFile), ".summary.sFam_repeats.txt")
outFile3 = paste0(sub("\\..*", "", tblFile), ".merged_repeats_unstr.gff3")
outFile4 = paste0(sub("\\..*", "", tblFile), ".merged_repeats_str.gff3")
##################################################################################################################################################
# I'll use the tbl file to get genome size
tbl = readLines(tblFile)
reps =  import(repFile)
reps$wdth = width(reps)
# Remove non-relevant classes
reps = reps[!grepl("^rRNA|snRNA|^tRNA|ARTEFACT", reps$type),]
reps$type = as.character(reps$type)
# split and reduce by type
# This will remove the overlaps between elements of the same super family
repList = split(reps, reps$type)
unstrandedList = repList[grepl("Unk|Sat|Low|Simple|MITE|SINE", names(repList)),]
strandedList = repList[!grepl("Unk|Sat|Low|Simple|MITE|SINE", names(repList)),]

# Note: stranded repeats can overlap just if they have a different strand 
# unstranded repeats however, should not have overlap
unstrandedReduced = lapply(unstrandedList, function(x){reduce(x, with.revmap=T, ignore.strand=T)})
strandedReduced = lapply(strandedList, function(x){reduce(x, with.revmap=T, ignore.strand=F)})

# Now I'll sort by genomic coverage to get my hierarchy order
unstrandedOrd = lapply(unstrandedReduced, function(x){sum(width(x))}) %>%
  unlist() %>% sort() %>% rev()
reord = unstrandedOrd[grepl("Satellite|Simple_repeat|Low_complexity|Unknown", names(unstrandedOrd))]
unstrandedOrd = unstrandedOrd[! names(unstrandedOrd) %in% names(reord)]
unstrandedOrd = c(unstrandedOrd, reord)
unstrandedListRed = unstrandedReduced[names(unstrandedOrd)]
# I'll use the following order to assign overlaps between elements of different super families
print("Unstranded repeats have this hierarchy")
print(names(unstrandedListRed))

# Now the strandedList
strandedOrd = lapply(strandedReduced, function(x){sum(width(x))}) %>%
  unlist() %>% sort() %>% rev()
strandedListRed = strandedReduced[names(strandedOrd)]
print("Stranded repeats have this hierarchy")
print(names(strandedListRed))
##################################################################################################################################################
# Retrieve  rep_name with revmap based on the longest element of the overlap
strandedListRedcls = list()
for (i in names(strandedListRed)) {
  print(i)
  type = i
  fam = strandedListRed[[i]]
  revmap = mcols(strandedListRed[[i]])$revmap
  fam$type = extractList(mcols(strandedList[[i]])$type, revmap)
  fam$wdth = extractList(mcols(strandedList[[i]])$wdth, revmap)
  fam$rep_name = extractList(mcols(strandedList[[i]])$ID, revmap)
  uniques = fam[unlist(lapply(fam$revmap, length)) == 1,]
  dups = fam[unlist(lapply(fam$revmap, length)) > 1,]
  dups$type = unlist(Map('[', dups$type, which.max(dups$wdth)))
  dups$rep_name = unlist(Map('[', dups$rep_name, which.max(dups$wdth)))
  redFam = c(uniques, dups)
  redFam$revmap = NULL
  redFam$wdth = NULL
  redFam$type = as.character(redFam$type)
  redFam$rep_name = as.character(redFam$rep_name)
  strandedListRedcls[[i]] = redFam
}
# Same for unstranded
unstrandedListRedcls = list()
for (i in names(unstrandedListRed)) {
  print(i)
  type = i
  fam = unstrandedListRed[[i]]
  revmap = mcols(unstrandedListRed[[i]])$revmap
  fam$type = extractList(mcols(unstrandedList[[i]])$type, revmap)
  fam$wdth = extractList(mcols(unstrandedList[[i]])$wdth, revmap)
  fam$rep_name = extractList(mcols(unstrandedList[[i]])$ID, revmap)
  uniques = fam[unlist(lapply(fam$revmap, length)) == 1,]
  dups = fam[unlist(lapply(fam$revmap, length)) > 1,]
  dups$type = unlist(Map('[', dups$type, which.max(dups$wdth)))
  dups$rep_name = unlist(Map('[', dups$rep_name, which.max(dups$wdth)))
  redFam = c(uniques, dups)
  redFam$revmap = NULL
  redFam$wdth = NULL
  redFam$type = as.character(redFam$type)
  redFam$rep_name = as.character(redFam$rep_name)
  unstrandedListRedcls[[i]] = redFam
}
# To do not mess up previous versions I'll use the same variables 
strandedList = strandedListRedcls
unstrandedList = unstrandedListRedcls

##################################################################################################################################################
# The next step is to remove overlaps between super families
# I'll use findOverlaps to obtain the rep name after trim with setdiff
allWidth = do.call("sum",sapply(c(strandedList, unstrandedList), width))
print("Removing overlaps, according to sorted types")
allSegments = strandedList[[1]]
for (type in names(strandedList)[-1]) {
  print(type)
  trimmedSegs = GenomicRanges::setdiff(strandedList[[type]], allSegments)
  if (length(trimmedSegs) > 0) {
    trimmedSegs$type = type
    overs = findOverlaps(trimmedSegs, strandedList[[type]])
    trimmedSegs$rep_name = strandedList[[type]][subjectHits(overs)]$rep_name
    allSegments = c(allSegments, trimmedSegs)
  }
}
for (type in names(unstrandedList)) {
  print(type)
  fam =  unstrandedList[[type]]
  trimmedSegs <- GenomicRanges::setdiff(fam, allSegments, ignore.strand=TRUE)
  if (length(trimmedSegs) > 0) {
    trimmedSegs$type = type
    overs = findOverlaps(trimmedSegs, fam)
    trimmedSegs$rep_name = fam[subjectHits(overs)]$rep_name
    allSegments = c(allSegments, trimmedSegs)
  }
}

# Format allsegements
allSegments$source = "RepeatMasker"
allSegments$score = NA
allSegments$phase = NA
allSegments = data.frame(allSegments)

allSegments = allSegments[,c("seqnames", "start", "end", "strand", "source",
                             "type", "score", "phase", "rep_name")]
repsByTypeRed = makeGRangesFromDataFrame(allSegments, keep.extra.columns = T)

# Calculate difference in sum of all widths after overlap reduction
remWidth = round((allWidth-sum(width(repsByTypeRed)))/allWidth*100,2)
print(paste0("Removed ",remWidth,"% of annotated bases due to overlaps in sorted types"))

##################################################################################################################################################
# Now lets add the class feature
cats= c("DNA", "RC", "LINE", "SINE", "LTR", "Unknown", "Satellite")
repsByTypeRed$class =  "Other_repeat"
for (i in cats) {
  repsByTypeRed$class = ifelse(grepl(i, repsByTypeRed$type), i, repsByTypeRed$class)
}
# Do not forget PLE
repsByTypeRed$class = ifelse(grepl("PLE", repsByTypeRed$type), "LINE", repsByTypeRed$class)
# Give fromat
repsByTypeRed$class = factor(repsByTypeRed$class, levels = c(cats,"Other_repeat"))
repsByTypeRed$type = factor(repsByTypeRed$type)
repsByTypeRed$wdth = width(repsByTypeRed)
##################################################################################################################################################
repsByTypeRed2 = split(repsByTypeRed, repsByTypeRed$class)
# Make a table by class
Mb = unlist(lapply(repsByTypeRed2, function(x){sum(width(reduce(x)))/1e6}))
Cp = unlist(lapply(repsByTypeRed2, function(x){length(x)}))
GS = strsplit(tbl[4], " ") %>% unlist() %>%  stri_remove_empty
GS = as.numeric(GS[3])/1e6
GP = (Mb/GS)*100
RP = (Mb/sum(Mb))*100
df = data.frame(Mb, Cp, RP, GP, GS)
total = colSums(df)
df = rbind(df, "Total" = total)
df$GS = GS
cat = rownames(df)
df = cbind(cat, df)
rownames(df) = NULL
# Save table
write.table(df, outFile1, col.names = T, row.names = F, quote = F, sep = "\t")
##################################################################################################################################################
# Now I need a table by family
# Let's count using disjoin; here overlappingranges in diffrenet strand will have 
# two type classificactions 
repsByTypeRed3 = disjoin(repsByTypeRed, with.revmap=T, ignore.strand=T)
# Let's assing this bases based on the length of overlapping elements
# Here the longest element will have the overlapping bases
# This is just to count % in the genome
revmap = mcols(repsByTypeRed3)$revmap
# So I need type and width
repsByTypeRed3$type = extractList(mcols(repsByTypeRed)$type, revmap)
repsByTypeRed3$wdth = extractList(mcols(repsByTypeRed)$wdth, revmap)
# Now let's divide elements that do not have this kind of issue
uniques = repsByTypeRed3[unlist(lapply(repsByTypeRed3$revmap, length)) == 1,]
# And the ones with overlaps
dups = repsByTypeRed3[unlist(lapply(repsByTypeRed3$revmap, length)) > 1,]
# Using Map, I'll get the longest element and I'll assing teh respectiv type
dups$type = unlist(Map('[', dups$type, which.max(dups$wdth)))
# Now collapse
redFam = c(uniques, dups)
# Now I can split to count
redFam = split(redFam, redFam$type)
redFamList = split(repsByTypeRed, repsByTypeRed$type)
# Build the table
options(scipen = 999)
Mb = unlist(lapply(redFamList, function(x){sum(width(reduce(x)))/1e6}))
Cp = unlist(lapply(redFamList, function(x){length(x)}))
GS = strsplit(tbl[4], " ") %>% unlist() %>%  stri_remove_empty
GS = as.numeric(GS[3])/1e6
GP = (Mb/GS)*100
RP = (Mb/sum(Mb))*100
df = data.frame(Mb, Cp, RP, GP, GS)
total = colSums(df)
df = rbind(df, "Total" = total)
df$GS = GS
cat = rownames(df)
df = cbind(cat, df)
rownames(df) = NULL
# Save table
write.table(df, outFile2, col.names = T, row.names = F, quote = F, sep = "\t")
##################################################################################################################################################
# Add extra columns
repsByTypeRed$wdth = NULL
repsByTypeRed$ID = paste0(repsByTypeRed$type,":",1:length(repsByTypeRed))

# Order rows and columns
repsByTypeRed_df = as.data.frame(repsByTypeRed)
repsByTypeRed_df = repsByTypeRed_df[,c("seqnames", "start", "end", "width", "strand",
                    "source", "type", "score", "phase", "ID", "rep_name","class")]
repsByTypeRed = makeGRangesFromDataFrame(repsByTypeRed_df, keep.extra.columns = T)
repsByTypeRed = repsByTypeRed[order(seqnames(repsByTypeRed), start(repsByTypeRed)),]

# Save the last object
export(repsByTypeRed, outFile3, "gff3")
##################################################################################################################################################
# Add strand info to type and class
repsByTypeRed$type = as.character(repsByTypeRed$type)
repsByTypeRed$class = as.character(repsByTypeRed$class)
repsByTypeRed$type = ifelse(!repsByTypeRed$type %in% names(unstrandedList), 
                            paste0(repsByTypeRed$type,"_S"), 
                        repsByTypeRed$type)
repsByTypeRed$class = ifelse(!repsByTypeRed$type %in% names(unstrandedList),
                             paste0(repsByTypeRed$class,"_S"), 
                            repsByTypeRed$class)

# Create antisense version
As_repsByTypeRed = repsByTypeRed[grepl("_S", repsByTypeRed$type),]
As_repsByTypeRed$type = sub("_S", "_As", As_repsByTypeRed$type)
As_repsByTypeRed$class = sub("_S", "_As", As_repsByTypeRed$class)
strand(As_repsByTypeRed) = ifelse(strand(As_repsByTypeRed) == "+","-","+")

# Merge and create new IDs
repStrand = c(repsByTypeRed, As_repsByTypeRed)
repStrand = repStrand[order(seqnames(repStrand), start(repStrand)),]
repStrand$ID = paste0(repStrand$type, ":", 1:length(repStrand))
export(repStrand, outFile4, "gff3")
##################################################################################################################################################
