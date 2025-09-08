pacman::p_load(optparse, rtracklayer, plyranges, dplyr, doParallel, foreach, stringr)
list.files()

option_list = list(
  make_option(c("-l", "--ltr_gff3"), type = "character", default = NULL,
              help = "RepeatMasker annotation with splited LTRs"),
  make_option(c("-d", "--merged_gff"), type = "character", default = NULL,
              help = "Full repreat annotation after assign overlaps (merged_repeats_unstr)"),
  make_option(c("-o", "--out_gff"), type = "character", default = NULL,
              help = "Out file name of soloLTR annotation"),
  make_option(c("-t", "--threads"), type = "numeric", default = 1,
              help = "Number of threads")
)

# Parse command line arguments
opt = parse_args(OptionParser(option_list = option_list))

# Check if input file is provided
# Validate input arguments
if (is.null(opt$ltr_gff3)) {
  stop("Error: No LTR GFF3 file provided. Use --ltr_gff3 to specify the input file.")
}

if (is.null(opt$merged_gff)) {
  stop("Error: No merged repeat annotation provided. Use --merged_gff to specify the input file.")
}

if (is.null(opt$out_gff)) {
  stop("Error: No out file name provided. Use --out_gff to specify the input file.")
}

if (is.null(opt$threads)) {
  stop("Error: No num threads specified")
}

# This script aims to identify soloLTRs
# After using self-blast I was able to identify the LTR region of most of my models
# Thus, I can distinguish between the LTR region and the internal region.  
# Also I have used TEsorte to look for TEs with coding potential

# To identify soloLTRs I should first remove false positives form my LTR annotation
# For this I will use my RepeatCraft annotation

# Let's start reading the input annotations
# RepeatMasker annotation using LTR-CDS-LTR nomenclature
ltr_gff = import(opt$ltr_gff3)
# Merged annotation
rep_ann = import(opt$merged_gff)
# Define the out file
outFile =  opt$out_gff
# N threads
THREADS = opt$threads

numCores = THREADS  # Leave one core free
cl = makeCluster(numCores)
registerDoParallel(cl)

# Only thoes non-LTR repeats
oth_ann = rep_ann[!grepl("LTR", rep_ann$type),]
# Now remove if there is an overlap
# First get false positives (just in case)
fp_ltr = ltr_gff[findOverlaps(ltr_gff,oth_ann, type = "within") %>% queryHits() %>% unique()]
# And remove 
ltr_gff = ltr_gff[!1:length(ltr_gff) %in% c(findOverlaps(ltr_gff,oth_ann, type = "within", minoverlap = 50) %>% queryHits())]

# also remove hits with DIRS
cand_solo = ltr_gff[ltr_gff$type == "LTR",]
cand_solo = cand_solo[!grepl("DIRS", cand_solo$class),]
#cand_solo = sample(cand_solo, 1000)

# Get cadidates 
csoloList = split(cand_solo, cand_solo$ID)


distance_threshold = 1000
filter_candidates = foreach(i = names(csoloList), .packages = c("GenomicRanges", "plyranges")) %dopar% {
  cat("Processing:", i, "\n")
  tmp.x = csoloList[[i]]  # Extract candidate
  tmp.comp = ltr_gff[ltr_gff$ID != tmp.x$ID,]  # Remove candidate from the whole list
  
  # Find neighbours
  dst_down = join_follow_upstream(tmp.x, tmp.comp)
  dst_up = join_precede_downstream(tmp.x, tmp.comp)
  down_neigh = tmp.comp[tmp.comp$ID == dst_down$ID.y,]
  up_neigh = tmp.comp[tmp.comp$ID == dst_up$ID.y,]
  tmp.neigbours = c(down_neigh, up_neigh)
  
  # If there are no neighbors, keep tmp.x by default
  if (length(tmp.neigbours) == 0) return(tmp.x)
  
  # If the LTR region is within 2 CDS of different category
  if (sum(tmp.neigbours$LTR_type == "int") == 2 & 
      sum(tmp.neigbours$rep_name != tmp.x$rep_name) == 2) return(tmp.x)
  
  # Estimate the distance
  tmp.neigbours$distance = GenomicRanges::distance(tmp.x, tmp.neigbours)
  
  # If neighbors are too far (>10kb), keep tmp.x
  if (sum(tmp.neigbours$distance > 10e3) == 2) return(tmp.x)
  
  # If neigbour is from the same family and close enough and is int discrt
  if (any(tmp.neigbours$LTR_type == "int" & 
          tmp.neigbours$rep_name == tmp.x$rep_name & 
          tmp.neigbours$distance <= 200, na.rm = TRUE)) {
    return(NULL)  # Discard candidate
  }
  
  # Same for putative trims
  if (any(tmp.neigbours$LTR_type == "LTR" & 
          tmp.neigbours$rep_name == tmp.x$rep_name & 
          tmp.neigbours$distance <= 200, na.rm = TRUE)) {
    return(NULL)  # Discard candidate
  }
  
  # Extract metadata as data frame
  neigbour_metadata = mcols(tmp.neigbours)
  
  # Check filtering conditions
  distance_condition = neigbour_metadata$distance < distance_threshold
  rep_name_condition = neigbour_metadata$rep_name == tmp.x$rep_name
  
  # If at least one condition is TRUE for any neighbor, discard tmp.x
  if (any(distance_condition | rep_name_condition, na.rm = TRUE)) {
    return(NULL)
  } else {
    return(tmp.x)
  }
}
stopCluster(cl)  # Stop the cluster

# Remove NULL values
filter_candidates = filter_candidates[!sapply(filter_candidates, is.null)]
names(filter_candidates) = NULL
filter_candidates = do.call(c, filter_candidates)
filter_candidates = filter_candidates[order(seqnames(filter_candidates), start(filter_candidates)),]
filter_candidates$color = "#078F1A"
filter_candidates$type = "soloLTR"
filter_candidates$LTR_type = "soloLTR"
filter_candidates$ID = gsub("LTR", "soloLTR", filter_candidates$ID)
target_str = filter_candidates$segment_ID %>% strsplit("_") %>% 
  lapply(function(x){x[length(x)]}) %>% unlist() %>% 
  gsub(":.*", "", .)
#filter_candidates$ID = str_replace(filter_candidates$ID,":", paste0("_", target_str, ":"))
filter_candidates$igv = paste0(as.character(seqnames(filter_candidates)), ":", 
                               start(filter_candidates), "-",end(filter_candidates))


width(filter_candidates) %>% hist(main="soloLTR length")

export(filter_candidates, outFile)
