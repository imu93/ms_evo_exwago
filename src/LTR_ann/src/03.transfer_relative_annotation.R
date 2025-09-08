# This script aims to split LTR segments based on their structure
# considering the LTR region and the Internal/Coding region
pacman::p_load(optparse, rtracklayer, plyranges, dplyr, stringr, future.apply)
# Parse arguments
option_list = list(
  make_option(c("-l", "--ltr_rm"),
              type = "character",
              help = "RepeatMasker annotation with split LTRs"),
  make_option(c("-s", "--seg_gff"), 
              type = "character",
              help = "Segmented repeat annotation"),
  make_option(c("-o", "--out_gff"),
              type = "character",
              help = "Output filename for structured LTRs"),
  make_option(c("-t", "--threads"),
              type = "numeric",
              help = "Number of threads")
)
opt = parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$ltr_rm), !is.null(opt$seg_gff), !is.null(opt$out_gff))

# Import GFF3 files
ltr_gff = import(opt$ltr_rm)
rep_ann = import(opt$seg_gff)

#ltr_gff = "../../../annotations/nippostrongylus_brasiliensis/nippostrongylus_brasiliensis.genomic.LTR_segements_final.fa.gff3"
#rep_ann = "../../../annotations/nippostrongylus_brasiliensis/nippostrongylus_brasiliensis.PRJNA994163.genomic.segment_final.gff3.gz"
#ltr_gff = import(ltr_gff)
#rep_ann = import(rep_ann)

# Fix encoded IDs on the repeatmasker files
ltr_gff$ID = as.character(seqnames(ltr_gff)) %>% str_replace("%2f", "/")

# Keep only LTR elements from the segemented annotation
ltr_ann = rep_ann[grepl("LTR", rep_ann$type)]

# Parallel setup
future::plan(multisession, workers = opt$threads)

# Process each LTR element in parallel
ids = unique(ltr_ann$ID)
n_ltr_lst = future_lapply(ids, function(sq) {
  # For each LTR I need:
  # The original segement
  #sq = "LTR/Pao_rnd-1_family-371_As:191218"
  tmp.repann = ltr_ann[ltr_ann$ID == sq]
  # All their possible new annotations
  tmp.gffann = ltr_gff[ltr_gff$ID == sq]
  # Now If no new annotations, save the original
  # One kind of these TEs could be very small fragemnts after assigning overlapping bases
  if (length(tmp.gffann) < 1) return(tmp.repann)
  # Assign coordinates based on strand-aware logic
  
  # Firts I need the info from the original element
  # Strand
  feat_strand = as.character(strand(tmp.repann))
  # Start
  rep_start = start(tmp.repann)
  # End
  rep_end = end(tmp.repann)
  # And the same for the new annotations
  feat_start = start(tmp.gffann)
  rep_width = width(tmp.gffann)
  # Now this is the core of the script
  if (feat_strand == "+") {
    # If strand + the start is the original start of the inital copy plus 1 
    # which will be the srat of the RM annotation and -1 to fix this
    # And the end will be new start plus the width of the new annotation
    n_start = rep_start + feat_start - 1
    n_end = n_start + rep_width - 1
  } else if (feat_strand == "-") {
    # Now a bit more confusing beacuse the ranges has the same numbers in + and - 
    # but this is the solution if strand == -
    # The new end will be the original end - the new start
    # And again just the width for the start
    n_end = rep_end - feat_start + 1
    n_start = n_end - rep_width + 1
  } else {
    stop("Unrecognized strand: ", feat_strand)
  }
  # Consistent seqname and strand
  n_sqn = rep(as.character(seqnames(tmp.repann)), length(tmp.gffann))
  n_strand = rep(feat_strand, length(tmp.gffann))
  # Same strand as the original segement
  n_strand =  rep(as.character(strand(tmp.repann)), length(tmp.gffann))
  # Build the new GRanges
  n_ltr_gff = GRanges(seqnames = n_sqn, IRanges(start = n_start, end = n_end), strand = n_strand)
  # Add meta-columns
  n_ltr_gff$source = "RepeatMasker"
  n_ltr_gff$type = "Undetermined"
  n_ltr_gff$type = ifelse(grepl("CDS", tmp.gffann$rep_name), "int", n_ltr_gff$type)
  n_ltr_gff$type = ifelse(grepl("LTR", tmp.gffann$rep_name), "LTR", n_ltr_gff$type)
  n_ltr_gff$score = NA
  n_ltr_gff$phase = NA
  n_ltr_gff$SW_score = tmp.gffann$SW_score
  n_ltr_gff$prec_div = tmp.gffann$perc_div
  n_ltr_gff$perc_del = tmp.gffann$per_del
  n_ltr_gff$perc_ins = tmp.gffann$perc_ins
  n_ltr_gff$rep_name = tmp.gffann$rep_name
  # Fix the rep name
  n_ltr_gff$rep_name = sub("_LTR|_CDS", "", n_ltr_gff$rep_name)
  n_ltr_gff$ID = tmp.gffann$ID
  
  # Add the class
  n_ltr_gff$class =  sub(":.*", "", n_ltr_gff$ID) %>% sub("_rnd.*", "", .)  
  n_ltr_gff$str = strsplit(n_ltr_gff$ID, "_") %>% lapply(function(x){x[length(x)]}) %>% unlist() %>% sub(":.*", "", .)
  return(n_ltr_gff)
})

# Merge and finalize
fix_ltr_gff = do.call(c, n_ltr_lst)
fix_ltr_gff$LTR_type = fix_ltr_gff$type
fix_ltr_gff$color = case_when(
  fix_ltr_gff$type == "LTR" ~ "#009C17",
  fix_ltr_gff$type == "int" ~ "#E59C17",
  TRUE ~ "#545454"
)
fix_ltr_gff$segment_ID = fix_ltr_gff$ID
fix_ltr_gff$ID = paste0(fix_ltr_gff$type, "_", fix_ltr_gff$str, ":", seq_along(fix_ltr_gff))

# Export to GFF3
export(fix_ltr_gff, opt$out_gff, format = "gff3")