# This script aims to identify full length LTR rtreotransposons and
# extract their LTR and in ternal region

# I will use the deffinition of complete  based on TEsorter
pacman::p_load(rtracklayer, plyranges, dplyr, optparse)
list.files()

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
  make_option(c("-t", "--tesorter"),
              type = "character",
              help = "Number of threads")
)

opt = parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$ltr_rm), !is.null(opt$seg_gff), !is.null(opt$out_gff))

inLTR  = opt$ltr_rm
inSegments = opt$seg_gff
tsTab= opt$tesorter

#inLTR = "../../../annotations/nippostrongylus_brasiliensis/nippostrongylus_brasiliensis.ltr_transf.gff3"
#inSegments = "../../../annotations/nippostrongylus_brasiliensis/nippostrongylus_brasiliensis.PRJNA994163.genomic.segment_final.gff3.gz"
#tsTab = "../../../annotations/nippostrongylus_brasiliensis/nippostrongylus_brasiliensis.genomic.LTR_segements_final.fa.rexdb-metazoa.cls.tsv"

# Read files
# Fiest the segemented annotation
inSegments = import(inSegments)
# RepeatMasker annotation with the splited LTR library
inRmGff = import(inLTR)

# TEsorter results using rexdb
TEsorter = read.delim(tsTab)

# Extract complete LTRs
complete_LTRs = TEsorter[TEsorter$Complete == "yes",]
fl_ltrs = inSegments[match(complete_LTRs$X.TE, inSegments$ID),]
# Subset by ovelraps only thoes with within
fl_over =  subsetByOverlaps(inRmGff, fl_ltrs, type = "within")
# Lets get the LTR region
fl_over_LTR = fl_over[fl_over$LTR_type == "LTR",]

# And the internal 
fl_over_int = fl_over[fl_over$LTR_type == "int",]

# Now lets eddit type
fl_over_LTR$type = paste0(fl_over_LTR$type, "_fl")
fl_over_int$type = paste0(fl_over_int$type, "_fl")
# Same for LTR type
fl_over_LTR$LTR_type = fl_over_LTR$type
fl_over_int$LTR_type = fl_over_int$type

# Now just merge
fl =  c(fl_over_LTR, fl_over_int)
export(fl, opt$out_gff, "gff3")
