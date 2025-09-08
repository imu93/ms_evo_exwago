if (!require("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(dplyr, optparse, rtracklayer, stringr)
# Define colors
cols <- c(
  "LINE" = "#78DB3F",
  "SINE" = "#FFB200",
  "DNA|RC" = "#00C3DB",
  "LTR" = "#00AB14",
  "Low_complexity" = "#d1d1e0",
  "Satellite" = "#ff99ff",
  "Simple_repeat" = "#8686ac",
  "Unknown" = "#D2F102"
)

# Function to parse RepeatMasker output to GFF3 format with colors
parse_repeatmasker_to_gff3 = function(input_file) {
  data = readLines(input_file) 
  data = data[4:length(data)]
  asterisk = ifelse(grepl("\\*$", data), "*", NA) # asterisk indicate overlaps between elements 
  data = str_replace(data,"\\*$", "")
  data = strsplit(data, "\\s+") %>% lapply(function(x){x[x != ""]})
  data =  do.call(rbind, data) %>% as.data.frame() 
  data$overlap = asterisk
  colnames(data) = c("SW_score", "perc_div", "per_del", "perc_ins", "seqnames",
                     "start", "end", "ext1", "strand", "rep_name","class", 
                      "subject_start", "subject_end", "ext2", "ID", "Overlap")
  
  data$strand =  ifelse(data$strand == "C", "-", data$strand)
  gff3 = makeGRangesFromDataFrame(data, keep.extra.columns = T)
  gff3$color = "Other"
  for (i in names(cols)) {
    tmp.col = cols[[i]] 
    gff3$color = ifelse(grepl(i, gff3$class), tmp.col, gff3$color)
  }

  return(gff3)
}

# Command line options
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Input RepeatMasker .out file")
)

# Parse command line arguments
opt = parse_args(OptionParser(option_list = option_list))

# Check if input file is provided
if (is.null(opt$input)) {
  stop("Error: No input file provided. Use --input to specify the input file.")
}

# Parse RepeatMasker output to GFF3 with colors
gff3 = parse_repeatmasker_to_gff3(opt$input)

# Write GFF3 file
output_file = sub(".out", ".gff3", opt$input)
export(gff3, output_file, format = "gff3")

cat("GFF3 file created:", output_file, "\n")
