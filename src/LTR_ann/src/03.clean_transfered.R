# This script aims to clean transfered regions of spurious RepeatMasker predictions
pacman::p_load(rtracklayer, plyranges)
args = commandArgs(trailingOnly = T)
trans_ann = import(args[1]) # transfered annotation
rep_ann = import(args[2]) # segmented annotation
# My logic is that I can have spurious hits of diffrent families compared to the
# original one in the DEA becuse how I collapsed the with repeatcraft and also 
# due to random hits given sequence similarity

# I expect this hits to be:
# 1. smaller than the original region
# 2. And divergent

# This I will use these critearia to filer
# First add length
trans_ann$wdt = width(trans_ann)
rep_ann$wdt = width(rep_ann)
rep_ann = rep_ann[grepl("LTR", rep_ann$type),]
# If the hit is less than a quarter of the original region is a possible candidate
# But this does not deffine that this is a bad hit
rep_ann$tird_part = rep_ann$wdt/4
new_ann = join_overlap_intersect_within_directed(trans_ann, rep_ann)
# So if it is also from an other familiy this is a true red flag
inter = new_ann[new_ann$rep_name.x != new_ann$rep_name.y]
# Now filter by divergence
fp = inter[inter$wdt.x <= inter$tird_part & inter$prec_div >= 15,]
# And subset 
clean_transf = trans_ann[match(setdiff(trans_ann$ID, fp$ID.x), trans_ann$ID),]
# Export
export(clean_transf, args[1], "gff3")