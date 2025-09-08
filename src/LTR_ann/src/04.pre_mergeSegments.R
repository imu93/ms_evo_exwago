# This script is post solo and fl_perdictions
# It aims to remove possible overaps between soloLTRs and full length LTRs
pacman::p_load(rtracklayer)
# Import files
args = commandArgs(trailingOnly = T)
soloFile = args[1]
flLTRFile = args[2]
soloLTR = import(soloFile)
flLTR = import(flLTRFile)
# Remove overlaping soloLTRs
soloLTR = soloLTR[!soloLTR %over% flLTR,]
# Over write the original soloLTR annotation
export(soloLTR,soloFile, "gff3")
