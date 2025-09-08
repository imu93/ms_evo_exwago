pacman::p_load(Biostrings)
args=commandArgs(trailingOnly = TRUE)
inFile=args[1]
fa = readDNAStringSet(inFile)
names(fa) =  sub(" .*", "", names(fa))

for (i in names(fa)) {
  sq_n = sub("#", "_", i)
  outFile = paste(sub("\\/.*$", "", sq_n), ".fa", sep = "")
  outFile = sub("\\?", "", outFile)
  outFile = sub(" ", "", outFile)
  f = DNAStringSet(fa[[i]])
  names(f) = sub("\\#", "__", i)
  writeXStringSet(f, outFile)
}
