## Scripts related to structural-based LTR annotation
1.ltr_prep.sh master scrip for LTR annotation
  Requirements: curated repeat library, genome, segmented annotation.

2. **src/** contains all the required scripts to classify LTR segements
    ├── 01.get_ltr_from_lib.R : retrieve only LTR families from curated library
    ├── 01.split_lib.R : get individual models
    ├── 01.split_LTR_models.R : split LTR retrotranspsosns based on self-blast
    ├── 02.get_LTR_fasta.R : get sequences
    ├── 03.clean_transfered.R : process RepeatMasker results
    ├── 03.TEsorter.sh : Run TEsorter for LTR protein annotation
    ├── 03.transfer_relative_annotation.R : Transfers LTR segements to structures
    ├── 04.flLTRFinder.R : identify fl LTRs based on TEsorter
    ├── 04.merge_Segements.R : merge new LTRs with the rest of the genome annotation
    ├── 04.pre_mergeSegments.R : prepares new annotation to be merged, just some checkpoints 
    ├── 04.soloLTRFinder_parallel.R :  looks for soloLTRs based on their surrounding LTRs, considering family, structure and distance
    └── rm2gff3.R : parse RM.out and builds gff3
 
