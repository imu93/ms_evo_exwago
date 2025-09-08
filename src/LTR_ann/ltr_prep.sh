#!/bin/bash
# Fisrt I need to define some variables
srcPath=$(pwd)
srcPath=${srcPath}/src
lib=$1 # repeat library
genome=$2 # genome assembly
segments=$3 # segenmented annotation
r_threads=$4 # N therads for R
id=$(echo $genome |sed 's/\..*//')

let RM_THREADS=12/4

if [ -z "$lib" ]; then echo "ERROR: No library provided"; exit 1; fi
if [ -z "$genome" ]; then echo "ERROR: No genome provided"; exit 1; fi
if [ -z "$segments" ]; then echo "ERROR: No segmentation file provided"; exit 1; fi
if [ -z "$r_threads" ]; then echo "ERROR: No thread count provided"; exit 1; fi

# Subset only LTRs
Rscript $srcPath/01.get_ltr_from_lib.R $lib
ltr_lib=$(echo *ltr_cons.fa)

# Check point
if [ -e "$ltr_lib" ]; then
    echo "Process $ltr_lib"
else
    echo "ERROR: there is not ltr lib in this directory."
    sleep 5
    exit 1
fi

echo "#########################  Split repeat Library #################################"
mkdir lib
Rscript $srcPath/01.split_lib.R $ltr_lib
mv *LTR.fa lib/

# Now run TE-Aid and produce self dot plots
teaid_path=$(which TE-Aid)
if [ -z "$lib" ]; then
    echo "TE-Aid is not availabe"
    sleep 5
    exit 1
fi

mkdir self_blast
mkdir pdf
mkdir index
echo "#########################  TE-aid for selfblast tables  #################################"

for file in lib/*.fa; do
TE-Aid -q $file -g $genome -m 200 -t -o ./
done
mv *.pdf pdf/
mv *pairs.txt self_blast/
rm *blastp.out
rm *blastn.out
rm *orfs.fasta
rm *orftable.txt
mv *.fa.n* index/


echo "#########################  Build LTR library with structural information  #################################"
#Preppare library and inputs for RepeatMasker
Rscript $srcPath/01.split_LTR_models.R $ltr_lib
Rscript $srcPath/02.get_LTR_fasta.R $genome $segments

echo "#########################  Reannotate segements unsig RepeatMasker  #################################"
# Run RepeatMasker using the LTR segements
LTR_LIBRARY=$(echo *ltr_library.fa)
LTR_ANNOTATION=$(echo *genomic.LTR_segements_final.fa)
RM_DIR=${id}_repeatmasker
mkdir $RM_DIR
RepeatMasker -pa $RM_THREADS -nolow -cutoff 400 -lib $LTR_LIBRARY $LTR_ANNOTATION
# I will also implement a cutoff just in case
mv *.tbl $RM_DIR/
mv *.masked $RM_DIR/
mv *cat.gz $RM_DIR/
echo "#########################  Run TEsorter to look for full-length elements  #################################"
# And TEsorter
bash $srcPath/03.TEsorter.sh $r_threads
ln -s ${id}_TEsorter_results/*rexdb-metazoa.cls.tsv $(pwd)

echo "#########################  GFF3 from RM outFile  #################################"

LTR_ANN=$(echo *LTR_segements_final.fa.out)

# Format into GFF3
Rscript $srcPath/rm2gff3.R -i ${LTR_ANN}
mv ./*.out  $RM_DIR/

RM_ANN=$(echo *LTR_segements_final.fa.gff3)

echo "#########################  Transfer regions  #################################"
# Transfer regions
Rscript $srcPath/03.transfer_relative_annotation.R -l ${RM_ANN} -s $segments  -o ${id}.ltr_transf.gff3 -t $r_threads
Rscript $srcPath/03.clean_transfered.R ${id}.ltr_transf.gff3 $segments
# Note: this filter based on divergence and fragmentation.
# RM can produce divergent hits of similar families, thus I'll clean.

echo "#########################  LTR structural clssification  #################################"
# Now look for soloLTRs and fl enelements
TEsorter_tab=$(echo *rexdb-metazoa.cls.tsv)
echo "#########################  FL  LTR classification  #################################"
Rscript $srcPath/04.flLTRFinder.R -l ${id}.ltr_transf.gff3 -s $segments -o ${id}.fl_LTR_struc.gff3 -t $TEsorter_tab

echo "#########################  soloLTR annotation  #################################"
Rscript $srcPath/04.soloLTRFinder_parallel.R -l ${id}.ltr_transf.gff3 -d $segments -o ${id}.soloLTR.gff3 -t $r_threads

Rscript $srcPath/04.pre_mergeSegments.R ${id}.soloLTR.gff3 ${id}.fl_LTR_struc.gff3
deTab=$(echo *rRNA_deTbl.txt)
Rscript $srcPath/04.merge_Segements.R ${id}.soloLTR.gff3 ${id}.fl_LTR_struc.gff3 ${id}.ltr_transf.gff3 $deTab
