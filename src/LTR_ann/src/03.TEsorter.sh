#!/bin/bash

thr=$1
for file in *segements_final.fa;do
ID=$(echo $file | sed "s/\..*//")
TEsorter $file -p $thr -cov 60 -eval 1e-5 --hmm-database rexdb-metazoa
outDir=${ID}_TEsorter_results
mkdir -p $outDir
chmod 775 $outDir
mv *rex* $outDir

cd $outDir
outRT=${ID}.rexdb.dom.RT.faa
outENDO=${ID}.rexdb.dom.ENDO.faa
outPROT=${ID}.rexdb.dom.PROT.faa
outGAG=${ID}.rexdb.dom.GAG.faa
outINT=${ID}.rexdb.dom.INT.faa
outRH=${ID}.rexdb.dom.RH.faa
outHEL2=${ID}.rexdb.dom.HEL2.faa


cat *metazoa.dom.tsv | grep -P "\-RT\t" | get_record.py -i *metazoa.dom.faa -o $outRT -t fasta
cat *metazoa.dom.tsv | grep -P "\-ENDO\t" | get_record.py -i *metazoa.dom.faa -o $outENDO -t fasta
cat *metazoa.dom.tsv | grep -P "\-PROT\t" | get_record.py -i *metazoa.dom.faa -o $outPROT -t fasta
cat *metazoa.dom.tsv | grep -P "\-GAG\t" | get_record.py -i *metazoa.dom.faa -o $outGAG -t fasta
cat *metazoa.dom.tsv | grep -P "\-INT\t" | get_record.py -i *metazoa.dom.faa -o $outINT -t fasta
cat *metazoa.dom.tsv | grep -P "\-RH\t" | get_record.py -i *metazoa.dom.faa -o $outRH -t fasta
cat *metazoa.dom.tsv | grep -P "\-HEL2\t" | get_record.py -i *metazoa.dom.faa -o $outHEL2 -t fasta

cat $outRT $outENDO $outPROT $outGAG $outINT $outRH $outHEL2 $outTPase > ${ID}_repeat.dom.faa

cd ..
done
