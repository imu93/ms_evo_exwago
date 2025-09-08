#!/bin/bash
#$ -V
#$ -N LTR_finder
#$ -cwd
#$ -j yes
#$ -o mr.$JOB_ID.log
#$ -pe smp64 20


SCRIPT_DIR=/home/isaac/Software/merge_repeats
threads=10
genome=$(echo *.genomic.fa)
rmOut=$(echo *.fa.out)
table=$(echo *.tbl)
gffFile=$(echo *.gff)
spe=$(echo $genome | sed 's/\..*//')
genSize=$(sed -n '4p' $table | rev | cut -f1,1 -d ':' | rev | sed 's/ bp.*//g; s/ //g')

echo "####Create_dirs####"
mkdir LtrFinder/
cd LtrFinder/
ln -s ../${genome} ./
ln -s ../${rmOut} ./
ln -s ../${table} ./
ln -s ../${gffFile} ./

echo "####Ed_genome_ID####"
# Edit genome ID
Rscript ${SCRIPT_DIR}/ed_genome_id.R $genome
genometmp=$(echo *.tmp.fa)

#LTR_finder_parallel
echo "####LTR_finder####"
${SCRIPT_DIR}/LTR_FINDER_parallel/LTR_FINDER_parallel -seq $genometmp -threads $threads
rm -rf $genometmp 
cat ${gffFile} | sed '1,2d; s/Target=/Target "Motif:/g' | awk '{OFS="\t"}{gsub(/$/,"\"", $10); print}' |sed 's/""/"/'|sed 's/Target\t/Target /g; 1s/^/##gff-version 2\n##date xxxx-xx-xx\n/; s/\(.*\)\t/\1 /; s/\(.*\)\t/\1 /' > ${spe}.out.gff2
sed -i '/##sequence-region.*/d' ${spe}.out.gff2
mkdir ${spe}LtrFinder/
cp *fa.finder.combine.gff3 ${spe}LtrFinder/
cd ${spe}LtrFinder/
gffDir=$(pwd)
inGFF3=$(echo *.finder.combine.gff3)
cd ..
Rscript ${SCRIPT_DIR}/create_cfg.R ${gffDir}/ $inGFF3


#Repeatcraft
echo "####rc####"
python ${SCRIPT_DIR}/repeatCraft/repeatcraft.py -r ${spe}.out.gff2 -u $rmOut -c ${SCRIPT_DIR}/repeatCraft/example/repeatcraft.cfg -o ${spe} -m loose
#sort -s -k1,1 -k7,7 ${spe}.rmerge.gff > ${spe}.rmerge.gff.sorted
#filter-gff overlap -v --progress -d -s 1 -t -c length -a length ${spe}.rmerge.gff.sorted ${spe}.rmerge.gff.filtered

