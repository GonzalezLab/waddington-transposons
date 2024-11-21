#!/bin/bash

module load cesga/system miniconda3/22.11.1-1
conda activate waddington

# Run

DIR="/home/csic/gcy/mcz/scratch/Waddington/FIX"

cd $DIR

types="D907_A D907_N Mix_A Mix_N Mix2_A Mix2_N"

awk '$3 == "gene" '  ${DIR}/referenceGenome/dmel-all-r6.31.gtf | grep -E "^2L|^2R|^3L|^3R|^X" | sort -k1,1 -k4,4n -k5,5n > ${DIR}/referenceGenome/dmel-all-r6.31.gene.gtf

mkdir ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/

for type in $types; do
> ${DIR}/ANALYSES/post_analysis_R01/denovoTEresults/denovoTE.infoIntersect.$type.tab
type2="$type"
if [[ "$type" == "Mix2_A" ]]; then
type2="Mix_A"
elif [[ "$type" == "Mix2_N" ]]; then
type2="Mix_N"
fi

i=1

# Read each line in the input file
while read insertion; do

# Save the current insertion to a temporary BED file
echo "$insertion" > ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.bed 

# Run bedtools window with a 1kb window
bedtools window -a ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.bed \
-b ${DIR}/referenceGenome/dmel-all-r6.31.gene.gtf \
-w 1000 > ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.intersect 
#cat ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.intersect 

# Check if the bedtools window output is empty
if [ ! -s ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.intersect  ]; then
# If empty, run bedtools closest instead
bedtools closest -a ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.bed \
-b ${DIR}/referenceGenome/dmel-all-r6.31.gene.gtf -d > ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.window 
gene=$(cat ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.window | cut -f 16 | cut -f 2 -d' ' | tr -d '"' | tr -d ';')
distance=$(cat ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.window | cut -f 17 ) 

if grep -q "$gene" $DIR/list_genes_DE/genes_*_${type2}.lst; then
DEG=1
else
DEG=0
fi

result=$(echo "${gene}:${distance}bp ($DEG)")
else
genes=$(cat ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.intersect  | cut -f 16 | cut -f 2 -d' ' | tr -d '"' | tr -d ';')
result=""

while read gene
do
grep -w $gene ${DIR}/referenceGenome/dmel-all-r6.31.gene.gtf > ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.gene.$gene.gff
distance=$(bedtools closest -a ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.bed -b ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.gene.$gene.gff -d | cut -f 17 )

DEG=""

# Loop through all matching files
for file in $DIR/list_genes_DE/genes_*_${type2}.lst; do
    if grep -q "$gene" "$file"; then
        # Extract the condition from the filename and append to DEG
        case=$(basename "$file" | sed -E "s/genes_([A-Za-z]+_[A-Za-z]+)_${type2}\.lst/\1/")
        if [[ -n "$DEG" ]]; then
            DEG="$DEG, $case"
        else
            DEG="$case"
        fi
    fi
done

# If DEG is still empty, set it to "0"
if [[ -z "$DEG" ]]; then
    DEG=0
fi

result=$(echo "${gene}:${distance}bp ($DEG);$result")
done <<< "$genes"
fi

H3K27ac=""

bedtools closest \
-a ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp/$type.denovoTE.freq_0.1.$i.bed \
-b ${DIR}/ANALYSES/post_analysis_R01/H3K27ac_allMixLinesMerged.fix.bed -d \
> ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp//$type.denovoTE.freq_0.1.$i.closest_H3K27ac.bed
H3K27ac=$(cat ${DIR}/ANALYSES/post_analysis_R01/denovoTE/tmp//$type.denovoTE.freq_0.1.$i.closest_H3K27ac.bed | cut -f8-11  | awk -F'\t' '{print $1 ":" $2 "-" $3 "=" $4}')



# Chr Start End TEName PopName Freq Gene Distance
infoIns=$(echo "$insertion" | cut -f 1-5,7)
echo -e "$infoIns\t$result\t$H3K27ac" >> ${DIR}/ANALYSES/post_analysis_R01/denovoTEresults/denovoTE.infoIntersect.$type.tab

# Increment the counter
i=$(( i + 1 ))

done < ${DIR}/ANALYSES/post_analysis_R01/denovoTE/$type.denovoTE.freq_0.1.bed

done 






# cruzar con genes 1kb
# mas lejos si no hay (closest)
# cruzar DEG
# cruzar TE con enhancer
# mirar que gen es

