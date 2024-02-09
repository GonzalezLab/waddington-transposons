#!/bin/bash

#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --job-name=process_frq
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH -o logs/denovoTE%a.out
#SBATCH -e logs/denovoTE%a.err  
#SBATCH --exclude=mr-03-[01-26],mr-03-28,mr-04-00

module load Miniconda3/4.7.10

source activate anopheles

module load R/3.6.0-foss-2018b

# Run

# sbatch --array=1-3 process_freq_all.sh

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

# de novo TE
# they can't be in the other lines nor C nor P0 nor other-selected (but in the selected they can)

types='D907_A:D907_P0|D907_N|D907_C|Mix_C|Mix2_C|Mix_P0|Mix2_P0|Mix_N|Mix2_N D907_N:D907_P0|D907_A|D907_C|Mix_C|Mix2_C|Mix_P0|Mix2_P0|Mix_A|Mix2_A Mix_A:Mix_P0|Mix_N|Mix_C|D907_C|Mix2_C|D907_P0|Mix2_P0|D907_N|Mix2_N Mix_N:Mix_P0|Mix_A|Mix_C|D907_C|Mix2_C|D907_P0|Mix2_P0|D907_A|Mix2_A Mix2_A:Mix2_P0|Mix2_N|Mix2_C|D907_C|Mix_C|D907_P0|Mix_P0|D907_N|Mix_N Mix2_N:Mix2_P0|Mix2_A|Mix2_C|D907_C|Mix_C|D907_P0|Mix_P0|D907_A|Mix_A'

#types='D907_A:D907_P0|D907_N|D907_C:Mix D907_N:D907_P0|D907_A|D907_C:Mix Mix_A:Mix_P0|Mix_N|Mix_C:D907|Mix2 Mix_N:Mix_P0|Mix_A|Mix_C:D907|Mix2 Mix2_A:Mix2_P0|Mix2_N|Mix2_C:D907|Mix Mix2_N:Mix2_P0|Mix2_A|Mix2_C:D907|Mix'



for type in $types; do
class=$(echo $type | cut -f1 -d':')
exclude=$(echo $type | cut -f2 -d':')
#typeExclude=$(echo $type | cut -f3 -d':')
echo $class $exclude

grep "$class" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "$exclude" >  ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.unique.bed

> ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.unique.freq.bed
> ${DIR}/ANALYSES/post_analysis/$class.denovo.TEMP2.tab
> ${DIR}/ANALYSES/post_analysis/$class.denovo.PopTE2.tab

while read insertions
do
chr=$(echo "$insertions" | cut -f1)
start=$(echo "$insertions" | cut -f2)
end=$(echo "$insertions" | cut -f3)
TE=$(echo "$insertions" | cut -f4)
awk -v TE="$TE" -v chr="$chr" -v start="$start" -v end="$end" ' $3 == chr && $4 == start && $5 == end && $2 == TE ' ${DIR}/ANALYSES/post_analysis/frequency_insertions_all_${class}.tab | awk -v OFS='\t' '{print $3, $4, $5, $2, $1, $6, $7}' >> ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.unique.freq.bed
done < ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.unique.bed

sed -i "s/TEMP2\t$/TEMP2\tND/g" ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.unique.freq.bed
sed -i "s/PopTE2\t$/PopTE2\tND/g" ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.unique.freq.bed


excludeParse=$(echo "$exclude" | tr '|' ' ')
> ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.exclude.tmp.bed
for exclude in $excludeParse
do
#echo $exclude
cat $DIR/ANALYSES/common_insertions/$exclude/$exclude.teinsertion.eu.popte2.temp2.bed >> ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.exclude.tmp.bed
done
cat ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.exclude.tmp.bed | sort -k1,1 -k2,2n > ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.exclude.bed 

bedtools window -v -w 20 -a <(cat ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.unique.freq.bed ) -b ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.exclude.bed  > ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed
#bedtools window -v -w 20 -a ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.unique.freq.bed -b ${DIR}/ANALYSES/post_analysis/tmp/denovo/$class.exclude.bed | grep -v "$typeExclude" > ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed

sed -i "s/${class}=//g" ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed
nTEMP2=$(grep "TEMP2" ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed | wc -l )
nPopTE2=$(grep "PopTE2" ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed | wc -l )
nL_TEMP2=$(awk ' $7 > 0.1 ' ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed | grep "TEMP2" |grep -vw "ND" | wc -l)
nL_PopTE2=$(awk ' $7 > 0.1 ' ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed | grep "PopTE2" | grep -vw "ND" | wc -l)

nM_TEMP2=$(awk ' $7 > 0.3 ' ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed | grep "TEMP2" |grep -vw "ND" | wc -l)
nM_PopTE2=$(awk ' $7 > 0.3 ' ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed | grep "PopTE2" |grep -vw "ND" | wc -l)

nH_TEMP2=$(awk ' $7 > 0.5 ' ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed | grep "TEMP2" |grep -vw "ND" | wc -l)
nH_PopTE2=$(awk ' $7 > 0.5 ' ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed |grep "PopTE2" |grep -vw "ND" |  wc -l)
echo -e "$class\tTEMP2\t$nTEMP2\t$nL_TEMP2\t$nM_TEMP2\t$nH_TEMP2" >> ${DIR}/ANALYSES/post_analysis/$class.denovo.TEMP2.tab
echo -e "$class\tPopTE2\t$nPopTE2\t$nL_PopTE2\t$nM_PopTE2\t$nH_PopTE2" >> ${DIR}/ANALYSES/post_analysis/$class.denovo.PopTE2.tab
echo -e "$class\tTEMP2\t$nTEMP2\t$nL_TEMP2\t$nM_TEMP2\t$nH_TEMP2" > ${DIR}/ANALYSES/post_analysis/$class.denovo.tab
echo -e "$class\tPopTE2\t$nPopTE2\t$nL_PopTE2\t$nM_PopTE2\t$nH_PopTE2" >> ${DIR}/ANALYSES/post_analysis/$class.denovo.tab

#echo -e "$class\t$n"
cd ${DIR}/ANALYSES/post_analysis/
Rscript distribution_denovo_frequency.R $class

cat ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.bed | sort -k7 -r  > ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.sort.bed

grep "TEMP2" ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.sort.bed > ${DIR}/ANALYSES/post_analysis/denovoTE/$class.denovoTE.sort.TEMP2.bed

done


