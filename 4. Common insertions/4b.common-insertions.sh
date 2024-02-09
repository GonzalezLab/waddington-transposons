#!/bin/bash

#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --job-name=common_insertions
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH -o logs/common-insertions.out
#SBATCH -e logs/common-insertions.err  

module load Miniconda3/4.7.10

source activate anopheles

# Run

# sbatch --array=1-12 4b.common-insertions.sh

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

bedops --everything \
${DIR}/ANALYSES/common_insertions/*_*/*.teinsertion.eu.popte2.freq.temp2.reliable.bed | sort -k1,1 -k2,2n | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $5"="$10}' \
> ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.bed

mergeDistance=20

bedtools merge -i ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.bed \
 -c 4,5,6,4,5,6,10 -d ${mergeDistance} -o distinct,distinct,distinct,count_distinct,count_distinct,count_distinct,distinct -delim ';' > \
${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.collapse.bed

awk ' $9 == 1 ' ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.collapse.bed > ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed 

while read sample 
do
mkdir -p ${DIR}/ANALYSES/common_insertions/reliable/$sample

grep -w "$sample" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | cut -f1,2,3,4,5,6,8,10 > ${DIR}/ANALYSES/common_insertions/reliable/$sample/${sample}.teinsertion.eu.popte2.freq.temp2.reliable.tmp.bed

awk -v sample="$sample" '{
split($NF, a, ";");
for (i in a) {
split(a[i], b, "=");
if (b[1] == sample) {
  print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" b[2];
  break;
}
}
}' ${DIR}/ANALYSES/common_insertions/reliable/$sample/${sample}.teinsertion.eu.popte2.freq.temp2.reliable.tmp.bed > ${DIR}/ANALYSES/common_insertions/reliable/$sample/${sample}.teinsertion.eu.popte2.freq.temp2.reliable.bed

done < <(cut -f1  $DIR/samples_processed.lst)

types=$(cut -f1  $DIR/samples_processed.lst | cut -f2 -d'_' | sort -u)

for type in $types
do
grep "Mix_$type" ${DIR}/ANALYSES/common_insertions/reliable/"Mix_$type"/"Mix_$type".teinsertion.eu.popte2.freq.temp2.reliable.bed | grep "Mix2_$type" | cut -f 1,2,3,4,6,8 > ${DIR}/ANALYSES/post_analysis/freq/"Mix_$type".teinsertion.eu.popte2.freq.temp2.reliable.bed 
grep "Mix2_$type" ${DIR}/ANALYSES/common_insertions/reliable/"Mix2_$type"/"Mix2_$type".teinsertion.eu.popte2.freq.temp2.reliable.bed | grep "Mix_$type" | cut -f 1,2,3,4,6,8 > ${DIR}/ANALYSES/post_analysis/freq/"Mix2_$type".teinsertion.eu.popte2.freq.temp2.reliable.bed 
done


classes=$(cut -f1  $DIR/samples_processed.lst | cut -f1 -d'_' | sort -u)

for type in $types
do
for class in $classes
do
echo -ne "${class}_$type\t"
grep "${class}_$type" ${DIR}/ANALYSES/common_insertions/reliable/"${class}_$type"/"${class}_$type".teinsertion.eu.popte2.freq.temp2.reliable.bed | grep "${class}_$type" | wc -l 

done
done

wc -l ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed

awk '$8 == 1' ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | wc -l
awk '$8 == 12' ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | wc -l

awk '$8 > 1 && $8 < 12' ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | wc -l

cat ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | cut -f5 | grep -E "(Mix2_.*){4}" | wc -l

cat ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | cut -f5 | grep -E "(Mix_.*){4}" | wc -l

cat ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | cut -f5 | grep -E "(D907_.*){4}" | wc -l

for type in $types
do
for class in $classes
do
echo -n "${class}_$type D907_$type "
grep "${class}_$type" ${DIR}/ANALYSES/common_insertions/reliable/"${class}_$type"/"${class}_$type".teinsertion.eu.popte2.freq.temp2.reliable.bed | grep "D907_$type" | wc -l 
echo -n "${class}_$type Mix_$type "
grep "${class}_$type" ${DIR}/ANALYSES/common_insertions/reliable/"${class}_$type"/"${class}_$type".teinsertion.eu.popte2.freq.temp2.reliable.bed | grep "Mix_$type" | wc -l
echo -n "${class}_$type Mix2_$type "
grep "${class}_$type" ${DIR}/ANALYSES/common_insertions/reliable/"${class}_$type"/"${class}_$type".teinsertion.eu.popte2.freq.temp2.reliable.bed | grep "Mix2_$type" | wc -l

done
done

# TEs that are in A, not in P0 nor control nor N
grep "D907_A" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "D907_P0|D907_N|D907_C"  | wc -l

grep "D907_A" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "D907_P0|D907_N|D907_C"  > ${DIR}/ANALYSES/post_analysis/genes/D907_A.bed


# TEs that are in N, not in P0 nor control nor A
grep "D907_N" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "D907_P0|D907_A|D907_C" | wc -l
grep "D907_N" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "D907_P0|D907_A|D907_C" > ${DIR}/ANALYSES/post_analysis/genes/D907_N.bed


# TEs that were in P0 and were lost in A (present in C and N)
grep "D907_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "D907_A" | grep -E "D907_N|D907_C" | wc -l
grep "D907_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "D907_A" | grep -E "D907_N|D907_C" > ${DIR}/ANALYSES/post_analysis/genes/D907_A_lost.bed


# TEs that were in P0 and were lost in N (present in C and A)
grep "D907_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "D907_N" | grep -E "D907_A|D907_C" | wc -l 
grep "D907_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "D907_N" | grep -E "D907_A|D907_C" > ${DIR}/ANALYSES/post_analysis/genes/D907_N_lost.bed


# TEs that are in A, not in P0 nor control nor N
grep "Mix_A" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "Mix_P0|Mix_N|Mix_C" | wc -l

grep "Mix_A" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "Mix_P0|Mix_N|Mix_C" > ${DIR}/ANALYSES/post_analysis/genes/Mix_A.bed


# TEs that are in N, not in P0 nor control nor A
grep "Mix_N" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "Mix_P0|Mix_A|Mix_C" | wc -l

grep "Mix_N" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "Mix_P0|Mix_A|Mix_C" > ${DIR}/ANALYSES/post_analysis/genes/Mix_N.bed

# TEs that were in P0 and were lost in A (present in C and N)
grep "Mix_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "Mix_A" | grep -E "Mix_N|Mix_C"| wc -l
grep "Mix_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "Mix_A" | grep -E "Mix_N|Mix_C" > ${DIR}/ANALYSES/post_analysis/genes/Mix_A_lost.bed


# TEs that were in P0 and were lost in N (present in C and A)
grep "Mix_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "Mix_N" | grep -E "Mix_A|Mix_C"| wc -l

grep "Mix_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "Mix_N" | grep -E "Mix_A|Mix_C" > ${DIR}/ANALYSES/post_analysis/genes/Mix_N_lost.bed


# TEs that are in A, not in P0 nor control nor N
grep "Mix2_A" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "Mix2_P0|Mix2_N|Mix2_C" | wc -l

grep "Mix2_A" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "Mix2_P0|Mix2_N|Mix2_C"  > ${DIR}/ANALYSES/post_analysis/genes/Mix2_A.bed


# TEs that are in N, not in P0 nor control nor A
grep "Mix2_N" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "Mix2_P0|Mix2_A|Mix2_C" | wc -l

grep "Mix2_N" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -vE "Mix2_P0|Mix2_A|Mix2_C" > ${DIR}/ANALYSES/post_analysis/genes/Mix2_N.bed 


# TEs that were in P0 and were lost in A (present in C and N)
grep "Mix2_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "Mix2_A" | grep -E "Mix2_N|Mix2_C"| wc -l

grep "Mix2_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "Mix2_A" | grep -E "Mix2_N|Mix2_C" > ${DIR}/ANALYSES/post_analysis/genes/Mix2_A_lost.bed 


# TEs that were in P0 and were lost in N (present in C and A)
grep "Mix2_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "Mix2_N" | grep -E "Mix2_A|Mix2_C"| wc -l

grep "Mix2_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep -v "Mix2_N" | grep -E "Mix2_A|Mix2_C" > ${DIR}/ANALYSES/post_analysis/genes/Mix2_N_lost.bed 



# # TEs that are in A, not in P0 nor control nor N
# grep "Mix2_A" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep "Mix_A" | grep -vE "Mix2_P0|Mix2_N|Mix2_C|Mix_P0|Mix_N|Mix_P0" 

# # TEs that are in N, not in P0 nor control nor A
# grep "Mix2_N" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep "Mix_N" | grep -vE "Mix2_P0|Mix2_A|Mix2_C|Mix_P0|Mix_A|Mix_P0" 

# # TEs that were in P0 and were lost in A (present in C and N)
# grep "Mix2_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep "Mix_P0" | grep -vE "Mix2_A|Mix_A" | grep -E "Mix2_N|Mix2_C|Mix_N|Mix_C"

# # TEs that were in P0 and were lost in N (present in C and A)
# grep "Mix2_P0" ${DIR}/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | grep "Mix_P0" | grep -vE "Mix2_N|Mix_N" | grep -E "Mix2_A|Mix2_C|Mix_A|Mix_C"

types="A N"
classes=$(cut -f1  $DIR/samples_processed.lst | cut -f1 -d'_' | sort -u)

for class in $classes
do
for type in $types
do
bedtools window -a ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.bed -b ${DIR}/referenceGenome/dmel-all-r6.31.gtf > ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.intersect_1kb.bed

if [[ $class == "Mix2" ]]; then
class2="Mix"
elif [[ $class == "Mix" ]]; then
class2="Mix"
elif [[ $class == "D907" ]]; then
class2="D907"
fi

nTes=$(cat  ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.bed | wc -l )
nGenes=$(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u | wc -l)
wingUR=$(comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_UR_wing_${class2}_${type}.lst | sort -u ) | wc -l )
echo "$class $type - wing UR"
comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_UR_wing_${class2}_${type}.lst | sort -u )


wingDR=$(comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_DR_wing_${class2}_${type}.lst| sort -u ) | wc -l )
echo "$class $type - wing DR"
comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_DR_wing_${class2}_${type}.lst| sort -u ) 


pupalUR=$(comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_UR_pupal_${class2}_${type}.lst| sort -u ) | wc -l )
echo "$class $type - pupal UR"
comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_UR_pupal_${class2}_${type}.lst| sort -u )

pupalDR=$(comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_DR_pupal_${class2}_${type}.lst| sort -u ) | wc -l )
echo "$class $type - pupal DR"
comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_DR_pupal_${class2}_${type}.lst| sort -u ) 

bedtools window -a ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.bed -b ${DIR}/referenceGenome/dmel-all-r6.31.gtf > ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.intersect_1kb.bed
nTesLost=$(cat ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.bed  | wc -l )
nGenesLost=$(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u | wc -l)
wingURlost=$(comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_UR_wing_${class2}_${type}.lst| sort -u ) | wc -l )
echo "$class $type - wing UR  - lost"
comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_UR_wing_${class2}_${type}.lst| sort -u ) 


wingDRlost=$(comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_DR_wing_${class2}_${type}.lst| sort -u ) | wc -l )
echo "$class $type - wing DR - lost"
comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_DR_wing_${class2}_${type}.lst| sort -u )

pupalURlost=$(comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_UR_pupal_${class2}_${type}.lst| sort -u ) | wc -l )
echo "$class $type - pupal UR - lost"
comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_UR_pupal_${class2}_${type}.lst| sort -u )

pupalDRlost=$(comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_DR_pupal_${class2}_${type}.lst| sort -u ) | wc -l )
echo "$class $type - pupal DR  - lost"
comm -12 <(cut -f19 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_lost.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_DR_pupal_${class2}_${type}.lst| sort -u ) 

echo -e "${class}_$type\t$nTes\t$nGenes\t$wingUR\t$wingDR\t$pupalUR\t$pupalDR\t$nTesLost\t$nGenesLost\t$wingURlost\t$wingDRlost\t$pupalURlost\t$pupalDRlost"

done
done


