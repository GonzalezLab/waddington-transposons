#!/bin/bash

#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --job-name=popte2
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH -o logs/processing-popte2_%a.out
#SBATCH -e logs/processing-popte2_%a.err  

module load Miniconda3/4.7.10

source activate anopheles

# Run

# sbatch --array=1-12 3a.process-popte2.sh

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/samples_processed.lst))

# heterochromatin="${DIR}/referenceGenome/heterochromatin.bed"

# awk -v sample="$sample" -v OFS='\t' '{print $2, $3-1, $3, $5, sample, $4 }' ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions | sort -k1,1 -k2,2n > ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions.bed

# bedtools intersect \
# -a ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions.bed \
# -b ${heterochromatin} > \
# ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions.eu.bed

# # add family
# > ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions.family.eu.bed

# while read TE_insertion
# do 
# TE=$(echo "$TE_insertion" | cut -f4) 
# TE_family=$(awk -v TE="$TE" ' $1 == TE ' ${DIR}/TE_annotation/TE_library_family.csv | cut -f8 ) # we grep the consensus name in the TE.list file to look for the family
# echo -e "$TE_insertion\t$TE_family\tPoPoolationTE2" >> ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions.family.eu.bed
# done < ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions.eu.bed

# add frequency
> ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions.family.freq.eu.bed

while read TE_insertion
do 
chr=$(echo "$TE_insertion" | cut -f1) 

TE=$(echo "$TE_insertion" | cut -f4) 
pos=$(echo "$TE_insertion" | cut -f3) 

freq=$(awk -v TE="$TE" -v pos="$pos" -v chr="$chr" ' $2 == chr && $5 == TE && $3 == pos ' ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions | cut -f 9)
echo -e "$TE_insertion\t$freq" >> ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions.family.freq.eu.bed
done < ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions.family.eu.bed

