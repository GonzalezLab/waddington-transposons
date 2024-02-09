#!/bin/bash

#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --job-name=temp
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH -o logs/processing-temp2_%a.out
#SBATCH -e logs/processing-temp2_%a.err  

module load Miniconda3/4.7.10

source activate anopheles

# Run

# sbatch --array=1-12 3b.process-temp2.sh

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/samples_processed.lst))

rm -rf ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion/tmpTEMP2

# heterochromatin boundaries
heterochromatin="${DIR}/referenceGenome/heterochromatin.bed"

# we keep insertions inside
awk -v sample="$sample" -v OFS='\t' '{print $1, $2, $3, $4, sample, $6 }' ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion/${sample}.insertion.bed | tail -n +2  | sort -k1,1 -k2,2n | sed "s/:[^\t]*//"  > ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion/${sample}.teinsertion.bed

bedtools intersect \
-a ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion/${sample}.teinsertion.bed \
-b ${heterochromatin} > \
${DIR}/ANALYSES/TEMP2/output/${sample}/insertion/${sample}.teinsertion.eu.bed

# recover reference insertions
awk -v sample="$sample" -v OFS='\t' '{print $1, $2, $3, $4, sample, "." }' ${DIR}/ANALYSES/TEMP2/output/${sample}/absence/${sample}.absence.refined.bp.summary | tail -n +2 | sort -k1,1 -k2,2n  > ${DIR}/ANALYSES/TEMP2/output/${sample}/absence/${sample}.absence.bed

bedtools intersect \
-a ${DIR}/ANALYSES/TEMP2/output/${sample}/absence/${sample}.absence.bed \
-b ${heterochromatin} > \
${DIR}/ANALYSES/TEMP2/output/${sample}/absence/${sample}.absence.eu.bed

# all reference insertions in eu
bedtools intersect \
-a $DIR/ANALYSES/TEMP2/input_TEMP2/dmel-chr-r6.31.fasta.out.bed \
-b ${heterochromatin} > \
$DIR/ANALYSES/TEMP2/input_TEMP2/dmel-chr-r6.31.fasta.eu.out.bed

# keep insertions that are in reference not in absence: they are reference
bedtools intersect \
-v \
-a $DIR/ANALYSES/TEMP2/input_TEMP2/dmel-chr-r6.31.fasta.eu.out.bed  \
-b ${DIR}/ANALYSES/TEMP2/output/${sample}/absence/${sample}.absence.eu.bed | awk -v sample="$sample" -v OFS='\t' '{ print $1, $2, $3, $4, sample, "." }' >\
${DIR}/ANALYSES/TEMP2/output/${sample}/absence/${sample}.refinsertion.eu.bed


# we expect that a few insertions overlap
bedtools intersect \
-a ${DIR}/ANALYSES/TEMP2/output/${sample}/absence/${sample}.refinsertion.eu.bed \
-b ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion/${sample}.teinsertion.eu.bed

# we join the two insertions: de novo and reference
cat ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion/${sample}.teinsertion.eu.bed \
${DIR}/ANALYSES/TEMP2/output/${sample}/absence/${sample}.refinsertion.eu.bed | sort -k1,1 -k2,2n > ${DIR}/ANALYSES/TEMP2/output/${sample}/${sample}.teinsertion.eu.bed

# add family
> ${DIR}/ANALYSES/TEMP2/output/${sample}/${sample}.teinsertion.family.eu.bed

while read TE_insertion
do 
TE=$(echo "$TE_insertion" | cut -f4) 
TE_family=$(awk -v TE="$TE" ' $1 == TE ' ${DIR}/TE_annotation/TE_library_family.csv | cut -f8 ) # we grep the consensus name in the TE.list file to look for the family
echo -e "$TE_insertion\t$TE_family\tTEMP2" >> ${DIR}/ANALYSES/TEMP2/output/${sample}/${sample}.teinsertion.family.eu.bed 
done < ${DIR}/ANALYSES/TEMP2/output/${sample}/${sample}.teinsertion.eu.bed
