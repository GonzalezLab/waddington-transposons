#!/bin/bash

#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --job-name=popte2
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH -o logs/ppileup_%a.out
#SBATCH -e logs/ppileup_%a.err  

module load Miniconda3/4.7.10

source activate anopheles

module load PoPoolationTE2/v1.10.03

# Run

# sbatch --array=1-12 3.PopTE2-ppileup.sh
DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/samples_processed.lst))

mkdir -p ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/

# ppileup
java -Xmx32g -jar $EBROOTPOPTE/popte2-v1.10.03.jar ppileup \
--bam ${DIR}/ANALYSES/PoPoolationTE2/map/${sample}/map/${sample}.sort.bam \
--map-qual 15 \
--hier $DIR/TE_annotation/TE_hierarchy.txt \
--output ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.ppileup.gz  \
--detailed-log

# identifySignatures
java -Xmx32G -jar $EBROOTPOPTE/popte2-v1.10.03.jar identifySignatures \
--ppileup ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.ppileup.gz \
--mode separate \
--output ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.signatures \
--min-count 2 \
--detailed-log

# frequency
java -Xmx32G -jar $EBROOTPOPTE/popte2-v1.10.03.jar frequency \
--ppileup ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.ppileup.gz \
--signature ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.signatures \
--output ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.freqsigs \
--detailed-log

# pairupSignatures
java -Xmx32G -jar $EBROOTPOPTE/popte2-v1.10.03.jar pairupSignatures \
--signature ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.freqsigs  \
--ref-genome ${DIR}/ANALYSES/PoPoolationTE2/input_PoPoolationTE2/dmel.r6.31.TE.merged.fa \
--hier ${DIR}/TE_annotation/TE_hierarchy.txt \
--min-distance -200 \
--max-distance 300 \
--output ${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions \
--detailed-log \
--output-detail medium

