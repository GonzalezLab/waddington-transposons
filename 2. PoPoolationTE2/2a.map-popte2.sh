#!/bin/bash

#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --job-name=bwa_waddington
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH -o logs/bwa_se2pe_%a.out
#SBATCH -e logs/bwa_se2pe_%a.err  

module load Miniconda3/4.7.10

source activate anopheles

module load PoPoolationTE2/v1.10.03

# sbatch --array=1-12 2.map.sh

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/samples_processed.lst))

genome="$DIR/ANALYSES/PoPoolationTE2/input_PoPoolationTE2/dmel.r6.31.TE.merged.fa"

#  Map reads to the TE-combined-reference
mkdir -p ${DIR}/ANALYSES/PoPoolationTE2/map/${sample}/map

bwa mem \
-t 12 -M \
${genome} \
${DIR}/DATA_processed_fix/${sample}/${sample}_id_trim_fix_1.fq.gz > \
${DIR}/ANALYSES/PoPoolationTE2/map/${sample}/map/${sample}_id_1.sam 

bwa mem \
-t 12 -M \
${genome} \
${DIR}/DATA_processed_fix/${sample}/${sample}_id_trim_fix_2.fq.gz > \
${DIR}/ANALYSES/PoPoolationTE2/map/${sample}/map/${sample}_id_2.sam 

# Restore paired-end information with PoPoolationTE2 se2pe
java -Xmx32g -jar $EBROOTPOPTE/popte2-v1.10.03.jar se2pe \
--fastq1 ${DIR}/DATA_processed_fix/${sample}/${sample}_id_trim_fix_1.fq.gz \
--fastq2 ${DIR}/DATA_processed_fix/${sample}/${sample}_id_trim_fix_2.fq.gz \
--bam1 ${DIR}/ANALYSES/PoPoolationTE2/map/${sample}/map/${sample}_id_1.sam  \
--bam2 ${DIR}/ANALYSES/PoPoolationTE2/map/${sample}/map/${sample}_id_2.sam  \
--sort \
--detailed-log \
--output ${DIR}/ANALYSES/PoPoolationTE2/map/${sample}/map/${sample}.sort.bam


