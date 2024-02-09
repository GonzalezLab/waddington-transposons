#!/bin/bash

#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --job-name=temp2
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH -o logs/temp_abs_%a.out
#SBATCH -e logs/temp_abs_%a.err  

module load TEMP2/0.1.4-foss-2016b
# module load Miniconda3/4.7.10
# source activate anopheles
# PRESENCE MODULE
# module load BEDTools/2.27.1-foss-2016b
# ABSENCE MODULE
module load BEDTools/2.25.0-foss-2016b
module load BWA/0.7.17-GCC-9.3.0
module load BioPerl/1.7.7-GCCcore-9.3.0 

# sbatch --array=1-12 2b.map-temp2.sh

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/samples_processed.lst))

genome="${DIR}/referenceGenome/dmel-chr-r6.31.fasta"
genomeIndex="${DIR}/referenceGenome/dmel-chr-r6.31.fasta"
repeatmasker="$DIR/ANALYSES/TEMP2/input_TEMP2/dmel-chr-r6.31.fasta.out.bed"
genome2bit="$DIR/ANALYSES/TEMP2/input_TEMP2/dmel-chr-r6.31.2bit"

read1=$(echo "${DIR}/DATA_processed/${sample}/${sample}_id_trim_1.fq.gz")
read2=$(echo "${DIR}/DATA_processed/${sample}/${sample}_id_trim_2.fq.gz")

bam=$(echo ${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.sorted.bam)

#mkdir ${DIR}/ANALYSES/TEMP2/input_TEMP2/insertSize

# picard CollectInsertSizeMetrics \
#          I="${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.sorted.bam " \
#          R="${genome}" \
#          O="${DIR}/ANALYSES/TEMP2/input_TEMP2/insertSize/${sample}.insertSize.txt" \
#          H="${DIR}/ANALYSES/TEMP2/input_TEMP2/insertSize/${sample}.insertSize.pdf"

# mkdir -p ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion
#mkdir ${DIR}/ANALYSES/TEMP2/output/${sample}/absence

rm -rf ${DIR}/ANALYSES/TEMP2/output/${sample}/absence
mkdir -p ${DIR}/ANALYSES/TEMP2/output/${sample}/absence


float=$(cat "${DIR}/ANALYSES/TEMP2/input_TEMP2/insertSize/${sample}.insertSize.txt" | head -n8  | tail -n1 | awk  '{print $6}' )
size=$(echo "($float+0.5)/1" | bc)

# TEMP2 insertion \
#  -l ${read1} \
#  -r ${read2} \
#  -i ${bam} \
#  -g ${genome} \
#  -R ${DIR}/TE_annotation/consensuses_curated_v4.fasta \
#  -t ${repeatmasker} \
#  -o ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion \
#  -p ${sample} \
#  -c 12 \
#  -f ${size} \
#  -m 5 

TEMP2 absence \
-i ${bam} \
-o ${DIR}/ANALYSES/TEMP2/output/${sample}/absence \
-r ${repeatmasker} \
-t ${genome2bit} \
-f ${size} \
-c 12 \
-x 30
