#!/bin/bash

#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --job-name=temp2
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH -o logs/temp_%a.out
#SBATCH -e logs/temp_%a.err  
#SBATCH --exclude=mr-03-[01-26],mr-03-28,mr-04-00

#module load Miniconda3/4.7.10
#source activate anopheles

# sbatch --array=1-12 2a.map-temp2.sh

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/samples_processed.lst))

# RepeatMasker format
# rmsk2bed < ${DIR}/referenceGenome/RepeatMasker/dmel-chr-r6.31.fasta.out | cut -f1-6  > $DIR/ANALYSES/TEMP2/input_TEMP2/dmel-chr-r6.31.fasta.out.bed

# faToTwoBit ${DIR}/referenceGenome/dmel-chr-r6.31.fasta $DIR/ANALYSES/TEMP2/input_TEMP2/dmel-chr-r6.31.2bit

#mkdir -p ${DIR}/referenceGenome/index
#cd ${DIR}/referenceGenome/index
#bwa index ${DIR}/referenceGenome/dmel-chr-r6.31.fasta

genome="${DIR}/referenceGenome/dmel-chr-r6.31.fasta"
genomeIndex="${DIR}/referenceGenome/dmel-chr-r6.31.fasta"
repeatmasker="$DIR/ANALYSES/TEMP2/input_TEMP2/dmel-chr-r6.31.fasta.out.bed"
genome2bit="$DIR/ANALYSES/TEMP2/input_TEMP2/dmel-chr-r6.31.2bit"

read1=$(echo "${DIR}/DATA_processed/${sample}/${sample}_id_trim_1.fq.gz")
read2=$(echo "${DIR}/DATA_processed/${sample}/${sample}_id_trim_2.fq.gz")

# header=$(zcat "$read1" | head -n 1)
# id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
# sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
# echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"

# #Map

# mkdir -p ${DIR}/ANALYSES/TEMP2/map/${sample}/map/

# bwa mem \
# -t 12 \
# -v 3 \
# -Y \
# -T 20 \
# -R $(echo "@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA") \
# $genomeIndex \
# $read1 $read2 > ${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.sam 

# # #Compress
# samtools view -S -b ${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.sam  > ${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.bam 
# rm -f ${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.sam 

# #Sorting
# samtools sort -O bam -o ${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.sorted.bam  ${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.bam 
# rm -f ${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.bam 

# #Index 
# samtools index ${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.sorted.bam  


bam=$(echo ${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.sorted.bam)

# mkdir ${DIR}/ANALYSES/TEMP2/input_TEMP2/insertSize

# picard CollectInsertSizeMetrics \
#           I="${DIR}/ANALYSES/TEMP2/map/${sample}/map/${sample}.sorted.bam" \
#           R="${genome}" \
#           O="${DIR}/ANALYSES/TEMP2/input_TEMP2/insertSize/${sample}.insertSize.txt" \
#           H="${DIR}/ANALYSES/TEMP2/input_TEMP2/insertSize/${sample}.insertSize.pdf"

rm -rf ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion
mkdir -p ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion

float=$(cat "${DIR}/ANALYSES/TEMP2/input_TEMP2/insertSize/${sample}.insertSize.txt" | head -n8  | tail -n1 | awk  '{print $6}' )
size=$(echo "($float+0.5)/1" | bc)

module load TEMP2/0.1.4-foss-2016b
# PRESENCE MODULE
# module load BWA/0.7.17-GCC-9.3.0
module load BEDTools/2.27.1-foss-2016b
module load SAMtools/1.12-GCCcore-8.2.0

TEMP2 insertion \
 -l ${read1} \
 -r ${read2} \
 -i ${bam} \
 -g ${genome} \
 -R ${DIR}/TE_annotation/consensuses_curated_v4.fasta \
 -t ${repeatmasker} \
 -o ${DIR}/ANALYSES/TEMP2/output/${sample}/insertion \
 -p ${sample} \
 -c 8 \
 -f ${size} \
 -m 5 

