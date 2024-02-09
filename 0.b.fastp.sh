#!/bin/bash

#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --job-name=bwa_waddington
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH -o logs/fastp_%a.out
#SBATCH -e logs/fastp_%a.err  

module load Miniconda3/4.7.10

source activate anopheles

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/samples_processed.lst))

read1=$(echo "${DIR}/DATA_processed/${sample}/${sample}_1.fq.gz")
read2=$(echo "${DIR}/DATA_processed/${sample}/${sample}_2.fq.gz")

gunzip $read1
gunzip $read2

sed 's/ 1:N:0:\(.*\)$/\#\1\/1/g' ${DIR}/DATA_processed/${sample}/${sample}_1.fq > ${DIR}/DATA_processed/${sample}/${sample}_id_1.fq 
sed 's/ 2:N:0:\(.*\)$/\#\1\/2/g' ${DIR}/DATA_processed/${sample}/${sample}_2.fq > ${DIR}/DATA_processed/${sample}/${sample}_id_2.fq 

gzip ${DIR}/DATA_processed/${sample}/${sample}_id_1.fq 
gzip ${DIR}/DATA_processed/${sample}/${sample}_id_2.fq 

rm -f  ${DIR}/DATA_processed/${sample}/${sample}_id_1.fq 
rm -f  ${DIR}/DATA_processed/${sample}/${sample}_id_2.fq 
rm -f  ${DIR}/DATA_processed/${sample}/${sample}_1.fq 
rm -f  ${DIR}/DATA_processed/${sample}/${sample}_2.fq

echo $read1 $read2 $genome ${sample}

fastp
mkdir ${DIR}/DATA_processed/${sample}/report_fastp
fastp -w 12 -h ${DIR}/DATA_processed/${sample}/report_fastp/${sample}.html \
-j ${DIR}/DATA_processed/${sample}/report_fastp/${sample}.json \
-i ${DIR}/DATA_processed/${sample}/${sample}_id_1.fq.gz \
-I ${DIR}/DATA_processed/${sample}/${sample}_id_2.fq.gz \
-o ${DIR}/DATA_processed/${sample}/${sample}_id_trim_1.fq.gz \
-O ${DIR}/DATA_processed/${sample}/${sample}_id_trim_2.fq.gz \
-l 20 -q 20