#!/bin/bash

#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --job-name=common_insertions
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH -o logs/common-insertions_%a.out
#SBATCH -e logs/common-insertions_%a.err  

module load Miniconda3/4.7.10

source activate anopheles

# Run

# sbatch --array=1-12 4.common-insertions.sh

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/samples_processed.lst))

mkdir -p ${DIR}/ANALYSES/common_insertions/${sample}

bedops --everything \
${DIR}/ANALYSES/PoPoolationTE2/output/${sample}/${sample}.teinsertions.family.freq.eu.bed <(awk '{print $0 "\t" "NA"}' ${DIR}/ANALYSES/TEMP2/output/${sample}/${sample}.teinsertion.family.eu.bed) | sort -k1,1 -k2,2n \
> ${DIR}/ANALYSES/common_insertions/${sample}/${sample}.teinsertion.eu.popte2.freq.temp2.bed

mergeDistance=20

bedtools merge -i ${DIR}/ANALYSES/common_insertions/${sample}/${sample}.teinsertion.eu.popte2.freq.temp2.bed \
 -c 4,5,7,8,7,8,9,9 -d ${mergeDistance} -o distinct,distinct,distinct,distinct,count_distinct,count_distinct,distinct,count_distinct -delim ';' > \
${DIR}/ANALYSES/common_insertions/${sample}/${sample}.teinsertion.eu.popte2.freq.temp2.collapse.bed

awk ' $8 == 1 && $9 == 2 && $11 <= 2' ${DIR}/ANALYSES/common_insertions/${sample}/${sample}.teinsertion.eu.popte2.freq.temp2.collapse.bed | cut -f 1-10 > ${DIR}/ANALYSES/common_insertions/${sample}/${sample}.teinsertion.eu.popte2.freq.temp2.reliable.bed 

sed -i "s/;NA//g" ${DIR}/ANALYSES/common_insertions/${sample}/${sample}.teinsertion.eu.popte2.freq.temp2.reliable.bed 


